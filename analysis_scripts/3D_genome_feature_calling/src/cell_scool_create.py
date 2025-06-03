#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import cooler
import pandas as pd
import tabix
from tqdm import tqdm
from cooler import create_scool
from collections import Counter,defaultdict
import pybedtools
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",default=1,help="number of threads to use")
    p.add_argument('--chrom_size',type=str,dest="chrom_size",help="chromosome size file")
    p.add_argument('--blacklist',type=str,dest="blacklist",help="blacklist tabix file indicating the regions that should be excluded from the analysis")
    p.add_argument('--resolution',type=int,dest="resolution",help="resolution of the cooler file, default is 10000 bp")
    p.add_argument('--input_pair_folder',type=str,dest="input_pair_folder",help="input pair folder containing the pair files")
    p.add_argument('--output_prefix',type=str,dest="output_prefix",help="output prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


########################
# utility functions
########################

def is_in_blacklist(black_list_tbi, chrom, pos):
    try:
        for record in black_list_tbi.query(chrom, pos, pos + 1):
            return True
    except:
        return False


def get_chrom_offsets(bins_df):
    chrom_offset = {chrom: bins_df[bins_df['chrom'] == chrom].index[0]
                    for chrom in bins_df['chrom'].cat.categories}
    return chrom_offset


########################
# main function
########################

def process_cells(cell_batch_table, thread, blacklist, chrom_size_path, resolution, output_prefix): # parse the chromosome size file
    chrom_sizes = pd.read_csv(
        chrom_size_path,
        sep='\t',
        index_col=0,
        header=None).squeeze(axis=1)
    bins_df = cooler.util.binnify(chrom_sizes, resolution)
    chrom_offset = get_chrom_offsets(bins_df) 
    # split the cell_batch_table into batches
    chunk_dicts = defaultdict(dict)
    batch_n = max(1, thread - 2)
    for i, (cell, path) in enumerate(cell_batch_table.items()):
        chunk_dicts[i // batch_n][cell] = path
    # prepare output path
    output_path = f"{output_prefix}.scool"
    # output_path exists, remove it
    if os.path.exists(output_path):
        subprocess.run(['rm', '-f', output_path], check=True)
    # 
    if thread == 1:
        for batch_index, chunk_dict in chunk_dicts.items():
            batch_output = f"{output_prefix}.batch_{batch_index}.scool"
            print(f"Processing batch {batch_index}")
            process_chunk(chunk_dict, bins_df, chrom_offset, chrom_size_path, blacklist, resolution, output_path = batch_output)
            cell_pixel_dict = {}
            with pd.HDFStore(batch_output, mode='r') as hdf:
                for cell_id in hdf.keys():
                    cell_id = cell_id[1:]
                    cell_pixel_dict[cell_id] = hdf[cell_id]
            create_scool(output_path,
                         bins=bins_df,
                         cell_name_pixels_dict=cell_pixel_dict,
                         ordered=True,
                         mode='a')
            subprocess.run(['rm', '-f', batch_output], check=True)
    elif thread > 1:
        with ProcessPoolExecutor(max_workers=thread) as executor:
            futures = {}
            for batch_index, chunk_dict in chunk_dicts.items():
                batch_output = f"{output_prefix}.batch_{batch_index}.scool"
                future = executor.submit(
                    process_chunk,
                    chunk_dict,
                    bins_df,
                    chrom_offset,
                    chrom_size_path, 
                    blacklist,
                    resolution,
                    output_path = batch_output)
                futures[future] = batch_output
            for future in as_completed(futures):
                batch_output = futures[future]
                future.result()
                # dump batch result into cool
                cell_pixel_dict = {}
                with pd.HDFStore(batch_output, mode='r') as hdf:
                    for cell_id in hdf.keys():
                        cell_id = cell_id[1:]  # remove '/'
                        cell_pixel_dict[cell_id] = hdf[cell_id]
                create_scool(output_path,
                             bins=bins_df,
                             cell_name_pixels_dict=cell_pixel_dict,
                             ordered=True,
                             mode='a')
                subprocess.run(['rm', '-f', batch_output], check=True)
    else:
        print("Error: thread number should be greater than 0")
    print(f"Output file: {output_path}")


def process_chunk(cell_path_dict, bins_df, chrom_offset, chrom_size_path, blacklist, resolution, output_path = None):
    """
    this function process a chunk of cells' pair files can save the result into a single scool file
    """
    with pd.HDFStore(output_path, mode='w') as hdf:
        # filter read pairs
        for cell_id, pair_file in cell_path_dict.items():
            # save the result
            hdf[cell_id] = process_single_cell(cell_id, pair_file, bins_df, chrom_offset, chrom_size_path, blacklist, resolution)
    

def process_single_cell(cell_id, pair_file, bins_df, chrom_offset, chrom_size_path, blacklist, resolution):
    # filter the read pairs
    contacts = filter_single_read_pair_file(pair_file, chrom_offset, chrom_size_path, min_pos_dist = 1000, blacklist_file = blacklist, remove_duplicate = True)
    # convert the read pairs into bins
    bin1_id= (contacts['pos1'] // resolution) + contacts['chrom1'].map(chrom_offset)
    bin2_id = (contacts['pos2'] // resolution) + contacts['chrom2'].map(chrom_offset)
    counter = Counter()
    for x, y in zip(bin1_id, bin2_id):
        if x > y:
            x, y = y, x
        counter[(x, y)] += 1
    pixel_df = pd.Series(counter).reset_index()
    pixel_df.columns = pd.Index(['bin1_id', 'bin2_id', 'count'])
    pixel_df = pixel_df.sort_values(['bin1_id',
                                         'bin2_id']).reset_index(drop=True)
    return pixel_df


def filter_single_read_pair_file(read_pair_file, chrom_offset, chrom_size_path, min_pos_dist = 1000, blacklist_file = None, remove_duplicate = True):
    """
    Load read pair in format chr6    257496  chr6    2684661
    - Skip read pairs that the distance between the two positions is less than min_pos_dist
    - Skip read pairs that either of the two positions is in the blacklist
    Convert the read pair position into bin position
    - The bin position is calculated by pos // resolution
    """
    # using pandas to load the file
    if '.gz' in read_pair_file:
        contacts = pd.read_csv(read_pair_file, sep = "\t", header = None, compression = 'gzip', names = ['chrom1', 'pos1', 'chrom2', 'pos2'], dtype = {'chrom1': str, 'pos1': int, 'chrom2': str, 'pos2': int}, comment = '#')
    else:
        contacts = pd.read_csv(read_pair_file, sep = "\t", header = None, names = ['chrom1', 'pos1', 'chrom2', 'pos2'], dtype = {'chrom1': str, 'pos1': int, 'chrom2': str, 'pos2': int}, comment = '#')
    N_pair = contacts.shape[0]
    # drop chroms that are not in the chrom_size
    contacts = contacts[contacts['chrom1'].isin(chrom_offset.keys()) & contacts['chrom2'].isin(chrom_offset.keys())]
    # drop duplicate rows
    if remove_duplicate:
        contacts.drop_duplicates(inplace = True)
    # min_pos_dist for intra-chromosomal contacts
    if min_pos_dist > 0:
        contacts = contacts[abs(contacts['pos1'] - contacts['pos2']) >= min_pos_dist | (contacts['chrom1'] != contacts['chrom2'])]
    # remove blacklist regions
    if blacklist_file is not None:
        if '.gz' in blacklist_file:
            blacklist_bed_df = pd.read_csv(blacklist_file, sep='\t', index_col=None, header=None, compression = 'gzip')
        else:
            blacklist_bed_df = pd.read_csv(blacklist_file, sep='\t', index_col=None, header=None)
        blacklist_bed_df = blacklist_bed_df[blacklist_bed_df.iloc[:, 0].isin(chrom_offset)].copy()
        blacklist_bed = pybedtools.BedTool.from_dataframe(blacklist_bed_df).sort(g=chrom_size_path)
        left_bed_df = contacts[['chrom1', 'pos1', 'pos1']].reset_index().iloc[:, [1, 2, 3, 0]]
        right_bed_df = contacts[['chrom2', 'pos2', 'pos2']].reset_index().iloc[:, [1, 2, 3, 0]]
        left_bed_df.columns = ['chrom', 'start', 'end', 'id']
        right_bed_df.columns = ['chrom', 'start', 'end', 'id']
        contact_bed_df = pd.concat([left_bed_df, right_bed_df])
        contact_bed = pybedtools.BedTool.from_dataframe(contact_bed_df).sort(g=chrom_size_path)
        bad_contacts = contact_bed.intersect(blacklist_bed, wa=True, u=True).to_dataframe()
        if bad_contacts.shape[0] > 0:
            bad_index = bad_contacts['name'].unique()
        else:
            # no bad contacts
            bad_index = set() 
        #black_list_tbi = tabix.open(blacklist_file)
        #bad_index = []
        ## check the blacklist by each row
        #for index, row in tqdm(contacts.iterrows()):
        #    chrom1 = row['chrom1']
        #    pos1 = row['pos1']
        #    chrom2 = row['chrom2']
        #    pos2 = row['pos2']
        #    if is_in_blacklist(black_list_tbi, chrom1, pos1) or is_in_blacklist(black_list_tbi, chrom2, pos2):
        #        bad_index.append(index)
        # drop the bad index
        contacts.drop(bad_index, inplace = True)
    # copy into a new dataframe
    final_contacts = contacts[['chrom1', 'pos1', 'chrom2', 'pos2']].copy()
    N_pair_filtered = final_contacts.shape[0]
    # log how many read pairs are filtered
    print(f"{read_pair_file}: Filtered read pairs: {N_pair - N_pair_filtered} ({(N_pair - N_pair_filtered) / N_pair * 100:.2f}%)") 
    # return the result
    return final_contacts
 

def main():
    global args
    args = parse_arg() 
    # list the pair files in the input folder
    cell_batch_table = {}
    for filename in os.listdir(args.input_pair_folder): 
        if '.pair' not in filename:
            continue
        pair_file = os.path.join(args.input_pair_folder, filename)
        cell_id = filename.split(".")[0]
        cell_id = cell_id.replace(",", "_comma_") # fix the comma in the cell id
        cell_id = cell_id.replace("-", "_dash_")
        cell_batch_table[cell_id] = pair_file 
    print(f"Total {len(cell_batch_table)} files to process in the input folder {args.input_pair_folder}")
    # process cells
    process_cells(cell_batch_table, args.thread, args.blacklist, args.chrom_size, args.resolution, args.output_prefix)
    print("All done!")

    
if __name__=="__main__":
    main()

