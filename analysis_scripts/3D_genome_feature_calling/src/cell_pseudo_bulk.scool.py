#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import pandas as pd
from tqdm import tqdm
import cooler
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",help="number of threads to use",default=8)
    p.add_argument('--scool_file',type=str,dest="scool_file",help="scool file containing the raw read pair files for each cell")
    p.add_argument('--cell_anno',type=str,dest="cell_anno",help="cell annotation file")
    p.add_argument('--genome_size',type=str,dest="genome_size",help="chromosome size file")
    p.add_argument('--output_folder',type=str,dest="output_folder",help="output folder")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def report_cell_count(cell_anno, anno_col, out_folder):
    # group by projid and anno_col to count the number of cellso
    cell_count = cell_anno.groupby(['projid', anno_col]).size().reset_index(name='count')
    # write the cell count to a file
    out_file = os.path.join(out_folder, 'cell_count.txt')
    if os.path.exists(out_file):
        print("Warning: file exists: %s" % out_file)
        return True
    else:
        with open(out_file, 'w') as fout:
            print("projid\t%s\tcount" % anno_col, file = fout)
            for index, row in cell_count.iterrows():
                fout.write("{}\t{}\t{}\n".format(row['projid'], row[anno_col], row['count']))


def main():
    global args
    args = parse_arg()
    # load cell annotation file
    cell_anno = pd.read_csv(args.cell_anno, sep="\t", header=0) 
    # add major cell type to the cell annotation file
    if 'major' not in cell_anno.columns:
        # split the sub column and take the first part as major cell type
        cell_anno['major'] = cell_anno['sub'].apply(lambda x: x.split(" ")[0])
    # fix the space in the sub cell type
    cell_anno['sub'] = cell_anno['sub'].apply(lambda x: x.replace(" ", "_").replace("/", "_"))
    # report the summary of the cell type into a file
    out_major_folder = os.path.join(args.output_folder, 'major')
    out_sub_folder = os.path.join(args.output_folder, 'sub')
    for folder in [out_major_folder, out_sub_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
    # count the number of cells in each major cell type
    report_cell_count(cell_anno, 'major', out_major_folder)
    report_cell_count(cell_anno, 'sub', out_sub_folder)
    # for each cell type in one patient, generate the pseudo bulk data
    processing_job = {'major':{}, 'sub': {}}
    # list the single-cell pairs in the scool file
    single_cell_list = cooler.fileops.list_scool_cells(args.scool_file)
    print("Total %d cells in the scool file" % len(single_cell_list))
    for file_path in single_cell_list:
        file_path = args.scool_file + "::" + file_path
        cell_id = file_path.split('/')[-1]
        if '_comma_' in cell_id:
            cell_id = cell_id.replace('_comma_', ',')
        if '_dash_' in cell_id:
            cell_id = cell_id.replace('_dash_', '-')
        # skip the cell that is not in the cell annotation file
        if cell_id not in cell_anno['cell_id'].values:
            continue
        # get patient id from cell annotation file
        patient_id = cell_anno[cell_anno['cell_id'] == cell_id]['projid'].values[0]
        disease_group = cell_anno[cell_anno['cell_id'] == cell_id]['disease_group'].values[0]
        key = "%s_patient_%s" % (disease_group, patient_id)
        # get the cell type of the cell
        cell_type_major = cell_anno[cell_anno['cell_id'] == cell_id]['major'].values[0]
        cell_type_sub = cell_anno[cell_anno['cell_id'] == cell_id]['sub'].values[0]
        if key not in processing_job['major']:
            processing_job['major'][key] = {}
        if key not in processing_job['sub']:
            processing_job['sub'][key] = {}
        if cell_type_major not in processing_job['major'][key]:
            processing_job['major'][key][cell_type_major] = []
        if cell_type_sub not in processing_job['sub'][key]:
            processing_job['sub'][key][cell_type_sub] = []
        # add the pair file to the processing job
        processing_job['major'][key][cell_type_major].append(file_path)
        processing_job['sub'][key][cell_type_sub].append(file_path)
    # do the jobs parallel
    pool = ProcessPoolExecutor(max_workers = args.thread)
    p_list = []
    for cell_label in ['major', 'sub']:
        # get the output folder
        for key in processing_job[cell_label]:
            # major -> AD_patient_patient_id
            out_folder = os.path.join(args.output_folder, cell_label, key)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            for cell_type in processing_job[cell_label][key]:
                # sub -> CD4_T_cell
                # get the chromosome size file
                chrom_size_file = args.genome_size
                # get the list of pairs
                pair_file_list = processing_job[cell_label][key][cell_type]
                p = pool.submit(process_single_cell_type, pair_file_list, cell_type, chrom_size_file, out_folder)
    # do the job
    for p in tqdm(as_completed(p_list), total = len(p_list)):
        p.result()
    # wait and close
    pool.shutdown(wait = True)
    # report to output 
    print("All finished")


def process_single_cell_type(pair_file_list, cell_type, chrom_size_file, out_folder):
    """
    parse the pair file list, remove pairs that in the blacklist, and generate the pseudo bulk data
    """
    out_cool_file = os.path.join(out_folder, cell_type + ".cool")
    print("Processing %s merge %d files into %s" % (cell_type, len(pair_file_list), out_cool_file)) 
    cooler.merge_coolers(out_cool_file, pair_file_list, mergebuf = 20000000)
    print("Finished %s" % cell_type)

    
if __name__=="__main__":
    main()

