#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 24 Apr 2025 06:07:33 PM

import os,sys,argparse
import pandas as pd
from tqdm import tqdm
import numpy as np
import gzip


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument("--folder", dest="folder", type=str, default=None, help="The folder containing the cell read pairs files")
    p.add_argument("--patient_id", dest="patient_id", type=str, default=None, help="The patient id")
    p.add_argument('--cell_white_list',type=str,dest="cell_white_list",help="A file contains the the white list of cells by cell id")
    p.add_argument('--cell_white_list_format',type=str,dest="cell_white_list_format",default="csv",choices=['cell_label','cell_id'],help="The format of the cell white list file")
    p.add_argument("--output_prefix", dest="output_prefix", type=str, default=None, help="The output file prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def build_Mattew_bin():
    """
    To count the number of cis (intra-chromosomal) contacts in each cell and bulk Hi-C data, we divided the contacts into 143 logarithmic bins, the first of which was for contacts that were separated by less than 1 kb. Each subsequent bin covered an exponent step of 0.125, using base 2.
    """
    bin_list = []
    for nn in range(144):
        start = 1000 * 2 ** (nn * 0.125)
        end = 1000 * 2 ** ((nn+1) * 0.125)
        bin_list.append((int(start), int(end)))
    # function to determine the bin index given the distance
    def find_bin_index(distance):
        for idx, (start, end) in enumerate(bin_list):
            if distance >= start and distance < end:
                return idx
    def find_bin_index_fast(distance):
        return int(np.log2(distance/1000) / 0.125)
    return bin_list, find_bin_index_fast

def build_Mattem_bin_simple():
    bin_list = [(1000, 128000), (128000, 4096000), (4096000, 65536000), (65536000, 262144000)]
    def find_bin_index(distance):
        for idx, (start, end) in enumerate(bin_list):
            if distance >= start and distance < end:
                return idx
    return bin_list, find_bin_index


def parse_white_list(filename, format = 'cell_label'):
    """
    "","Idents.obj."
    "hPFC-lib-1_A1,A7","Micro"
    "hPFC-lib-1_A1,C5","OPC"
    """
    if format == 'cell_label':
        cell_id2group = pd.read_csv(filename, sep = ',', header = 0)
        cell_id2group.columns = ['cell_id', 'cell_type'] 
        # convert to dictionary with cell_id as key and cell_type as value
        cell_id2group = dict(zip(cell_id2group['cell_id'], cell_id2group['cell_type']))
    elif format == 'cell_id':
        cell_id2group = {}
        with open(filename, 'r') as fin:
            for line in fin:
                if line.strip() == '':
                    continuee
                row = line.strip()
                cell_id = row[0].replace('"', '')
                cell_id2group[cell_id] = True
    else:
        print("Unknown format for cell white list file")
        exit(1)
    return cell_id2group



def main():
    global args
    args = parse_arg()
    # get cell list 
    # get the white list of cells
    if args.cell_white_list is not None:
        cell_white_list = parse_white_list(args.cell_white_list, args.cell_white_list_format)
        print("Total number of selected cells: ", len(cell_white_list))
    else:
        cell_white_list = None
        print("No cell white list is provided, use all the cells in the folder")
    # list all the files in the folder
    file_list = []
    for filename in os.listdir(args.folder):
        if filename.endswith(".pair") or filename.endswith(".pair.gz"):
            cell_id = os.path.basename(filename).replace('.pair.gz', '').replace('.pair', '')
            if cell_white_list is not None:
                if cell_id in cell_white_list:
                    file_list.append(os.path.join(args.folder, filename))
            else:
                file_list.append(os.path.join(args.folder, filename))
    print("Total number of filtered files: ", len(file_list))
    file_list = sorted(file_list)
    # determine the distance bin
    #bin_list_Mattew, find_bin_index_Mattew = build_Mattew_bin()
    bin_list_Mattew, find_bin_index_Mattew = build_Mattem_bin_simple()
    #print(bin_list_MUSIC)
    #print(bin_list_Mattew)
    # read the files and count the number of reads
    out_folder = os.path.dirname(args.output_prefix)
    print(out_folder)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    fout_cell_distance = open(args.output_prefix + ".by_chrom.txt", 'w')
    print("Start to process the files")
    print("cell_id\tsample_id\tdistance_group\tbin_idx\tbin_start\tbin_end\tchrom\tcount\tcumu_count", file = fout_cell_distance)
    chrom_white_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    for filename in tqdm(file_list):
        cell_id = filename.split("/")[-1].replace('.pair.gz', '').replace(".pair", "")  
        result = {}
        for chrom in chrom_white_list:
            result[chrom] = {'Mattew': np.zeros(4), 'read_count': 0, 'inter_chrom':0}
        if '.gz' in filename: 
            fin = gzip.open(filename, 'rt')
        else:
            fin = open(filename, 'r')
        for line in fin:
            row = line.strip().split()
            chrom_1 = row[0]
            pos_1 = int(row[1])
            chrom_2 = row[2]
            pos_2 = int(row[3])
            # skip the line if the chromosome is not in the white list
            if chrom_1 not in chrom_white_list or chrom_2 not in chrom_white_list:
                continue
            result[chrom_1]['read_count'] += 1
            result[chrom_2]['read_count'] += 1
            # check if the two reads are in the same chromosome
            if chrom_1 != chrom_2:
                result[chrom_1]['inter_chrom'] += 1
                result[chrom_2]['inter_chrom'] += 1 
            else:
                distance = abs(pos_2 - pos_1) 
                if distance > 1000:
                    bin_idx_Mattew = find_bin_index_Mattew(distance)
                    result[chrom_1]['Mattew'][bin_idx_Mattew] += 1
        fin.close()
        # output the result
        for chrom in chrom_white_list:
            cumu_count = 0
            for idx, bin_pair in enumerate(bin_list_Mattew):
                start = bin_pair[0]
                end = bin_pair[1]
                count = result[chrom]['Mattew'][idx]
                cumu_count += count
                print("%s\t%s\tMattew\t%d\t%d\t%d\t%s\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, chrom, count, cumu_count), file = fout_cell_distance)
            count = result[chrom]['inter_chrom'] 
            cumu_count += count        
            print("%s\t%s\tMattew\tinter\t-1\t-1\tinter\t%d\t%d" % (cell_id, args.patient_id, count, cumu_count), file = fout_cell_distance)

    

if __name__=="__main__":
    main()

