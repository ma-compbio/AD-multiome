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
    p.add_argument('--ctcf_peak',type=str,dest="ctcf_peak",help="The CTCF peak annotation file")
    p.add_argument("--output_prefix", dest="output_prefix", type=str, default=None, help="The output file prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def build_MUSIC_bin():
    """
    To derive the contact frequency versus genomic distances heatmap for all frontal cortex cells we first generated 149 genomic bins spanning from 5,000 bp to 150 megabase pairs, with the nth bin spanning the genomic region 
    """
    bin_list = []
    for nn in range(151)[1:]:
        start = 2 ** (np.log2(5000) + (nn-1) * (np.log2(150000000) - np.log2(5000)) / 150)
        end = 2 **(np.log2(5000) + nn * (np.log2(150000000) - np.log2(5000)) / 150)
        bin_list.append((int(start), int(end)))
    # function to determine the bin index given the distance
    def find_bin_index(distance):
        for idx, (start, end) in enumerate(bin_list):
            if distance >= start and distance < end:
                return idx
    def find_bin_index_fast(distance):
        return int((np.log2(distance) - np.log2(5000)) / (np.log2(150000000) - np.log2(5000)) * 150)
    return bin_list, find_bin_index_fast


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


def parse_ctcf_peak(filename, bin_size):
    table = {}
    with open(filename, 'r') as fin:
        for line in tqdm(fin):
            row = line.strip().split()
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            if table.get(chrom) is None:
                table[chrom] = {}
            bin_idx = start // bin_size
            bin_idx_end = end // bin_size
            for idx in range(bin_idx, bin_idx_end):
                table[chrom][idx] = True
    return table


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
    # parse the CTCF peak annotation to get the bin that overlapping with CTCF peaks
    CTCF_bin_size = 10000
    ctcf_anno = parse_ctcf_peak(args.ctcf_peak, CTCF_bin_size)
    print("CTCF peak annotation has been loaded")
    print(f"CTCF bin size is {CTCF_bin_size}")
    total_bin_w_motif = 0
    for chrom in ctcf_anno:
        total_bin_w_motif += len(ctcf_anno[chrom])
    print("Total number of bins with CTCF motif: ", total_bin_w_motif)
    # determine the distance bin
    bin_list_MUSIC, find_bin_index_MUSIC = build_MUSIC_bin()
    bin_list_Mattew, find_bin_index_Mattew = build_Mattew_bin()
    #print(bin_list_MUSIC)
    #print(bin_list_Mattew)
    # read the files and count the number of reads
    out_folder = os.path.dirname(args.output_prefix)
    print(out_folder)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    fout_cell_distance = open(args.output_prefix + ".pair_distance.ctcf.txt", 'w')
    print("Start to process the files")
    print("cell_id\tsample_id\tdistance_group\tbin_idx\tbin_start\tbin_end\tctcf_group\tcount\tcumu_count", file = fout_cell_distance)
    for filename in tqdm(file_list):
        cell_id = filename.split("/")[-1].replace('.pair.gz', '').replace(".pair", "")  
        result = {}
        for ctcf_group in ['ctcf-ctcf', 'ctcf-none', 'none-none']:
            result[ctcf_group] = {'MUSIC': np.zeros(151), 'Mattew': np.zeros(144), 'read_count': 0, 'inter_chrom':0}
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
            # skip chrM
            if chrom_1 == 'chrM' or chrom_2 == 'chrM':
                continue
            result[ctcf_group]['read_count'] += 1
            bin_1 = pos_1 // CTCF_bin_size
            bin_2 = pos_2 // CTCF_bin_size
            # overlap with CTCF peak
            overlap_ctcf_1 = ctcf_anno[chrom_1].get(bin_1, None) 
            overlap_ctcf_2 = ctcf_anno[chrom_2].get(bin_2, None)
            if overlap_ctcf_1 is not None and overlap_ctcf_2 is not None:
                ctcf_group = 'ctcf-ctcf'
            elif overlap_ctcf_1 is not None or overlap_ctcf_2 is not None:
                ctcf_group = 'ctcf-none'
            else:
                ctcf_group = 'none-none'
            # check if the two reads are in the same chromosome
            if chrom_1 != chrom_2:
                result[ctcf_group]['inter_chrom'] += 1
            else:
                distance = abs(pos_2 - pos_1) 
                if distance > 1000:
                    bin_idx_Mattew = find_bin_index_Mattew(distance)
                    result[ctcf_group]['Mattew'][bin_idx_Mattew] += 1
                    if distance > 5000 and distance < 150000000:
                        bin_idx_MUSIC = find_bin_index_MUSIC(distance)
                        result[ctcf_group]['MUSIC'][bin_idx_MUSIC] += 1
                    elif distance >= 150000000:
                        result[ctcf_group]['MUSIC'][150] += 1
        fin.close()
        # output the result
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_MUSIC):
            start = bin_pair[0]
            end = bin_pair[1]
            count_ctcf2ctcf  = result['ctcf-ctcf']['MUSIC'][idx]
            count_ctcf2none = result['ctcf-none']['MUSIC'][idx]
            count_none2none = result['none-none']['MUSIC'][idx]
            count = count_ctcf2ctcf + count_ctcf2none + count_none2none
            cumu_count += count
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tCTCF-CTCF\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_ctcf2ctcf, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tCTCF-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_ctcf2none, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_none2none, cumu_count), file = fout_cell_distance)
        count_ctcf2ctcf  = result['ctcf-ctcf']['MUSIC'][150]
        count_ctcf2none = result['ctcf-none']['MUSIC'][150]
        count_none2none = result['none-none']['MUSIC'][150] 
        count = count_ctcf2ctcf + count_ctcf2none + count_none2none
        cumu_count += count
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tCTCF-CTCF\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2ctcf, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tCTCF-Other\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)
        count_ctcf2ctcf  = result['ctcf-ctcf']['inter_chrom']
        count_ctcf2none = result['ctcf-none']['inter_chrom']
        count_none2none = result['none-none']['inter_chrom'] 
        count = count_ctcf2ctcf + count_ctcf2none + count_none2none
        cumu_count += count        
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tCTCF-CTCF\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2ctcf, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tCTCF-Other\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_Mattew):
            start = bin_pair[0]
            end = bin_pair[1]
            count_ctcf2ctcf  = result['ctcf-ctcf']['Mattew'][idx]
            count_ctcf2none = result['ctcf-none']['Mattew'][idx]
            count_none2none = result['none-none']['Mattew'][idx]
            count = count_ctcf2ctcf + count_ctcf2none + count_none2none
            cumu_count += count
            print("%s\t%s\tMattew\t%d\t%d\t%d\tCTCF-CTCF\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_ctcf2ctcf, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMattew\t%d\t%d\t%d\tCTCF-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_ctcf2none, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMattew\t%d\t%d\t%d\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_none2none, cumu_count), file = fout_cell_distance)
        count_ctcf2ctcf  = result['ctcf-ctcf']['inter_chrom']
        count_ctcf2none = result['ctcf-none']['inter_chrom']
        count_none2none = result['none-none']['inter_chrom'] 
        count = count_ctcf2ctcf + count_ctcf2none + count_none2none
        cumu_count += count        
        print("%s\t%s\tMattew\tinter\t-1\t-1\tCTCF-CTCF\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2ctcf, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMattew\tinter\t-1\t-1\tCTCF-Other\t%d\t%d" % (cell_id, args.patient_id, count_ctcf2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMattew\tinter\t-1\t-1\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)

    

if __name__=="__main__":
    main()

