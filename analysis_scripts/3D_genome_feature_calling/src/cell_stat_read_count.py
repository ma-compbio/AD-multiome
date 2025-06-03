#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
from tqdm import tqdm
import gzip
import numpy as np


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument("--folder", dest="folder", type=str, default=None, help="The folder containing the cell read pairs files")
    p.add_argument("--read_count_cutoff",dest="read_count_cutoff",type=int,default=30000,help="The read count cutoff to filter the cells")
    p.add_argument("--folder_filtered", dest="folder_filtered", type=str, default=None, help="The folder containing the cell read pairs files after the filtering by --read_count_cutoff")
    p.add_argument("--patient_id", dest="patient_id", type=str, default=None, help="The patient id")
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


def main():
    global args
    args = parse_arg()
    # list all the files in the folder
    file_list = []
    for filename in os.listdir(args.folder):
        if filename.endswith(".pair") or filename.endswith(".pair.gz"):
            file_list.append(os.path.join(args.folder, filename))
    print("Total number of files: ", len(file_list))
    file_list = sorted(file_list)
    # if read count cutoff set then filter the files
    if args.read_count_cutoff is not None:
        print("Filter the files by read count cutoff: ", args.read_count_cutoff)
        # check the folder
        if not os.path.exists(args.folder_filtered):
            os.makedirs(args.folder_filtered)
        print("The folder to store the filtered files: ", args.folder_filtered)
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
    fout_cell_count = open(args.output_prefix + ".cell_count.txt", 'w')
    fout_cell_distance = open(args.output_prefix + ".pair_distance.txt", 'w')
    print("Start to process the files")
    print("cell_id\tsample_id\tread_pair_count", file = fout_cell_count)
    print("cell_id\tsample_id\tdistance_group\tbin_idx\tbin_start\tbin_end\tcount\tcumu_count", file = fout_cell_distance)
    for filename in tqdm(file_list):
        cell_id = filename.split("/")[-1].replace('.pair.gz', '').replace(".pair", "")  
        result = {'MUSIC': np.zeros(151), 'Mattew': np.zeros(144), 'read_count': 0, 'inter_chrom':0}
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
            result['read_count'] += 1
            # check if the two reads are in the same chromosome
            if chrom_1 != chrom_2:
                result['inter_chrom'] += 1
            else:
                distance = abs(pos_2 - pos_1) 
                if distance > 1000:
                    bin_idx_Mattew = find_bin_index_Mattew(distance)
                    result['Mattew'][bin_idx_Mattew] += 1
                    if distance > 5000 and distance < 150000000:
                        bin_idx_MUSIC = find_bin_index_MUSIC(distance)
                        result['MUSIC'][bin_idx_MUSIC] += 1
                    elif distance >= 150000000:
                        result['MUSIC'][150] += 1
        fin.close()
        # output the result
        print("%s\t%s\t%d" % (cell_id, args.patient_id, result['read_count']), file = fout_cell_count)
        # filter the cell by read count
        if args.read_count_cutoff is not None and result['read_count'] < args.read_count_cutoff:
            continue
        # copy the file to the folder
        if args.read_count_cutoff is not None:
            new_filename = os.path.join(args.folder_filtered, cell_id + ".pair")
            os.system("cp %s %s" % (filename, new_filename))
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_MUSIC):
            start = bin_pair[0]
            end = bin_pair[1]
            count  = result['MUSIC'][idx]
            cumu_count += count
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\t%d\t%d" % (cell_id, args.patient_id, result['MUSIC'][150], result['MUSIC'][150] + cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\t%d\t%d" % (cell_id, args.patient_id, result['inter_chrom'], result['inter_chrom'] + cumu_count), file = fout_cell_distance)
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_Mattew):
            start = bin_pair[0]
            end = bin_pair[1]
            count = result['Mattew'][idx]
            cumu_count += count
            print("%s\t%s\tMattew\t%d\t%d\t%d\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMattew\tinter\t-1\t-1\t%d\t%d" % (cell_id, args.patient_id, result['inter_chrom'], result['inter_chrom'] + cumu_count), file = fout_cell_distance)


if __name__=="__main__":
    main()

