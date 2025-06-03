#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 06 Feb 2025 05:46:12 PM

import os,sys,argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",help="number of threads",default=16)
    p.add_argument('--cool_folder',type=str,dest="cool_folder",help="folder containing saddle diginited file")
    p.add_argument('--cell_anno',type=str,dest="cell_anno",help="cell annotation file")
    p.add_argument('--read_pair_folder',type=str,dest="read_pair_folder",help="raw read pair folder")
    p.add_argument('--output',type=str,dest="output",help="output file")
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
    return bin_list


# function to determine the bin index given the distance
def find_bin_index_Mattew(bin_list, distance):
    for idx, (start, end) in enumerate(bin_list):
        if distance >= start and distance < end:
            return idx


def main():
    global args
    args = parse_arg()
    # load cell annotation
    cell_anno = pd.read_csv(args.cell_anno, sep="\t", header=0)
    cell_type2cell_id = dict()
    for index, row in cell_anno.iterrows():
        cell_type = row.loc["sub"]
        cell_type = cell_type.replace(" ", "_").replace("/", "_")
        if cell_type not in cell_type2cell_id:
            cell_type2cell_id[cell_type] = {'AD': dict(), 'CT': dict()}
        cell_id = row.loc["cell_id"]
        disease_group = row.loc["disease_group"]
        patient_id = str(row.loc["projid"])
        cell_type2cell_id[cell_type][disease_group][(patient_id, cell_id)] = row
    print(f"{len(cell_type2cell_id)} cell types loaded")
    for cell_type in sorted(cell_type2cell_id.keys()):
        for disease_group in sorted(cell_type2cell_id[cell_type].keys()):
            print(f"{disease_group} {cell_type}: {len(cell_type2cell_id[cell_type][disease_group])} cells")
    # build bin list
    bin_list_Mattew = build_Mattew_bin()
    # process the cell types
    with open(args.output, 'w') as fout:
        print("disease_group\tcell_type\tbin_idx\tbin_start\tbin_end\tscore_1\tscore_2\tcount", file = fout)
        if args.thread > 1:
            pool = ProcessPoolExecutor(max_workers = args.thread)
            futures = {}
        for cell_type in sorted(cell_type2cell_id.keys()):
            cell_list_AD = sorted(cell_type2cell_id[cell_type]['AD'].keys())
            cell_list_CT = sorted(cell_type2cell_id[cell_type]['CT'].keys())
            # get the saddle digitized results
            saddle_digitized_file_AD = os.path.join(args.cool_folder, 'AD', f"{cell_type}.saddle.100kb.digitized.tsv")
            saddle_digitized_file_CT = os.path.join(args.cool_folder, 'CT', f"{cell_type}.saddle.100kb.digitized.tsv")
            if not os.path.exists(saddle_digitized_file_AD):
                print(f"{saddle_digitized_file_AD} not found")
                exit(1)
            if not os.path.exists(saddle_digitized_file_CT):
                print(f"{saddle_digitized_file_CT} not found")
                exit(1)
            if args.thread > 1:
                futures[pool.submit(process_chunks, cell_list_AD, bin_list_Mattew, find_bin_index_Mattew, cell_type, 'AD', saddle_digitized_file_AD, args.read_pair_folder)] = True
                futures[pool.submit(process_chunks, cell_list_CT, bin_list_Mattew, find_bin_index_Mattew, cell_type, 'CT', saddle_digitized_file_CT, args.read_pair_folder)] = True
            else:
                summary_table_AD, disease_group, cell_type = process_chunks(cell_list_AD, bin_list_Mattew, find_bin_index_Mattew, cell_type, 'AD', saddle_digitized_file_AD, args.read_pair_folder)
                print(f"Processing {disease_group} {cell_type} done")
                report(summary_table_AD, disease_group, cell_type, bin_list_Mattew, fout)
                summary_table_CT, disease_group, cell_type = process_chunks(cell_list_CT, bin_list_Mattew, find_bin_index_Mattew, cell_type, 'CT', saddle_digitized_file_CT, args.read_pair_folder)
                print(f"Processing {disease_group} {cell_type} done")
                report(summary_table_CT, disease_group, cell_type, bin_list_Mattew, fout)
        if args.thread > 1:
            # wait for all the futures to complete
            for future in as_completed(futures):
                summary_table, disease_group, cell_type = future.result()
                print(f"Processing {disease_group} {cell_type} done")
                report(summary_table, disease_group, cell_type, bin_list_Mattew, fout)
            # wait and close
            pool.shutdown(wait = True)
    #
    print(f"Done")
    

def report(summary_table, disease_group, cell_type, bin_list_Mattew, fout):
    for bin_idx_Mattew in sorted(summary_table.keys()):
        bin_start, bin_end = bin_list_Mattew[bin_idx_Mattew]
        for compartment_pair in sorted(summary_table[bin_idx_Mattew].keys()):
            count = summary_table[bin_idx_Mattew][compartment_pair]
            score_1, score_2 = compartment_pair
            print(f"{disease_group}\t{cell_type}\t{bin_idx_Mattew}\t{bin_start}\t{bin_end}\t{score_1}\t{score_2}\t{count}", file = fout)

def process_chunks(cell_list, bin_list_Mattew, find_bin_index_Mattew, cell_type, disease_group, saddle_digitized_file, read_pair_folder):
    """
    go through each read pairs
    stratify read pairs into bins based on the distance between two interacting loci
    and get the AB compartment percentile score for each loci
    """
    bin_size = 100000
    # load saddle digitized results
    bin2compartment = parse_saddle_digitized(saddle_digitized_file, bin_size)
    # main function to process the chunks
    result = {}
    #print(f"Processing {disease_group} {cell_type} with {len(cell_list)} cells")
    N = 0
    for cell in tqdm(cell_list):
        N += 1
        #if N % 3 == 0:
        #    break
        patient_id, cell_id = cell
        raw_pair_file_gz = os.path.join(read_pair_folder, str(patient_id), f"{cell_id}.pair.gz")
        raw_pair_file = os.path.join(read_pair_folder, str(patient_id), f"{cell_id}.pair")
        if not os.path.exists(raw_pair_file) and not os.path.exists(raw_pair_file_gz):
            print(f"{raw_pair_file} not found")
            continue
        # parse the raw pair file to get the distance between two interacting loci
        fin = gzip.open(raw_pair_file_gz, 'rt') if os.path.exists(raw_pair_file_gz) else open(raw_pair_file, 'r')
        for line in fin:
            row = line.strip().split('\t')
            chrom_1 = row[0]
            pos_1 = int(row[1])
            chrom_2 = row[2]
            pos_2 = int(row[3])
            # skip chrM
            if chrom_1 in ['chrM'] or chrom_2 in ['chrM']:
                continue
            # skip if the two loci are not on the same chromosome
            if chrom_1 != chrom_2:
                continue
            # get the distance between two interacting loci
            distance = abs(pos_2 - pos_1) 
            if distance > 1000:
                bin_idx_Mattew = find_bin_index_Mattew(bin_list_Mattew, distance)
            else:
                continue
            # get the AB compartment percentile score for each loci
            bin_1 = pos_1 // bin_size
            bin_2 = pos_2 // bin_size
            compartment_1 = bin2compartment[chrom_1][bin_1]
            compartment_2 = bin2compartment[chrom_2][bin_2]
            if compartment_1 == -1 or compartment_2 == -1:
                continue
            # update the summary table
            if result.get(bin_idx_Mattew, None) is None:
                result[bin_idx_Mattew] = {}
            compartment_pair = (compartment_1, compartment_2)
            if result[bin_idx_Mattew].get(compartment_pair, None) is None:
                result[bin_idx_Mattew][compartment_pair] = 0
            result[bin_idx_Mattew][compartment_pair] += 1
            compartment_pair_reverse = (compartment_2, compartment_1)
            if result[bin_idx_Mattew].get(compartment_pair_reverse, None) is None:
                result[bin_idx_Mattew][compartment_pair_reverse] = 0
            result[bin_idx_Mattew][compartment_pair_reverse] += 1
        fin.close()
    return result, disease_group, cell_type


def parse_saddle_digitized(saddle_digitized_file, bin_size):
    """
    parse the saddle digitized results
    """
    bin2compartment = dict()
    with open(saddle_digitized_file, 'r') as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            row = line.strip().split('\t')
            if row[0] == "chrom":
                continue
            chrom = row[0]
            start = int(row[1])
            start_bin = start // bin_size
            compartment_score = int(row[3])
            if bin2compartment.get(chrom, None) is None:
                bin2compartment[chrom] = dict()
            bin2compartment[chrom][start_bin] = compartment_score
    return bin2compartment


if __name__=="__main__":
    main()

