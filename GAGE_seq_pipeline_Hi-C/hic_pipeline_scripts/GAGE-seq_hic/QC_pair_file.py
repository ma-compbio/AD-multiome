#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import gzip
from tqdm import tqdm

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--pair_file',type=str,required=True,help="The pair file")
    p.add_argument('--mode',type=str,dest='mode',choices=['pair_type','pair_range', 'dup_level'],required=True,nargs="+",help="The mode of the pair file")
    p.add_argument('--label',type=str,dest='label',help="The label of the pair file")
    p.add_argument('--output',type=str,dest="output",help="The output prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def parse_pair_for_qc(pair_file, distance_group, do_pair_type=True, do_pair_range=True):
    # prep summary table
    result_pair_type = {}
    result_pair_range = {'inter':0, 'intra':dict((group,0) for group in distance_group)}
    result_pair_dup_level = {'inter':{}, 'intra': dict((group, {}) for group in distance_group)}
    for group in distance_group:
        result_pair_dup_level['intra'][group] = {}
    # if gzip
    if pair_file.endswith('.gz'):
        fin = gzip.open(pair_file,'rt')
    else:
        fin = open(pair_file)
    # read the file
    N = 1
    for line in fin:
        # the 8th columns is the pair type
        if N % 1000000 == 0:
            print(f"Processing {N} pairs", flush = True)
        row = line.strip().split('\t') 
        N += 1
        pair_id = row[0]
        chrom1 = row[1]
        pos1 = row[2]
        chrom2 = row[3]
        pos2 = row[4]
        pair_type = row[7]
        barcode = row[8]
        dup_count = row[9]
        if do_pair_type:
            if result_pair_type.get(pair_type, None) is None:
                result_pair_type[pair_type] = 0
            result_pair_type[pair_type] += 1
        # if pair type in ['UU', 'UR', 'RU', 'RR']
        if do_pair_range:
            if chrom1 == '!' or chrom2 == '!':
                continue
            if chrom1 != chrom2:
                result_pair_range['inter'] += 1
                if result_pair_dup_level['inter'].get(dup_count, None) is None: 
                    result_pair_dup_level['inter'][dup_count] = 0
                result_pair_dup_level['inter'][dup_count] += 1
            else:
                pair_group = 'intra'
                # for pair range only include intra
                if pair_type in ['UU', 'UR', 'RU', '..']: 
                    distance = abs(int(pos2) - int(pos1))
                    if distance < 1000:
                        distance_group = '<1kb'
                    elif distance >= 1000 and distance < 5000:
                        distance_group = '1kb-5kb'
                    elif distance >= 5000 and distance < 10000:
                        distance_group = '5kb-10kb'
                    elif distance >= 10000 and distance < 25000:
                        distance_group = '10kb-25kb'
                    elif distance >= 25000 and distance < 50000:
                        distance_group = '25-50kb'
                    elif distance >= 50000 and distance < 100000:
                        distance_group = '50kb-100kb'
                    elif distance >= 100000 and distance < 250000:
                        distance_group = '100kb-250kb'
                    elif distance >= 250000 and distance < 500000:
                        distance_group = '250kb-500kb'
                    elif distance >= 500000 and distance < 1000000:
                        distance_group = '500kb-1mb'
                    elif distance >= 1000000 and distance < 5000000:
                        distance_group = '1mb-5mb'
                    elif distance >= 5000000 and distance < 10000000:
                        distance_group = '5mb-10mb'
                    elif distance >= 10000000 and distance < 50000000:
                        distance_group = '10mb-50mb'
                    else:
                        distance_group = '>50mb'
                    # update
                    result_pair_range[pair_group][distance_group] += 1
                    if result_pair_dup_level[pair_group][distance_group].get(dup_count, None) is None:
                        result_pair_dup_level[pair_group][distance_group][dup_count] = 0
                    result_pair_dup_level[pair_group][distance_group][dup_count] += 1
    fin.close()
    return result_pair_type, result_pair_range, result_pair_dup_level


def main():
    global args
    args = parse_arg()
    if 'pair_type' in args.mode:
        do_pair_type = True
        print("Processing pair type")
    else:
        print("Skip pair type")
        do_pair_type = False
    if 'pair_range' in args.mode:
        do_pair_range = True
        print("Processing pair range")
    else:
        print("Skip pair range") 
        do_pair_range = False
    if 'dup_level' in args.mode:
        do_dup_level = True
        print("Processing dup level")
    else:
        print("Skip dup level")
        do_dup_level = False
    # process
    distance_group = ['<1kb', '1kb-5kb', '5kb-10kb', '10kb-25kb', '25-50kb', '50kb-100kb', '100kb-250kb', '250kb-500kb', '500kb-1mb', '1mb-5mb', '5mb-10mb', '10mb-50mb', '>50mb']
    result_pair_type, result_pair_range, result_pair_dup_level = parse_pair_for_qc(args.pair_file, distance_group, do_pair_type, do_pair_range)
    print("Finish processing")
    # report
    if do_pair_type:
        with open(args.output + '.qc.pair_type.tsv', 'w') as fout:
            print("label\tpair_type\tcount", file = fout)
            for pair_type, count in sorted(result_pair_type.items()):
                print(f"{args.label}\t{pair_type}\t{count}", file = fout)
        print("Finish pair type")
    if do_pair_range:
        with open(args.output + '.qc.pair_range.tsv', 'w') as fout:
            print("label\tpair_type\trange\tcount", file = fout)
            print(f"{args.label}\tinter\tinter\t{result_pair_range['inter']}", file = fout)
            for range_group in distance_group:
                print(f"{args.label}\tintra\t{range_group}\t{result_pair_range['intra'][range_group]}", file = fout)
        print("Finish pair range")
    if do_dup_level:
        with open(args.output + '.qc.pair_duplevel.tsv', 'w') as fout:
            print("label\tpair_type\trange\tdup_level\tcount", file = fout)
            for dup in sorted(result_pair_dup_level['inter'].keys()):
                print(f"{args.label}\tinter\tinter\t{dup}\t{result_pair_dup_level['inter'][dup]}", file = fout)
            for range_group in distance_group:
                for dup in sorted(result_pair_dup_level['intra'][range_group].keys()):
                    print(f"{args.label}\tintra\t{range_group}\t{dup}\t{result_pair_dup_level['intra'][range_group][dup]}", file = fout)
        print("Finish dup level")
    

if __name__=="__main__":
    main()

