#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import pysam

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bam_file',type=str,dest="bam_file",help="The bam file")
    p.add_argument('--label',type=str,dest="label",help="The label of the bam file")
    p.add_argument('--output',type=str,dest="output",help="output file prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def main():
    global args
    args = parse_arg()
    fin = pysam.AlignmentFile(args.bam_file,'rb')
    table = {'read_1': {'n_primary': 0, 'n_30', 'n_multi': 0, 'n_multi_30': 0, 'n_chimeric': 0, 'n_chimeric_30': 0, 'n_chimeric_all': 0, 'n_chimeric_30_all': 0, 'n_unique': 0, 'n_unique_30': 0}, 'read_2':{'n_primary': 0, 'n_30', 'n_multi': 0, 'n_multi_30': 0, 'n_chimeric': 0, 'n_chimeric_30': 0, 'n_chimeric_all': 0, 'n_chimeric_30_all': 0, 'n_unique': 0, 'n_unique_30': 0}}
    #
    for read in fin.fetch(until_eof=True):
        # primary or secondary
        is_primary = True
        if read.is_secondary:
            is_primary = False
        # read quality
        high_quality = True
        if read.mapping_quality < 30:
            high_quality = False
        # XA tag
        is_multi = False
        if read.has_tag('XA'):
            is_multi  = True
        # chimena
        is_chimeric = False
        if read.has_tag('SA'):
            is_chimeric = True
        # first or second in pair
        read_group= 'read_1'
        if read.is_read2:
            read_group = 'read_2'
        # assign to table 
        if is_primary:
            table[read_group]['n_primary'] += 1
            if high_quality:
                table[read_group]['n_30'] += 1
            if is_multi:
                table[read_group]['n_multi'] += 1
                if high_quality:
                    table[read_group]['n_multi_30'] += 1
            if is_chimeric:
                table[read_group]['n_chimeric'] += 1
                if high_quality:
                    table[read_group]['n_chimeric_30'] += 1
            if not is_chimeric:
                table[read_group]['n_unique'] += 1
                if high_quality:
                    table[read_group]['n_unique_30'] += 1
        # chimeric all
        if is_chimeric:
            table[read_group]['n_chimeric_all'] += 1
            if high_quality:
                table[read_group]['n_chimeric_30_all'] += 1
    fin.close()
    # report to output
    with open(args.output, 'w') as fout:
        print("label\tentry\tvalue")
        for read_group in ['read_1', 'read_2']:
            for label in ['n_primary', 'n_30', 'n_multi', 'n_multi_30', 'n_chimeric', 'n_chimeric_30', 'n_chimeric_all', 'n_chimeric_30_all', 'n_unique', 'n_unique_30']:
                print(f"{args.label}\t{label}\t{table[read_group][label]}", file = fout)
    
if __name__=="__main__":
    main()

