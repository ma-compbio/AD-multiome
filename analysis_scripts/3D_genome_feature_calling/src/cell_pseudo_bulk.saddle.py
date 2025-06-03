#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
import cooler
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",help="number of threads to use",default=8)
    p.add_argument('--folder',type=str,dest="folder",help="pseudo bulk foldr")
    p.add_argument('--cell_anno',type=str,dest="cell_anno",help="cell annotation file")
    p.add_argument('--cell_label',type=str,dest="cell_label",choices=['major', 'sub'],help="cell label, either major or sub")
    p.add_argument('--view_file',type=str,dest="view_file",help="genome view file")
    p.add_argument('--phasing_track',type=str,dest="phasing_track",help="phasing track file in 100kb resolution")
    p.add_argument('--use_trans',dest="use_trans",action="store_true",help="use trans mode to run saddle")
    p.add_argument('--only_AD_and_CT',dest="only_AD_and_CT",action="store_true",help="only use AD and CT to run saddle")
    p.add_argument('--script_out',type=str,dest="script_out",help="output script file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def main():
    global args
    args = parse_arg()
    # load cell annotation file
    cell_anno = pd.read_csv(args.cell_anno, sep="\t", header=0) 
    if 'major' not in cell_anno.columns:
        # split the sub column and take the first part as major cell type
        cell_anno['major'] = cell_anno['sub'].apply(lambda x: x.split(" ")[0])
    # fix the space in the sub cell type
    cell_anno['sub'] = cell_anno['sub'].apply(lambda x: x.replace(" ", "_").replace("/", "_"))
    # 
    cell_type_list = cell_anno[args.cell_label].unique()
    print("Cell type list: %s" % ','.join(cell_type_list))
    # list the folder in args.folder
    pool = ProcessPoolExecutor(max_workers = args.thread)
    p_list = []
    target_folder_list = [item for item in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, item))]
    for target_folder in target_folder_list:
        if args.only_AD_and_CT and target_folder not in ['AD', 'CT']:
            continue
        for cell_type in cell_type_list:
            cool_file = os.path.join(args.folder, target_folder, cell_type + '.mcool')
            if not os.path.exists(cool_file):
                print("Error: cool file not found: %s" % cool_file)
                print("skip")
                continue
            else:
                print("Verify cool file: %s" % cool_file)
                try:
                    c = cooler.Cooler(cool_file+'::resolutions/100000')
                except:
                    print("Error: cool file at 100kb resolution not valid: %s" % cool_file)
                    exit(1)
            # add cool file to the job list
            res = '100000'
            p = pool.submit(process_single_cool_file, cool_file, res, args.view_file)
            p_list.append(p)
    # do the job
    cmd_list = []
    for p in tqdm(as_completed(p_list), total = len(p_list)):
        cmd = p.result()
        cmd_list.append(cmd)
    # wait and close
    pool.shutdown(wait = True)
    # write the script
    with open(args.script_out, 'w') as fout:
        for cmd in sorted(cmd_list):
            print(cmd, file = fout)
    # report to output 
    print("All finished")
    

def process_single_cool_file(cool_file, res, view_file):
    eigs_file = cool_file.replace(".mcool", ".%dkb.cis.vecs.tsv" % (int(res)/1000))
    if args.use_trans:
        expected_file = cool_file.replace(".mcool", ".expected_trans.%dkb.tsv" % (int(res)/1000))
        if not os.path.exists(expected_file):
            cmd = "cooltools expected-trans %s::resolutions/%s --nproc 18 -o %s --view %s" % (cool_file, res, expected_file, view_file)
            print(cmd)
        output_file = cool_file.replace(".mcool", ".saddle_trans.%dkb" % (int(res)/1000))
        cmd = "cooltools saddle --contact-type trans --qrange 0.02 0.98 --fig png -o %s --view %s %s::resolutions/%s %s %s" % (output_file, view_file, cool_file, res, eigs_file, expected_file)
    else:
        expected_file = cool_file.replace(".mcool", ".expected_cis.%dkb.tsv" % (int(res)/1000))
        if not os.path.exists(expected_file):
            cmd = "cooltools expected-cis %s::resolutions/%s --nproc 18 -o %s --view %s" % (cool_file, res, expected_file, view_file)
            print(cmd)
        output_file = cool_file.replace(".mcool", ".saddle.%dkb" % (int(res)/1000))
        cmd = "cooltools saddle --qrange 0.02 0.98 --fig png -o %s --view %s %s::resolutions/100000 %s %s" % (output_file, view_file, cool_file, eigs_file, expected_file) 
    return cmd

if __name__=="__main__":
    main()
