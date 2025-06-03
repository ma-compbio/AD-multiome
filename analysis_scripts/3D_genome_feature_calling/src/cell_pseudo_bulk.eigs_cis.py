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
            p = pool.submit(process_single_cool_file, cool_file, res, args.view_file, args.phasing_track)
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
    

def process_single_cool_file(cool_file, res, view_file, phasing_track):
    output_file = cool_file.replace(".mcool", ".%dkb" % (int(res)/1000))
    cmd = "cooltools eigs-cis -o %s --view %s --phasing-track %s --n-eigs 1 %s::resolutions/100000" % (output_file, view_file, phasing_track, cool_file)
    return cmd

if __name__=="__main__":
    main()
