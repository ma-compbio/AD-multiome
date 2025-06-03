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
    p.add_argument('--blacklist',type=str,dest="blacklist",help="blacklist file")
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
            cool_file = os.path.join(args.folder, target_folder, cell_type + '.cool')
            if not os.path.exists(cool_file):
                print("Error: cool file not found: %s" % cool_file)
                print("skip")
                continue
            else:
                print("Verify cool file: %s" % cool_file)
                try:
                    c = cooler.Cooler(cool_file)
                except:
                    print("Error: cool file not valid: %s" % cool_file)
                    exit(1)
            # add cool file to the job list
            if args.cell_label == 'major':
                res_list = "10000,20000,50000,100000,250000,500000"
            elif args.cell_label == 'sub':
                res_list = "10000,50000,100000,250000,500000"
            p = pool.submit(process_single_cool_file, cool_file, args.blacklist, res_list)
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
    

def process_single_cool_file(cool_file, blacklist_file, res_list):
    mcool_file = cool_file.replace(".cool", ".mcool")
    cmd = "cooler zoomify -p 24 --resolutions %s --balance --balance-args '--blacklist %s --ignore-diags 0 --nproc 32' -o %s %s" % (res_list, blacklist_file, mcool_file, cool_file)
    return cmd

if __name__=="__main__":
    main()

