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
    p.add_argument('--use_trans',dest="use_trans",action="store_true",help="analyze saddle result from trans mode")
    p.add_argument('--output',type=str,dest="output",help="output script file")
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
    table = {}
    # get the AD saddle file
    # 100kb
    res = 100000
    with open(args.output, "w") as fout:
        print("disease_group\tcell_type\trow_idx\tcol_idx\tsaddle_value", file = fout)
        for cell_type in tqdm(cell_type_list):
            cool_file = os.path.join(args.folder, "AD", "%s.mcool" % cell_type)
            if args.use_trans:
                saddledump_file = cool_file.replace(".mcool", ".saddle_trans.%dkb.saddledump.npz" % (int(res)/1000))
            else:
                saddledump_file = cool_file.replace(".mcool", ".saddle.%dkb.saddledump.npz" % (int(res)/1000))
            if not os.path.exists(saddledump_file):
                print("Missing output file: %s" % saddledump_file)
            else:
                print(f"saddle file: {saddledump_file}")
            # parse the saddle file
            saddle = np.load(saddledump_file, allow_pickle = True)
            # get the saddle summary
            for row_idx in range(saddle['saddledata'].shape[0]):
                for col_idx in range(saddle['saddledata'].shape[1]):
                    value = saddle['saddledata'][row_idx, col_idx]
                    print("AD\t%s\t%d\t%d\t%.6f" % (cell_type, row_idx, col_idx, value), file=fout)
        for cell_type in tqdm(cell_type_list):
            cool_file = os.path.join(args.folder, "CT", "%s.mcool" % cell_type)
            if args.use_trans:
                saddledump_file = cool_file.replace(".mcool", ".saddle_trans.%dkb.saddledump.npz" % (int(res)/1000))
            else:
                saddledump_file = cool_file.replace(".mcool", ".saddle.%dkb.saddledump.npz" % (int(res)/1000))
            if not os.path.exists(saddledump_file):
                print("Missing output file: %s" % saddledump_file)
            else:
                print(f"saddle file: {saddledump_file}")
            # parse the saddle file
            saddle = np.load(saddledump_file, allow_pickle = True)
            # get the saddle summary
            for row_idx in range(saddle['saddledata'].shape[0]):
                for col_idx in range(saddle['saddledata'].shape[1]):
                    value = saddle['saddledata'][row_idx, col_idx]
                    print("CT\t%s\t%d\t%d\t%.6f" % (cell_type, row_idx, col_idx, value), file=fout)
        

    # report to output 
    print("All finished")


if __name__=="__main__":
    main()
