#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import numpy as np
from tqdm import tqdm


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--folder',type=str,dest="folder",help="The folder containing the pileup files",required=True)
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def main():
    global args
    args = parse_arg()
    # load the filename in the folder
    table = {}
    with open(args.output, 'w') as fout:
        print("cell_type\tcluster_id\tctcf_type\tloop_type\trow_idx\tcol_idx\tvalue", file=fout)
        for filename in tqdm(os.listdir(args.folder)):
            if filename.endswith(".npz"):
                pileup_file = os.path.join(args.folder,filename)
                # parse the file name to get different information
                cell_type, cluster_id, ctcf_type, loop_type = filename.replace('.npz', '').split('.')
                # load the npz file
                pile = np.load(pileup_file)
                pileup = pile['pileup']
                # convert the 2D matrix into long format
                for row_idx in range(pileup.shape[0]):
                    for col_idx in range(pileup.shape[1]):
                        print(f"{cell_type}\t{cluster_id}\t{ctcf_type}\t{loop_type}\t{row_idx}\t{col_idx}\t{pileup[row_idx, col_idx]:.6f}", file=fout)

    
if __name__=="__main__":
    main()

