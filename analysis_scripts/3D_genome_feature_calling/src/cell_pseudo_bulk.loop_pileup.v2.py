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
    p.add_argument('--loop_folder',type=str,dest="loop_folder",help="loop folder")
    p.add_argument('--loop_processed_folder',type=str,dest="loop_processed_folder",help="folder to save the processed loop files")
    p.add_argument('--output_folder',type=str,dest="output_folder",help="output folder")
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
    #
    # get the loop files
    loop_table = {}
    for filename in os.listdir(args.loop_folder):
        if filename.endswith(".bedpe") and 'CTCF' in filename:
            loop_file = os.path.join(args.loop_folder, filename)
        else:
            continue
        # process the loop file
        # load the loop file
        data_loop = pd.read_csv(loop_file, sep = "\t", header = None)
        data_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'score', 'support', 'loop_type']
        data_loop['loop_size'] = abs((data_loop.start2 + data_loop.end2)/2 - (data_loop.start1 + data_loop.end1)/2)
        # split the loop by size into short-, mid-, and long-range loops and whether the loop_type is CTCF
        short_range_ctcf = data_loop[(data_loop.loop_size < 128000) & (data_loop.loop_type != 'NA-NA')]
        short_range_other = data_loop[(data_loop.loop_size < 128000) & (data_loop.loop_type == 'NA-NA')] 
        mid_range_ctcf = data_loop[(data_loop.loop_size >= 128000) & (data_loop.loop_size < 4096000) & (data_loop.loop_type != 'NA-NA')]
        mid_range_other = data_loop[(data_loop.loop_size >= 128000) & (data_loop.loop_size < 4096000) & (data_loop.loop_type == 'NA-NA')]
        long_range_ctcf = data_loop[(data_loop.loop_size >= 4096000) & (data_loop.loop_type != 'NA-NA')]
        long_range_other = data_loop[(data_loop.loop_size >= 4096000) & (data_loop.loop_type == 'NA-NA')]
        whole_ctcf = data_loop[data_loop.loop_type != 'NA-NA']
        whole_other = data_loop[data_loop.loop_type == 'NA-NA']
        whole = data_loop
        # save the processed loop files
        output_prefix = os.path.join(args.loop_processed_folder, filename.replace('.bedpe', ''))
        loop_label = filename.replace('.bedpe', '').replace('loop_summit.', '')
        if not os.path.exists(args.loop_processed_folder):
            os.makedirs(args.loop_processed_folder)
        if not os.path.exists(output_prefix + '.whole_ctcf.bedpe'):
            whole_ctcf.to_csv(output_prefix + '.whole_ctcf.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.whole_other.bedpe'):
            whole_other.to_csv(output_prefix + '.whole_other.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.whole.bedpe'):
            whole.to_csv(output_prefix + '.whole.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.short_range_ctcf.bedpe'):
            short_range_ctcf.to_csv(output_prefix + '.short_range_ctcf.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.short_range_other.bedpe'):
            short_range_other.to_csv(output_prefix + '.short_range_other.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.mid_range_ctcf.bedpe'):
            mid_range_ctcf.to_csv(output_prefix + '.mid_range_ctcf.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.mid_range_other.bedpe'):
            mid_range_other.to_csv(output_prefix + '.mid_range_other.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.long_range_ctcf.bedpe'):
            long_range_ctcf.to_csv(output_prefix + '.long_range_ctcf.bedpe', sep = "\t", header = True, index = False)
        if not os.path.exists(output_prefix + '.long_range_other.bedpe'):
            long_range_other.to_csv(output_prefix + '.long_range_other.bedpe', sep = "\t", header = True, index = False)
        # 
        loop_table[loop_label] = {'whole': output_prefix + '.whole.bedpe',
                                  'whole_ctcf': output_prefix + '.whole_ctcf.bedpe',
                                  'whole_other': output_prefix + '.whole_other.bedpe',
                                  'short_range_ctcf': output_prefix + '.short_range_ctcf.bedpe',
                                  'short_range_other': output_prefix + '.short_range_other.bedpe',
                                  'mid_range_ctcf': output_prefix + '.mid_range_ctcf.bedpe',
                                  'mid_range_other': output_prefix + '.mid_range_other.bedpe',
                                  'long_range_ctcf': output_prefix + '.long_range_ctcf.bedpe',
                                  'long_range_other': output_prefix + '.long_range_other.bedpe'}
        ## report the number of loops in each category
        print("Loop file: %s" % loop_file)
        print("Short-range CTCF loops: %d" % short_range_ctcf.shape[0])
        print("Short-range other loops: %d" % short_range_other.shape[0])
        print("Mid-range CTCF loops: %d" % mid_range_ctcf.shape[0])
        print("Mid-range other loops: %d" % mid_range_other.shape[0])
        print("Long-range CTCF loops: %d" % long_range_ctcf.shape[0])
        print("Long-range other loops: %d" % long_range_other.shape[0])
    # list the folder in args.folder
    pool = ProcessPoolExecutor(max_workers = args.thread)
    p_list = []
    target_folder_list = [item for item in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, item)) and item in ['AD', 'CT']]
    if args.cell_label == 'major':
        cell2cluster = {'Astro': ['Astro', 'all'], 'Endo': ['Endo', 'all'], 'Exc': ['Exc', 'all'], 'Inh': ['Inh', 'all'], 'Micro': ['Micro', 'all'], 'Oligo': ['Oligo', 'all'], 'OPC': ['OPC', 'all'], 'VLMC': ['VLMC', 'all']}
    else:
        cell2cluster = {'Astro': ['Astro', 'all'], 'Endo': ['Endo', 'all'], 'Exc_L2_3_IT': ['Exc_L2_3_IT', 'Exc', 'all'], 'Exc_L4_IT': ['Exc_L4_IT', 'Ex', 'all'], 'Exc_L5_6_NP': ['Exc_L5_6_NP', 'Exc', 'all'], 'Exc_L5_ET': ['Exc', 'all'], 'Exc_L5_IT': ['Exc_L5_IT', 'Exc', 'all'], 'Exc_L6b': ['Exc_L6b', 'Exc', 'all'], 'Exc_L6_CT': ['Exc_L6_CT', 'Exc', 'all'], 'Exc_L6_IT_Car3': ['Exc_L6_IT_Car3', 'Exc', 'all'], 'Exc_L6_IT': ['Exc_L6_IT', 'Exc','all'], 'Inh_Chandelier': ['Inh_Chandelier', 'Inh', 'all'], 'Inh_Lamp5': ['Inh_Lamp5', 'Inh', 'all'], 'Inh_PAX6': ['Inh_PAX6', 'Inh', 'all'], 'Inh_Pvalb': ['Inh_Pvalb', 'Inh', 'all'], 'Inh_Sncg': ['Inh_Sncg', 'Inh', 'all'], 'Inh_Sst': ['Inh_Sst', 'Inh', 'all'], 'Inh_Vip': ['Inh_Vip', 'Inh', 'all'], 'Micro': ['Micro', 'all'], 'Oligo': ['Oligo', 'all'], 'OPC': ['OPC', 'all'], 'VLMC': ['VLMC', 'all']}
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
                    c = cooler.Cooler(cool_file+'::resolutions/10000')
                except:
                    print("Error: cool file not valid: %s" % cool_file)
                    exit(1)
            # add cool file to the job list
            if args.cell_label == 'major':
                res = '10000'
            elif args.cell_label == 'sub':
                res = '50000'
            out_target_folder = os.path.join(args.output_folder, target_folder)
            if not os.path.exists(out_target_folder):
                os.makedirs(out_target_folder)
            for cluster_label in loop_table:
                # skip if the cluster_label is not in the cell2cluster for this cell type
                if not cluster_label.split('.')[0] in cell2cluster.get(cell_type, []):
                    continue
                #
                for loop_type in loop_table[loop_label]:
                    loop_file = loop_table[cluster_label][loop_type]
                    p = pool.submit(process_single_cool_file, cool_file, res, args.view_file, loop_file, cluster_label, loop_type, out_target_folder)
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
    

def process_single_cool_file(cool_file, res, view_file, loop_file, cluster_label, loop_type_label, output_folder):
    expected_file = cool_file.replace(".mcool", ".expected_cis.%dkb.tsv" % (int(res)/1000))
    output_file = os.path.join(output_folder, "%s.%s.%s.npz" % (os.path.basename(cool_file).replace(".mcool", ""), cluster_label, loop_type_label))
    if loop_type_label == 'whole':
        cmd = "cooltools pileup --features-format BEDPE --nproc 18 -o %s --store-snips --view %s --flank 200000 --expected %s %s::resolutions/%s %s" % (output_file, view_file, expected_file, cool_file, res, loop_file)
    else:
        cmd = "cooltools pileup --features-format BEDPE --nproc 18 -o %s --view %s --flank 200000 --expected %s %s::resolutions/%s %s" % (output_file, view_file, expected_file, cool_file, res, loop_file)
    return cmd

if __name__=="__main__":
    main()

