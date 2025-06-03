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
    # get the patient list from cell_anno, patient list is the disease_group plut projid column 
    agg_table = cell_anno[['disease_group', 'gender', 'projid']].drop_duplicates().reset_index(drop = True)
    # verify patient specific cool file exist
    for index, row in agg_table.iterrows():
        disease_group = row['disease_group']
        projid = row['projid']
        patient_folder = os.path.join(args.folder, "%s_patient_%s" % (disease_group, projid))
        for cell_type in cell_type_list:
            cool_file = os.path.join(patient_folder, cell_type + '.cool')
            if not os.path.exists(cool_file):
                print("Error: cool file not found: %s" % cool_file)
                print("skip")
            else:
                print("Verify cool file: %s" % cool_file)
                try:
                    c = cooler.Cooler(cool_file)
                except:
                    print("Error: cool file not valid: %s" % cool_file)
                    exit(1)
    # add cool file to the merge list
    merge_gender = {}
    merge = {}
    for index, row in agg_table.iterrows():
        disease_group = row['disease_group']
        gender = row['gender']
        projid = row['projid']
        patient_folder = os.path.join(args.folder, "%s_patient_%s" % (disease_group, projid))
        for cell_type in cell_type_list:
            # coolfile
            cool_file = os.path.join(patient_folder, cell_type + '.cool')
            # skip if the cool file not exist
            if not os.path.exists(cool_file):
                print(f"skip {cool_file}")
                continue
            key_gender = "%s_%s" % (disease_group, gender)
            key = "%s" % disease_group
            if merge_gender.get(key_gender, None) is None:
                merge_gender[key_gender] = {}
                if not os.path.exists(os.path.join(args.folder, key_gender)):
                    os.mkdir(os.path.join(args.folder, key_gender))
            if merge.get(key, None) is None:
                merge[key] = {}
                if not os.path.exists(os.path.join(args.folder, key)):
                    os.mkdir(os.path.join(args.folder, key))
            #
            if merge_gender[key_gender].get(cell_type, None) is None:
                merge_gender[key_gender][cell_type] = []
            if merge[key].get(cell_type, None) is None:
                merge[key][cell_type] = []
            # add cool file to the list
            merge_gender[key_gender][cell_type].append(cool_file)
            merge[key][cell_type].append(cool_file)
    # do the jobs
    pool = ProcessPoolExecutor(max_workers = args.thread)
    p_list = []
    for key in tqdm(sorted(merge.keys())):
        for cell_type in merge[key]:
            cool_file_list = merge[key][cell_type]
            cool_out_file = os.path.join(args.folder, key, cell_type + '.cool')
            p = pool.submit(merge_single_cell_type, cool_out_file, cool_file_list)
            p_list.append(p)
    for key in tqdm(sorted(merge_gender.keys())):
        for cell_type in merge_gender[key]:
            cool_file_list = merge_gender[key][cell_type]
            cool_out_file = os.path.join(args.folder, key, cell_type + '.cool')
            p = pool.submit(merge_single_cell_type, cool_out_file, cool_file_list)
            p_list.append(p)
    # do the job
    for p in tqdm(as_completed(p_list), total = len(p_list)):
        p.result()
    # wait and close
    pool.shutdown(wait = True)
    # report to output 
    print("All finished")


def merge_single_cell_type(cool_out_file, cool_file_list):
    # merge
    cooler_cmd = "cooler merge %s %s" % (cool_out_file, " ".join(cool_file_list))
    os.system(cooler_cmd)
    
if __name__=="__main__":
    main()

