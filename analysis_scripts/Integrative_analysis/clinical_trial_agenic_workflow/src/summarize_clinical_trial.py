#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
from tqdm import tqdm
import json


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_json_folder',type=str,dest="input_json_folder",help="input json folder",default=None,required=True)
    p.add_argument('--output_summary',type=str,dest="output_summary",help="output summary file",default=None,required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def main():
    global args
    args = parse_arg()
    # list the json files in the input folder
    json_list = []
    for folder in tqdm(os.listdir(args.input_json_folder)):
        # check if it is a folder
        folder_path = os.path.join(args.input_json_folder, folder)
        trial_json_list = []
        if os.path.isdir(folder_path):
            net_id = folder
            # list json file in the folder
            for filename in os.listdir(folder_path):
                if filename.endswith(".json"):
                    run_id = filename.split(".")[-2]
                    trial_json_list.append((run_id, os.path.join(folder_path, filename)))
        else:
            continue
        # select the json file with the largest run id
        trial_json_list = sorted(trial_json_list, key=lambda x: int(x[0]), reverse=True)
        if len(trial_json_list) > 1:
            print(f"Processing {net_id} with {trial_json_list[0][1]}")
        # add to the json list
        json_list.append(trial_json_list[0][1])
    print(f"Total {len(json_list)} json files to process")
    # load all json files
    json_list = [json.load(open(json_file, 'r')) for json_file in json_list]
    # sort the json list by nct_id 
    json_list = sorted(json_list, key=lambda x: x['nct_id'])
    # process each json file
    with open(args.output_summary, 'w') as fout:
        print(f"Writing summary to {args.output_summary}")
        print("nct_id\tstatus\tlast_update_time\tphase\tstudy_type\tprimary_purpose\ttarget_category\tagent_type\tdrug", file = fout)
        for trial in tqdm(json_list):
            nct_id = trial['nct_id']
            status = trial['status'].lower()
            phase = ';'.join(trial['phase']).lower() if isinstance(trial['phase'], list) else trial['phase'].lower()
            study_type = trial['study_type'].lower()
            primary_purpose = trial['primary_purpose'].lower() if 'primary_purpose' in trial else 'NA'
            target_category = trial['target_category'] if trial['target_category'] != 'N/A' else 'NA'
            agent_type = trial['agent_type'] if trial['agent_type'] != 'N/A' else 'NA'
            # how to hide comma in drug name
            drug = ';'.join(trial['drug']).lower().replace(';', ' ') if len(trial['drug']) > 0 else 'NA' 
            print(f"{nct_id}\t{status}\t{trial['last_update_time']}\t{phase}\t{study_type}\t{primary_purpose}\t{target_category}\t{agent_type}\t{drug}", file = fout)
        
    
if __name__=="__main__":
    main()

