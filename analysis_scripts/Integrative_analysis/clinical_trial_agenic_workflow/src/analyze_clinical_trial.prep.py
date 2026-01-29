#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import json
from tqdm import tqdm
import re


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_prefix',type=str,dest='input_prefix',help='Input file prefix',default=None)
    p.add_argument('--output_prefix',type=str,dest='output_prefix',help='Output file prefix',default=None)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def main():
    global args
    args = parse_arg()
    trial_json_file = args.input_prefix + ".json"
    if not os.path.exists(trial_json_file):
        print("Error: Input file not exists!")
        exit(1)
    # Load the JSON file
    trial_list = []
    with open(trial_json_file,'r') as f:
        trial_list = json.load(f)
    print(f"Total {len(trial_list)} trials loaded from {trial_json_file}")
    # prepare the output file
    output_file_summary = args.output_prefix + ".summary.tsv"
    output_file_drug = args.output_prefix + ".drug.tsv"
    # for each trial, report drug and target gene in a long format
    fout_drug = open(output_file_drug,'w')
    print(f"nct_id\ttarget_category\tdrug_agent_type\tdrug_target_type\tdrug_name\tdrugbank_name\tdrugbank_id\tdrug_kingdom\tpathway", file=fout_drug)
    with open(output_file_summary,'w') as fout:
        print(f"nct_id\tstatus\tlast_update\tphase\tstudy_type\tprimary_purpose\ttarget_category\tdrug_agent_type\tdrug_target_type\tdrug_name\ttarget_gene\tconfidence_level\tconfidence_level_clean", file=fout)
        for trial in tqdm(trial_list):
            nct_id = trial.get("nct_id","NA")
            status = trial.get("status","NA")
            last_update = trial.get("last_update","NA")
            phase = trial.get("phase","NA")
            phase = phase if phase != 'na' else 'NA'
            study_type = trial.get("study_type","NA")
            primary_purpose = trial.get("primary_purpose","NA")
            target_category = trial.get("target_category","NA")
            drug_agent_type = trial.get("agent_type","NA")
            drug_target_type = trial.get("target_group","NA")
            drug_name = trial.get("drugs","NA")
            drug_object = trial.get("drug_objects", [])
            if len(drug_object) > 0:
                drug_id_list = []
                drug_name_list = []
                drug_kingdom_list = []
                drug_pathway_list = []
                for drug in drug_object:
                    drug_kingdom_list.append(drug["drug_kingdom"])
                    drug_pathway_list.extend(drug["pathway"])
                    drug_id_list.append(drug["drug_id"])
                    drug_name_list.append(drug["drug_name"])
                for drug_idx in range(len(drug_id_list)):
                    if len(drug_pathway_list) == 0:
                        print(f"{nct_id}\t{target_category}\t{drug_agent_type}\t{drug_target_type}\t{drug_name}\t{drug_name_list[drug_idx]}\t{drug_id_list[drug_idx]}\t{drug_kingdom_list[drug_idx]}\tNA", file=fout_drug)
                    for pathway_idx in range(len(drug_pathway_list)):
                        print(f"{nct_id}\t{target_category}\t{drug_agent_type}\t{drug_target_type}\t{drug_name}\t{drug_name_list[drug_idx]}\t{drug_id_list[drug_idx]}\t{drug_kingdom_list[drug_idx]}\t{drug_pathway_list[pathway_idx]}", file=fout_drug)
                #
            target_list = trial.get("targets",[])
            confidence_level = trial.get("confidence_level","NA")
            confidence_level = confidence_level if confidence_level not in [None, 'null'] else 'NA'
            # clean confidence level
            # remove content in parentheses
            confidence_level_clean = re.sub(r"\(.*?\)","",confidence_level).strip()
            # remove letters
            confidence_level_clean = re.sub(r"[a-zA-Z]","",confidence_level_clean).strip()
            # remove space
            confidence_level_clean = re.sub(r"\s+"," ",confidence_level_clean).strip()
            # remove gene_name;
            confidence_level_clean = re.sub(r"[a-zA-Z0-9]*:","",confidence_level_clean).strip()
            # remove ;
            confidence_level_clean = re.sub(r"[;]","",confidence_level_clean).strip()
            # remove anything after '/'
            confidence_level_clean = re.sub(r"/.*$","",confidence_level_clean).strip()
            # keep the first number
            confidence_level_clean = confidence_level_clean.split()[0] if len(confidence_level_clean.split()) > 0 else "NA"
            if len(target_list) == 0:
                continue
            # for each target, report a line
            for target_gene in target_list:
                gene_symbol = target_gene.get("gene_symbol","NA")
                if gene_symbol not in ['null', 'NA', None]:
                    print(f"{nct_id}\t{status}\t{last_update}\t{phase}\t{study_type}\t{primary_purpose}\t{target_category}\t{drug_agent_type}\t{drug_target_type}\t{drug_name}\t{gene_symbol}\t{confidence_level}\t{confidence_level_clean}", file=fout)
    fout_drug.close()
    
    
if __name__=="__main__":
    main()

