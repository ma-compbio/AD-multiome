#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 09 Oct 2025 03:50:38 PM

import os,sys,argparse
import json
from openai import OpenAI
from pydantic import BaseModel
from tqdm import tqdm
import random
import time
import pandas as pd
from thefuzz import fuzz
from thefuzz import process # for fuzzy string matching and comparison
import re


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_clinical_trial',type=str,dest="input_trial",help="clinical trial summary file")
    p.add_argument('--input_drugbank',type=str,dest="drugbank",help="drugbank parsed file from xml")
    p.add_argument('--input_trial2drug',type=str,dest="trial2durg",help="trial to drug relationship parsed from DrugBank")
    p.add_argument('--local_database_prefix',type=str,dest="local_db_prefix",help="local db prefix used to save information for drugs")
    p.add_argument('--input_openai_folder',type=str,dest="input_openai_folder",help="input folder to load trial's JSON results")
    p.add_argument('--output_openai_folder',type=str,dest="output_openai_folder",help="output folder to save OpenAI results")
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


class Drug(object):
    def __init__(self):
        self.drug_id = None # DrugBank ID
        self.drug_name = "NA" 
        self.drug_synonyms = [] # list of synonym name
        self.drug_kingdom = None # Organic Compounds, etc
        self.drug_superclass = None # Organic Acids, etc
        self.pathway = [] # pathway information
        self.targets = [] # list of targets 
    def add_pathway(self, pathways):
        """
        Add pathway information to the drug
        """
        # check if the pathway is already in the list
        for p in self.pathway:
            if p in pathways:
                return
        # if not, add it
        self.pathway.extend(pathways)
    def add_target(self, target):
        """
        Add a target object to the drug
        """
        # check if the target is already in the list
        for t in self.targets:
            if t.target_name == target.target_name:
                return
        # if not, add it
        self.targets.append(target)
    # enable JSON serialization
    def to_dict(self):
        return {
            'drug_id': self.drug_id,
            'drug_name': self.drug_name,
            'drug_synonyms': self.drug_synonyms,
            'drug_kingdom': self.drug_kingdom,
            'drug_superclass': self.drug_superclass,
            'pathway': self.pathway,
            'targets': [t.to_dict() for t in self.targets]
        }


class Target(object):
    def __init__(self):
        self.target_name = None
        self.target_type = None # protein,
        self.uniprot_id = None
        self.gene_symbol = None # gene symbol
        self.gene_name = None # gene name description
    # enable JSON serialization
    def to_dict(self):
        return {
            'target_name': self.target_name,
            'target_type': self.target_type,
            'uniprot_id': self.uniprot_id,
            'gene_symbol': self.gene_symbol,
            'gene_name': self.gene_name
        }
 

class Trial(object):
    def __init__(self, feature):
        self.nct_id = feature['nct_id']
        self.status = feature['status']
        self.last_update_time = feature['last_update_time']
        self.phase = feature['phase']
        self.study_type = feature['study_type']
        self.primary_purpose = feature['primary_purpose']
        self.target_category = feature['target_category'].capitalize() if feature['target_category'].strip() != 'NA' else "NA"
        self.agent_type = feature['agent_type']
        self.drugs = [drug for drug in feature['drug'].split(';') if drug.strip() != 'NA'] # list of drug names
        # associated Drug objects
        self.drug_objects = [] # list of Drug objects
        # 
        self.target_group = None # choose from DrugBank_direct, DrugBank_match, OpenAI, None
        self.targets = [] # list of Target objects associated with the trial
        self.confidence_level = None # confidence level from OpenAI if applicable
    def to_dict(self):
        return {
            'nct_id': self.nct_id,
            'status': self.status,
            'last_update_time': self.last_update_time,
            'phase': self.phase,
            'study_type': self.study_type,
            'primary_purpose': self.primary_purpose,
            'target_category': self.target_category,
            'agent_type': self.agent_type,
            'drugs': self.drugs,
            'drug_objects': [d.to_dict() for d in self.drug_objects],
            'target_group': self.target_group,
            'targets': [t.to_dict() for t in self.targets],
            'confidence_level': self.confidence_level
        }
    def build_target_from_drug(self):
        """
        When the drug information is available, build target information from the drug objects
        """
        for drug in self.drug_objects:
            for target in drug.targets:
                # check if the target is already in the list
                is_exist = False
                for t in self.targets:
                    if t.target_name == target.target_name:
                        is_exist = True
                        break
                if not is_exist:
                    self.targets.append(target)


class TargetGeneExtraction(BaseModel):
    gene_symbol: list[str]
    confidence_level: str 
    explanation: list[str]


def load_clinical_trial_summary(input_file):
    """
    load clinical trial summary file
    nct_id,status,last_update_time,phase,bstudy_type,primary_purpose,target_category,agent_type,drug
    NCT00000171,completed,2005-06-23,phase3,interventional,treatment,neuropsychiatric symptom improvement,O) Circadian Rhythm,melatonin (n-acetyl-5-methoxytryptamine)
    NCT00000172,completed,2005-06-23,phase3,interventional,treatment,cognitive enhancer,D) Neurotransmitter Receptors,galantamine
    """
    header = {}
    trial_list = []
    with open(input_file, 'r') as fin:
        for line in fin:
            row = line.strip().split('\t')
            if len(header) < 1:
                for nn, col in enumerate(row):
                    header[col] = nn
                continue
            # process each line
            trial = {}
            for col in header:
                trial[col] = row[header[col]]
            trial_list.append(Trial(trial))
    print(f"Loaded {len(trial_list)} clinical trials from {input_file}")
    return trial_list


def build_drug_local_db(drugbank_df, local_db_prefix):
    """
    Create a local database for drugs in drugbank 
    """
    drug_table = {}
    for index, row in tqdm(drugbank_df.iterrows()):
        # dg_id,dg_name,dg_synm,dg_kingdom,dg_superclass,dg_pathways,target_name,target_uniprot,target_gene_name,action,cell_loc
        #DB00001,Lepirudin,"[Leu1, Thr2]-63-desulfohirudin;Desulfatohirudin;Hirudin variant-1;Lepirudin;Lepirudin recombinant;R-hirudin",Organic Compounds,Organic Acids,Lepirudin Action Pathway,Prothrombin,P00734,F2,inhibitor,"Secreted, extracellular space"
        drug_id = row['dg_id']
        drug_name = row['dg_name']
        drug_synonyms = row['dg_synm'].split(';') if not pd.isna(row['dg_synm'])  else []
        drug_kingdom = row['dg_kingdom']
        drug_superclass = row['dg_superclass']
        drug_pathways = row['dg_pathways'].split(';') if not pd.isna(row['dg_pathways']) else []
        target_name = row['target_name']
        target_uniprot = row['target_uniprot']
        target_gene_name = row['target_gene_name']
        # create Drug object
        if drug_table.get(drug_id, None) is None:
            drug_obj = Drug()
            drug_obj.drug_id = drug_id
            drug_obj.drug_name = drug_name
            drug_obj.drug_kingdom = drug_kingdom
            drug_obj.drug_superclass = drug_superclass
            drug_obj.drug_synonyms = drug_synonyms
            drug_table[drug_id] = drug_obj
        drug_obj = drug_table[drug_id]
        # update pathyway & target information
        if len(drug_pathways) > 0:
            drug_obj.add_pathway(drug_pathways)
        if not pd.isna(target_name):
            target_obj = Target()
            target_obj.target_name = target_name
            target_obj.uniprot_id = target_uniprot if not pd.isna(target_uniprot) else None
            target_obj.gene_symbol = target_gene_name if not pd.isna(target_gene_name) else None
            drug_obj.add_target(target_obj)
    # save to local db file as json
    drug_table_json = {}
    for drug_id, drug in drug_table.items():
        drug_table_json[drug_id] = drug.to_dict()
    # dump to json file
    json.dump(drug_table_json, open(local_db_prefix + '.json', 'w'), indent=4)
    print(f"Built local drug database with {len(drug_table)} drugs and saved to {local_db_prefix}.json")
    return drug_table
        

def create_drug_name_fuzzy_search_list(drug_database):
    """
    Create a fuzzy search list for drug names and synonyms
    """
    drug_name_search_list = []
    name2drug = {}
    for drug_id, drug in drug_database.items():
        # add drug name
        drug_name_search_list.append(drug.drug_name.lower())
        name2drug[drug.drug_name.lower()] = drug.drug_id
        # add synonyms
        for syn in drug.drug_synonyms:
            if len(syn.strip()) <= 3:
                continue
            if syn.lower() in name2drug:
                continue
            drug_name_search_list.append(syn.lower())
            name2drug[syn.lower()] = drug.drug_id
    print(f"Created fuzzy search list with {len(drug_name_search_list)} drug names and synonyms")
    return drug_name_search_list, name2drug


def main():
    global args
    args = parse_arg()
    # load drugbank file
    drugbank_df = pd.read_csv(args.drugbank)
    print(f"Loaded {len(drugbank_df)} drugs from DrugBank: {args.drugbank}")
    drug_local_db = build_drug_local_db(drugbank_df, args.local_db_prefix)
    # create drug name fuzzy matching list 
    drug_name_search_list, drug_name2id = create_drug_name_fuzzy_search_list(drug_local_db)
    print(f"Drug name fuzzy search list created with {len(drug_name_search_list)} entries")
    # load input clinical trial summary file
    trial_list = load_clinical_trial_summary(args.input_trial)
    # load trial to drug relationship file
    trial2drug_df = json.load(open(args.trial2durg, 'r'))
    print(f"Loaded {len(trial2drug_df)} trial to drug relationships from {args.trial2durg}")
    #
    client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    if not os.path.exists(args.output_openai_folder):
        os.makedirs(args.output_openai_folder)
    # process each trial
    count_table = {}
    for trial in tqdm(trial_list):
        if trial.target_category.upper() == 'NA':
            count_table['No_target_category'] = count_table.get('No_target_category', 0) + 1
            continue
        # Senerio 1: DrugBank_direct
        if trial.nct_id in trial2drug_df:
            for drug_name in trial2drug_df[trial.nct_id]['Drug(s)']:
                drug_id = trial2drug_df[trial.nct_id]['Drug(s)'][drug_name]
                if drug_id in drug_local_db:
                    trial.drug_objects.append(drug_local_db[drug_id])
            if len(trial.drug_objects) > 0:
                trial.target_group = 'DrugBank_direct'
                count_table['DrugBank_direct'] = count_table.get('DrugBank_direct', 0) + 1
                trial.build_target_from_drug()
                continue
        # Senerio 2: DrugBank_match
        if len(trial.drugs) > 0:
            for drug_name in trial.drugs:
                # remove content in parentheses
                drug_name = re.sub(r'\(.*?\)', '', drug_name).strip()
                # fuzzy search
                best_match_list = process.extract(drug_name.lower(), drug_name_search_list, limit = 2, scorer=fuzz.token_sort_ratio)
                print(f"Trial {trial.nct_id} drug {drug_name} best matches: {best_match_list}")
                # count how many numeric letters in the drug_name
                num_numeric = sum(c.isdigit() for c in drug_name)
                fraction_numeric = num_numeric / len(drug_name) if len(drug_name) > 0 else 0
                if fraction_numeric > 0.5:
                    cutoff = 95
                else:
                    cutoff = 90
                for match_name, score in best_match_list:
                    if score >= cutoff:
                        drug_id = drug_name2id[match_name]
                        if drug_id in drug_local_db:
                            # check if the drug is already in the trial.drug_objects
                            is_exist = False
                            for d in trial.drug_objects:
                                if d.drug_id == drug_id:
                                    is_exist = True
                                    break
                            if not is_exist:
                                trial.drug_objects.append(drug_local_db[drug_id])
            if len(trial.drug_objects) > 0:
                trial.target_group = 'DrugBank_match'
                count_table['DrugBank_match'] = count_table.get('DrugBank_match', 0) + 1
                trial.build_target_from_drug()
            else:
                # Senerio 3: OpenAI
                trial.target_group = 'OpenAI'
                count_table['OpenAI'] = count_table.get('OpenAI', 0) + 1
                # load existing OpenAI result if exists
                trial_json_file = os.path.join(args.input_openai_folder, trial.nct_id, trial.nct_id + '.1.json')
                if not os.path.exists(trial_json_file):
                    print(f"OpenAI result file not found for trial {trial.nct_id}: {trial_json_file}")
                    continue
                trial_json = json.load(open(trial_json_file, 'r'))
                # get brief summary and agent type
                trial_title = trial_json.get('title', '')
                trial_drug = trial_json.get('drug', '')
                # prepare output json file
                trial_target_json_file = os.path.join(args.output_openai_folder, trial.nct_id + '.json')
                # if the output file already exists, load it
                if os.path.exists(trial_target_json_file):
                    trial_target_json = json.load(open(trial_target_json_file, 'r'))
                else:
                    print(f"Processing trial {trial.nct_id} with OpenAI, title: {trial_title}, drug: {trial_drug}")
                    try:
                        result = openai_get_target_gene_from_summary_react(client, {'nct_id': trial.nct_id, 'drug': trial_drug, 'title': trial_title})
                    except Exception as e:
                        print(f"Error processing trial {trial.nct_id} with OpenAI: {e}")
                        continue
                    # add result to trial_target_json
                    trial_json['gene_symbol'] = result['gene_symbol']
                    trial_json['confidence_level'] = result['confidence_level']
                    trial_json['explanation_target_gene'] = result['explanation_target_gene']
                    # save to output file
                    json.dump(trial_json, open(trial_target_json_file, 'w'), indent=4)
                    trial_target_json = json.load(open(trial_target_json_file, 'r'))
                    # sleep for a while to avoid rate limit
                    time.sleep(random.uniform(1, 10))
                # add target information to trial object 
                if 'gene_symbol' in trial_target_json and len(trial_target_json['gene_symbol']) > 0 and trial_target_json['gene_symbol'][0] != 'N/A':
                    for gene in trial_target_json['gene_symbol']:
                        target_obj = Target()
                        target_obj.gene_symbol = gene
                        trial.targets.append(target_obj)
                # add confidence level
                if 'confidence_level' in trial_target_json:
                    trial.confidence_level = trial_target_json['confidence_level']
        else:
            count_table['No_drug_info'] = count_table.get('No_drug_info', 0) + 1
    print(count_table)
    # output the result in json
    trial_output = []
    for trial in trial_list:
        trial_output.append(trial.to_dict())
    json.dump(trial_output, open(args.output + '.json', 'w'), indent=4)
    # output the result in tsv
    with open(args.output + '.tsv', 'w') as fout:
        fout.write("nct_id\tstatus\tlast_update_time\tphase\tstudy_type\tprimary_purpose\ttarget_category\tagent_type\tdrugs\tdrugbank_name\tdrug_target_group\tgene_symbols\tconfidence_level\n")
        for trial in trial_list:
            #
            #
            if trial.drug_objects is not None and len(trial.drug_objects) > 0:
                drugbank_name = ';'.join([d.drug_name for d in trial.drug_objects])
            else:
                drugbank_name = 'NA'
            #
            gene_symbols = ';'.join([t.gene_symbol for t in trial.targets if t.gene_symbol is not None])
            drugs = ';'.join(trial.drugs) if len(trial.drugs) > 0 else 'NA'
            # Capitalize the first letter
            # remove \x00 in drugs
            drugs = drugs.replace('\x00', '')
            trial_phase = trial.phase if trial.phase != '' else 'NA'
            trial_purpose = trial.primary_purpose if trial.primary_purpose != 'n\/a' else 'NA'
            confidence_level = 'NA'
            if hasattr(trial, 'confidence_level'):
                confidence_level = str(trial.confidence_level)
                confidence_level = re.sub(r'\(.*?\)', '', confidence_level).strip()
            fout.write(f"{trial.nct_id}\t{trial.status}\t{trial.last_update_time}\t{trial_phase}\t{trial.study_type}\t{trial_purpose}\t{trial.target_category}\t{trial.agent_type}\t{drugs}\t{drugbank_name}\t{trial.target_group if trial.target_group is not None else 'NA'}\t{gene_symbols if gene_symbols != '' else 'NA'}\t{confidence_level}\n")


def openai_get_target_gene_from_summary_react(client, trial):
    """
    The purpose of this function is to use openai to read the trial description and summarize the gene target of trial.
    """
    # use openai to summarize the trial
    #response = client.chat.completions.create(
    response = client.responses.parse(
        model = "gpt-5-mini",
        #model = "gpt-4o",
        tools = [{"type": "web_search"}],
        input = [
            {"role": "system", 
             "content": (
                 "You are an expert in clinical trial design and analysis. Your task is to summarize the target genes for the drug based on the provided description"
                 "Return the result in the format specified by `text_format`. If no target gene of drug can be inferred return 'N/A'. "
                 "Solve this using the ReAct (Reasoning and Acting) approach, structured as follows: "
                 "**Observation**: Analyze the provided drug description to identify any explicit or implied references to target genes, proteins, or associated pathways relevant to the drug’s mechanism of action."
                 "**Reason**: Reason through the information to determine the specific gene(s) targeted by the drug. Consider the drug’s mechanism (e.g., direct binding, mRNA degradation, pathway modulation) and any disease-specific context. If no gene targets are clear, evaluate whether the description suggests non-genetic targets or lacks sufficient detail."
                 "**Act**: Extract the official HGNC-approved gene names for the identified targets. If necessary, cross-reference with reliable sources (e.g., preclinical studies, clinical trial data, or pharmacological literature) to confirm targets. If no genes are identifiable, conclude with 'N/A'."
                 "**Result**: Present the target gene(s) in the format specified by text_format. If 'N/A', provide a brief justification for the lack of identifiable gene targets."
                 "Output Requirements: Ensure the response is concise, accurate, and based on the provided description or inferred from credible pharmacological knowledge."
                 "- Use HGNC-approved gene symbols (e.g., APP, TNF) for consistency."
                 "- Report confidence level from 1 to 5 with 5 the most confident."
                 "- If the description is ambiguous or lacks gene-specific information, clearly state the reasoning for returning 'N/A'."
             )},
              {"role": "user", "content": f"Extract the target gene of drug according to the following information: drug: {trial['drug']} and clinical trial title: {trial['title']}. If the drug is not known, use web search to find the drug information. Also output the web search results as explanation."}
        ],
        #temperature=0.1,
        text_format = TargetGeneExtraction,
        #response_format={"type": "json_schema", "schema": schema}
    )
    # print the response
    print("")
    print(">>>>>")
    print(response.output_parsed)
    print("<<<<<")
    print("")
    # parse the response into the TargetCategoryExtraction model
    result_parse = response.output_parsed
    # get the "parsed" field from the response
    result = {}
    result['gene_symbol'] = result_parse.gene_symbol
    result['confidence_level'] = result_parse.confidence_level
    result['explanation_target_gene'] = result_parse.explanation
    return result
    

if __name__=="__main__":
    main()

