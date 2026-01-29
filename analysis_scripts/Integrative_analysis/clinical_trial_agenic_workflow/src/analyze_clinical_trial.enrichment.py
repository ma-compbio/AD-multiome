#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 15 Oct 2025 05:19:50 PM

import os,sys,argparse
import pandas as pd
import json
from tqdm import tqdm
import random


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_drugbank_db',type=str,dest="input_drugbank_db",help="DrugBank database in json format")
    p.add_argument('--input_trial_summary',type=str,dest="input_trial_summary",help="")
    p.add_argument('--config_deg',type=str,dest="config_deg",help="")
    p.add_argument('--config_trial',type=str,dest="config_trial",help="")
    p.add_argument('--output_prefix',type=str,dest="output_prefix",help="")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def parse_drugbank(drugbank_df):
    """
    Parse the DrugBank database to extract drug-gene interactions.
    return a dictionary: drug2gene
    """
    drug2gene = {}
    for drug_id in tqdm(drugbank_df):
        drug = drugbank_df[drug_id]
        drug_name = drug.get('drug_name', 'NA')
        if drug_name == 'NA':
            continue 
        else:
            drug2gene[drug_id] = {'drug_name': drug_name, 'genes': set()}
        target_list = drug.get('targets', [])
        for target in target_list:
            gene_name = target.get('gene_symbol', 'NA')
            if gene_name not in ['NA', None]:
                drug2gene[drug_id]['genes'].add(gene_name.upper())
    return drug2gene


def load_deg_data(config_deg):
    """
    Load DEG data based on the config file into a list of dataframes.
    """
    # get the deg format
    deg_format = config_deg.get('deg_format', None)
    deg_pvalue_cutoff = config_deg.get('pvalue_cutoff', 0.05)
    deg_fc_cutoff = config_deg.get('fc_cutoff', None)
    deg_minimum_genes = config_deg.get('minimum_genes', 50)
    if deg_format is None:
        raise ValueError("deg_format must be specified in config_deg")
        exit(1)
    elif deg_format == 'MAST':
        deg_data_list = []
        deg_files = config_deg.get('deg_files', [])
        for deg_file in deg_files:
            deg_df = pd.read_csv(deg_file, sep=",", header=0, index_col=0)
            cell_type = deg_file.split('/')[-1].replace('-AD-vs-CT-MAST.csv', '').replace(' ', '_')
            # rename the columns
            deg_df['gene_symbol'] = deg_df.index
            # make gene symbol upper case
            deg_df['gene_symbol'] = deg_df['gene_symbol'].str.upper()
            deg_df['p_val_adj'] = deg_df['p_val_adj']
            deg_df['log2_fc'] = deg_df['avg_log2FC']
            # filter DEGs based on p_val_adj and log2_fc thresholds defined in config_deg
            if deg_fc_cutoff is None: # default
                deg_filtered = deg_df[(deg_df['p_val_adj'] < deg_pvalue_cutoff)]
            else:
                deg_filtered = deg_df[(deg_df['p_val_adj'] < deg_pvalue_cutoff) & (deg_df['log2_fc'].abs() > deg_fc_cutoff)]
            # skip this dataset if the number of DEGs is less than minimum_genes
            if deg_filtered.shape[0] < deg_minimum_genes:
                print(f"[Logging] Skipping DEG dataset for cell type {cell_type} due to insufficient DEGs ({deg_filtered.shape[0]} < {deg_minimum_genes})")
                continue
            # if maximum number of genes is specified, keep only the top N DEGs based on p_val_adj from smallest to largest
            if 'maximum_genes' in config_deg:
                maximum_genes = config_deg['maximum_genes']
                deg_filtered = deg_filtered.sort_values(by='p_val_adj', ).head(maximum_genes)
            # select relevant columns
            deg_final = deg_filtered[['gene_symbol', 'p_val_adj', 'log2_fc']].copy()
            deg_data_list.append({"label": cell_type, "data": deg_final})
            print(f"[Logging] Loaded DEG data for cell type {cell_type} with {deg_final.shape[0]} significant DEGs")
        return deg_data_list
    else:
        print(f"[Error] Unsupported deg_format: {deg_format}")
        exit(1)


def main():
    global args
    args = parse_arg()
    # load the DrubBank database
    drugbank_df = json.load(open(args.input_drugbank_db))
    drug2gene = parse_drugbank(drugbank_df)
    print(f"[Logging] Total {len(drug2gene)} drugs with gene information parsed from DrugBank database")
    # load DEG config
    config_deg = json.load(open(args.config_deg))
    data_deg_list = load_deg_data(config_deg)
    # load trial summary
    trial_summary_df = pd.read_csv(args.input_trial_summary, sep="\t", header=0)
    # load trial config
    config_trial = json.load(open(args.config_trial))
    # group trials by the columns specific in config_trial
    groupby_cols = config_trial.get('groupby_cols', [])
    # remove the space in groupby_cols in the dataframe
    for col in groupby_cols:
        if col in trial_summary_df.columns:
            trial_summary_df[col] = trial_summary_df[col].str.replace(' ', '_')
        else:
            print(f"[Error] Column {col} not found in trial summary dataframe")
            exit(1)
    trial_groups = trial_summary_df.groupby(groupby_cols)
    # set random seed for reproducibility
    random.seed(config_trial.get('random_seed', 202510))
    bootstrap_count = config_trial.get('bootstrap_count', 1000)
    # for each trial group
    result_gene_overlap = {}
    with open(args.output_prefix + ".enrichment_results.txt", 'w') as fout:
        print("group_label\tcell_type\tnum_target_genes\tnum_deg_genes\tnum_overlap_genes\tbg_idx\tbg_count\tp_value", file=fout)
        for group_name, group_df in trial_groups:
            if isinstance(group_name, str):
                group_label = group_name.replace(' ', ';')
            else:
                group_label = ";".join([str(x).replace(' ', '_') for x in group_name])
            print("")
            print(f"[Logging] Processing trial group: {group_label} with {group_df.shape[0]} trials") # get the target genes in this group
            target_gene_set = group_df['target_gene'].dropna().tolist()
            target_gene_set = set([gene.upper() for gene in target_gene_set])
            print(f"[Logging] Found {len(target_gene_set)} unique target genes in this trial group")
            # skip this group if the number of target genes is less than parameter defined in config_trial
            if len(target_gene_set) < config_trial.get('min_target_genes', 50):
                print(f"[Logging] Skipping trial group {group_label} due to insufficient target genes")
                continue
            # for each DEG dataset, perform enrichment analysis
            for deg_entry in data_deg_list:
                cell_type = deg_entry['label']
                deg_df = deg_entry['data']
                # get the number of overlapping genes between target genes and DEGs
                deg_gene_set = set(deg_df['gene_symbol'].tolist())
                print(f"[Logging]\t DEG dataset {cell_type} has {len(deg_gene_set)} significant DEGs after filtering")
                # for each deg, calculate overlap with DEGs
                obv_overlap_target_genes = target_gene_set.intersection(deg_gene_set)
                result_gene_overlap[(group_label, cell_type)] = deg_df[deg_df['gene_symbol'].isin(obv_overlap_target_genes)].copy() 
                # prepare results storage
                result = {'cell_type': cell_type, 'group_label': group_label, 'num_target_genes': len(target_gene_set), 'num_deg_genes': len(deg_gene_set), 'num_overlap_genes': len(obv_overlap_target_genes), 'num_overlap_gene_bg': []}
                # randomly sample drugs and calculate overlap with DEGs to build null distribution 
                print(f"[Logging]\t Observed overlap between target genes and DEGs: {result['num_overlap_genes']} genes")
                print(f"[Logging]\t Performing bootstrap with {bootstrap_count} iterations to calculate p-value")
                for nn in tqdm(range(bootstrap_count)):
                    # randomly sample drugs to get a gene set of the same size as target_gene_set
                    sampled_gene_set = make_sampled_gene_set(drug2gene, len(target_gene_set))
                    # calculate overlap with DEGs
                    obv_overlap_sampled_genes = sampled_gene_set.intersection(deg_gene_set)
                    result['num_overlap_gene_bg'].append(len(obv_overlap_sampled_genes))
                # calculate p-value, the proportion of bootstrap overlaps that are >= observed overlap
                bg_overlaps = result['num_overlap_gene_bg']
                p_value = sum([1 for x in bg_overlaps if x >= result['num_overlap_genes']]) / bootstrap_count
                for bg_idx in range(bootstrap_count):
                    print(f"{group_label}\t{cell_type}\t{result['num_target_genes']}\t{result['num_deg_genes']}\t{result['num_overlap_genes']}\t{bg_idx}\t{bg_overlaps[bg_idx]}\t{p_value:.6f}", file=fout) 
    # report the overlapping genes for each (group_label, cell_type)
    with open(args.output_prefix + ".overlapping_genes.txt", 'w') as fout:
        print("group_label\tcell_type\tgene_symbol\tp_val_adj\tlog2_fc", file=fout)
        for (group_label, cell_type), overlap_df in result_gene_overlap.items():
            for idx, row in overlap_df.iterrows():
                print(f"{group_label}\t{cell_type}\t{row['gene_symbol']}\t{row['p_val_adj']}\t{row['log2_fc']}", file=fout)


def make_sampled_gene_set(drug2gene, number_of_gene):
    """
    sample drugs randomly to get a gene set of the same size as number_of_gene
    """
    drug_id_list = list(drug2gene.keys())
    # shuffle the drug_id_list
    random.shuffle(drug_id_list)
    gene_set = set()
    while True:
        # pick random drugs
        for drug_id in drug_id_list:
            gene_set.update(drug2gene[drug_id]['genes'])
        if len(gene_set) >= number_of_gene:
            break
    return set(random.sample(list(gene_set), number_of_gene))
    

if __name__=="__main__":
    main()

