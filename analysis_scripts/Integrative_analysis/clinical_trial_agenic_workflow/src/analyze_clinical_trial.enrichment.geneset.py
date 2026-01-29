#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 15 Oct 2025 05:19:50 PM

import os,sys,argparse
import pandas as pd
import json
from tqdm import tqdm
import random
# add parallel computing
from concurrent.futures import ProcessPoolExecutor, as_completed


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",default=1,help="Number of threads to use")
    p.add_argument('--input_drugbank_db',type=str,dest="input_drugbank_db",help="DrugBank database in json format")
    p.add_argument('--input_trial_summary',type=str,dest="input_trial_summary",help="")
    p.add_argument('--config_geneset',type=str,dest="config_geneset",help="")
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


def load_geneset_data(config_geneset):
    """
    Load gene set vs erosion or mingling score
    """
    geneset_format = config_geneset.get('geneset_format', None)
    if geneset_format is None:
        raise ValueError("geneset_format must be specified in config_geneset")
        exit(1)
    elif geneset_format == 'custom':
        geneset_data_list = []
        geneset_files = config_geneset.get('geneset_files', [])
        for geneset_file in geneset_files:
            df = pd.read_csv(geneset_file, sep="\t", header=0)
            # group by gene_set, cor_score, group, and methhod
            geneset_group = df.groupby(['gene_set', 'cor_score', 'group', 'method'])
            # for each group, get the gene list
            for group_name, group_df in geneset_group:
                gene_set_name = group_name[0]
                cor_score = group_name[1]
                group_label = group_name[2]
                method = group_name[3]
                gene_list = group_df['gene_name'].tolist()
                geneset_data_list.append({'gene_set': gene_set_name, 'cor_score': cor_score, 'group': group_label, 'method': method, 'genes': gene_list})
        # print logging
        print(f"[Logging] Loaded {len(geneset_data_list)} gene sets from {len(geneset_files)} files")
        return geneset_data_list
    else:
        print(f"[Error] Unsupported geneset_format: {geneset_format}")
        exit(1)


def main():
    global args
    args = parse_arg()
    # load the DrubBank database
    drugbank_df = json.load(open(args.input_drugbank_db))
    drug2gene = parse_drugbank(drugbank_df)
    print(f"[Logging] Total {len(drug2gene)} drugs with gene information parsed from DrugBank database")
    # load Gene set config
    config_geneset = json.load(open(args.config_geneset))
    data_geneset_list = load_geneset_data(config_geneset)
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
    with open(args.output_prefix + ".enrichment_results.txt", 'w') as fout:
        print("group_label\tgene_set\tcor_score\tmethod\tgroup\tnum_target_genes\tnum_geneset_genes\tnum_overlap_genes\tbg_idx\tbg_count\tp_value", file=fout)
        pool = ProcessPoolExecutor(max_workers=args.thread)
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
            if args.thread == 1:
                # for each geneset dataset, perform enrichment analysis
                for geneset_entry in data_geneset_list:
                    # {'gene_set': gene_set_name, 'cor_score': cor_score, 'group': group_label, 'method': method, 'genes': gene_list}
                    geneset_name = geneset_entry['gene_set']
                    cor_score = geneset_entry['cor_score']
                    group = geneset_entry['group']
                    method = geneset_entry['method']
                    geneset_gene_set = set(geneset_entry['genes'])
                    # get the number of overlapping genes between target genes and genesets
                    print(f"[Logging] Analyzing Gene set: {geneset_name}, cor_score: {cor_score}, group: {group}, method: {method} with {len(geneset_gene_set)} genes")
                    # for each geneset, calculate overlap with geneset
                    obv_overlap_target_genes = target_gene_set.intersection(geneset_gene_set)
                    # skip if observed overlap is zero
                    if len(obv_overlap_target_genes) == 0:
                        print(f"[Logging] Skipping Gene set: {geneset_name} due to zero overlap with target genes")
                        # still print into the output with zero overlap
                        print(f"{group_label}\t{geneset_name}\t{cor_score}\t{method}\t{group}\t{len(target_gene_set)}\t{len(geneset_gene_set)}\t0\tNA\tNA\t1.000000", file=fout)
                        continue
                    # prepare results storage
                    result = {'gene_set': geneset_name, 'cor_score': cor_score, 'method': method, 'group': group, 'num_target_genes': len(target_gene_set), 'num_geneset_genes': len(geneset_gene_set), 'num_overlap_genes': len(obv_overlap_target_genes), 'num_overlap_gene_bg': []}
                    # randomly sample drugs and calculate overlap with genesets to build null distribution 
                    print(f"[Logging]\t Observed overlap between target genes and genes in the geneset: {result['num_overlap_genes']} genes")
                    print(f"[Logging]\t Performing bootstrap with {bootstrap_count} iterations to calculate p-value")
                    for nn in tqdm(range(bootstrap_count)):
                        # randomly sample drugs to get a gene set of the same size as target_gene_set
                        sampled_gene_set = make_sampled_gene_set(drug2gene, len(target_gene_set))
                        # calculate overlap with genesets
                        obv_overlap_sampled_genes = sampled_gene_set.intersection(geneset_gene_set)
                        result['num_overlap_gene_bg'].append(len(obv_overlap_sampled_genes))
                    # calculate p-value, the proportion of bootstrap overlaps that are >= observed overlap
                    bg_overlaps = result['num_overlap_gene_bg']
                    p_value = sum([1 for x in bg_overlaps if x >= result['num_overlap_genes']]) / bootstrap_count
                    # print summary of results
                    print(f"[Logging]\t P-value for enrichment: {p_value:.6f}")
                    for bg_idx in range(bootstrap_count):
                        print(f"{group_label}\t{geneset_name}\t{cor_score}\t{method}\t{group}\t{result['num_target_genes']}\t{result['num_geneset_genes']}\t{result['num_overlap_genes']}\t{bg_idx}\t{bg_overlaps[bg_idx]}\t{p_value:.6f}", file=fout)
            else:
                # use parallel computing
                p_list = []
                for geneset_entry in data_geneset_list:
                    geneset_name = geneset_entry['gene_set']
                    cor_score = geneset_entry['cor_score']
                    group = geneset_entry['group']
                    method = geneset_entry['method']
                    geneset_gene_set = set(geneset_entry['genes'])
                    # for each geneset, calculate overlap with geneset
                    obv_overlap_target_genes = target_gene_set.intersection(geneset_gene_set)
                    # skip if observed overlap is zero
                    if len(obv_overlap_target_genes) == 0:
                        print(f"[Logging] Skipping Gene set: {geneset_name} due to zero overlap with target genes")
                        # still print into the output with zero overlap
                        print(f"{group_label}\t{geneset_name}\t{cor_score}\t{method}\t{group}\t{len(target_gene_set)}\t{len(geneset_gene_set)}\t0\tNA\tNA\t1.000000", file=fout)
                        continue
                    p = pool.submit(compute_geneset_enrichment, group_label, target_gene_set, geneset_name, cor_score, method, group, geneset_gene_set, drug2gene, bootstrap_count)
                    p_list.append(p)
                # collect results
                for p in tqdm(as_completed(p_list), total=len(p_list)):
                    result_lines = p.result()
                    for line in result_lines:
                        print(line, file=fout)
        pool.shutdown(wait=True)
        print("All done!")


def compute_geneset_enrichment(group_label, target_gene_set, geneset_name, cor_score, method, group, geneset_gene_set, drug2gene, bootstrap_count):
    obv_overlap_target_genes = target_gene_set.intersection(geneset_gene_set)
    # prepare results storage
    result = {'gene_set': geneset_name, 'cor_score': cor_score, 'method': method, 'group': group, 'num_target_genes': len(target_gene_set), 'num_geneset_genes': len(geneset_gene_set), 'num_overlap_genes': len(obv_overlap_target_genes), 'num_overlap_gene_bg': []}
    # randomly sample drugs and calculate overlap with genesets to build null distribution
    for nn in range(bootstrap_count):
        # randomly sample drugs to get a gene set of the same size as target_gene_set
        sampled_gene_set = make_sampled_gene_set(drug2gene, len(target_gene_set))
        # calculate overlap with genesets
        obv_overlap_sampled_genes = sampled_gene_set.intersection(geneset_gene_set)
        result['num_overlap_gene_bg'].append(len(obv_overlap_sampled_genes))
    # calculate p-value, the proportion of bootstrap overlaps that are >= observed overlap
    bg_overlaps = result['num_overlap_gene_bg']
    p_value = sum([1 for x in bg_overlaps if x >= result['num_overlap_genes']]) / bootstrap_count
    # prepare result
    return [f"{group_label}\t{geneset_name}\t{cor_score}\t{method}\t{group}\t{result['num_target_genes']}\t{result['num_geneset_genes']}\t{result['num_overlap_genes']}\t{bg_idx}\t{bg_overlaps[bg_idx]}\t{p_value:.6f}" for bg_idx in range(bootstrap_count)]


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

