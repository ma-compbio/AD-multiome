#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--config',type=str,help="config file",required=True)
    p.add_argument('--out_folder',type=str,help="output folder",required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def parse_config(config_file):
    table = {}
    with open(config_file) as fin:
        for line in fin:
            if line.strip() == "" or line.strip().startswith("#"):
                continue
            row = line.strip().split()
            key = row[0]
            value = row[1]
            table[key] = value
    return table


def main():
    global args
    args = parse_arg()
    # load config file
    config = parse_config(args.config)
    # check parameters
    for required_key in ['bam_list', 'result_folder', 'library_id']:
        if required_key not in config:
            print(f"Error: {required_key} is required in config file")
            exit(1)
    # build the sbatch script
    config['out_prefix'] = f"hPFC-scRNA_merge_{config['library_id']}"
    # merged library folder
    out_script_folder = os.path.join(args.out_folder, config['out_prefix'])
    if not os.path.exists(out_script_folder):
        os.makedirs(out_script_folder)
    #
    out_script_file = os.path.join(out_script_folder, f"{config['library_id']}_featurecount.sbatch.sh")
    out_log_file = os.path.join(out_script_folder, f"{config['library_id']}.log_featurecount")
    result_folder = os.path.join(config['result_folder'], config['out_prefix'])
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    # get the merged nodup bam file
    merged_nodup_bam_file = os.path.join(result_folder, f"{config['library_id']}_human_hg38_nodup.bam")
    if not os.path.exists(merged_nodup_bam_file):
        print(f"Error: {merged_nodup_bam_file} not exists") 
        exit(1)
    # get the gene annotation file
    gene_gtf="/jet/home/yzhang38/bighive/shared/AD_project/pipeline_RNA/annotation/gencode.v47.basic.annotation.gtf"
    with open(out_script_file, 'w') as fout:
        print("#!/bin/bash", file=fout)
        print("", file=fout)
        print("#SBATCH -p RM-shared", file = fout)
        print("#SBATCH --time=2-00:00:00", file = fout)
        config['thread'] = 16
        print("#SBATCH --ntasks-per-node=%s" % (config['thread']), file = fout) 
        print("##SBATCH --mem=40Gb", file = fout)
        print("#SBATCH --job-name feature_count", file = fout)
        print("#SBATCH --output " + out_log_file + "_%J.txt", file = fout)
        print("", file = fout)
        print("# activate conda env", file = fout) 
        print("source activate ad_project", file = fout)
        print("# go to working folder", file = fout)
        print("cd /ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/", file = fout)
        print("", file = fout)
        print("# run featureCounts and report stat", file = fout)
        print(f"/ocean/projects/bio240015p/shared/AD_project/pipeline_RNA/pipeline_script/RNA_feature_count.sh -t {config['thread']} -b {merged_nodup_bam_file} -g {gene_gtf} -l {config['library_id']} -o {result_folder}", file = fout)
        print("", file = fout)
        #
        print('sbatch ' + out_script_file)
    
if __name__=="__main__":
    main()

