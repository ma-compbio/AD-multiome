#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 12 Oct 2024 12:23:06 AM

import os,sys,argparse
import gzip
import yaml
from rich import print

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--config', dest='config', help='config YAML file', required=True)
    p.add_argument('--job_list', dest='job_list', nargs="+", help='job list, if not set will create sbatch script for all jobs in the config file', required=False)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

##################################
# file format checking
##################################

def is_valid_fastq(fastq_file):
    """
    given a fastq file check if it is a valid fastq file
    This is a simple check, only check the first 40 lines
    """
    if '.gz' in fastq_file:
        fin = gzip.open(fastq_file, 'rt')
    else:
        fin = open(fastq_file, 'r')
    # read the first 40 lines:
    for i in range(40):
        line = fin.readline().strip()
        if i % 4 == 0:
            if not line.startswith('@'):
                return False
            if len(line) < 2:
                return False
        if i % 4 == 1:
            read_len = len(line)
            if read_len < 10:
                return False
            if not set(line).issubset(set('ACGTNactgn')):
                return False
        if i % 4 == 2:
            if not line.startswith('+'):
                return False
        if i % 4 == 3:
            qual_len = len(line)
            if qual_len != read_len:
                return False
    fin.close()
    return True

##################################
# config option checking
##################################

def option_check_global(config):
    """
    check the global parameters
    """
    config_section = 'global'
    required_folder = ['working_folder', 'script_folder', 'out_job_script_folder', 'out_result_folder']
    if config_section not in config:
        print("[red]Error: global not set, please check[/red]")
        exit(1)
    for key in required_folder:
        if key not in config[config_section]:
            print(f'[red]Error: {key} not set in config file, please check[/red]')
            exit(1)
        # check if folder exists
        if not os.path.exists(config[config_section][key]):
            print(f'[red]Error: {config[config_section][key]} does not exist, creating it then re-run[/red]')
            exit(1) 


def option_check_sample(config):
    """
    check the sample parameters 
    """
    config_section = 'sample'
    # required parameters
    if config_section not in config:
        print('[red]Error: sample not set, please check[/red]')
        exit
    if 'fastq_file_R1' not in config[config_section] or 'fastq_file_R2' not in config[config_section]:
        print('[red]Error: fastq_file_R1 or fastq_file_R2 not set, please check[/red]')
        exit(1)
    # fastq file exists
    fastq_R1 = config[config_section]['fastq_file_R1']
    fastq_R2 = config[config_section]['fastq_file_R2']
    for fastq_file in [fastq_R1, fastq_R2]:
        if not os.path.exists(fastq_file):
            print(f'[red]Error: {fastq_file} does not exist, please check[/red]')
            exit(1)
    # if fastq files are valided
    for fastq_file in [fastq_R1, fastq_R2]:
        if not is_valid_fastq(fastq_file):
            print(f'[red]Error: {fastq_file} is not a valid fastq file, please check[/red]')
            exit(1) 
    # out_prefix
    if 'out_prefix' not in config[config_section]:
        print('[red]Error: out_prefix not set, please check[/red]')
        exit(1)
    # create out folder
    out_result_folder = config['global']['out_result_folder']
    sample_result_folder = os.path.join(out_result_folder, config[config_section]['out_prefix'])
    # create the folder
    if not os.path.exists(sample_result_folder):
        os.makedirs(sample_result_folder)


def option_check_slurm(config, config_section):
    if config_section not in config:
        print(f'[red]Error: {config_section} not set, please check[/red]')
        exit(1)
    if 'slurm' not in config[config_section]:
        print('[red]Error: slurm not set, please check[/red]')
        config[config_section]['slurm'] = {}
        slurm_config = config[config_section]['slurm']
    else:
        slurm_config = config[config_section]['slurm']
    # check if partition
    if 'queue' not in slurm_config:
        print('partition not set, please check')
        print('use default partition: RM-shared')
        slurm_config['queue'] = 'RM-shared'
    if 'time' not in slurm_config:
        print('time not set, please check')
        print('use default time: 2:00:00')
        slurm_config['time'] = '2:00:00'
    if 'ntasks-per-node' not in slurm_config:
        print('ntasks-per-node not set, please check')
        print('use default ntasks-per-node: 1')
        slurm_config['ntasks-per-node'] = 1 
    # RM-shared mem is 1.8G per core
    if 'mem' not in slurm_config:
        print('mem not set, please check')
        print('use default mem: %dG' % (int(slurm_config['ntasks-per-node'] * 1.8)))
        slurm_config['mem'] = '%dG' % (int(slurm_config['ntasks-per-node'] * 1.8))
    # check if mem per core is exceeded  
    if slurm_config['queue'] == 'RM-shared':
        if int(slurm_config['mem'].replace('Gb', '')) > 1.8 * slurm_config['ntasks-per-node']:
            print("[red]Warning: mem per core is exceeded, reset it to 1.8G per core or use the RM queue[/red]")
            slurm_config['mem'] = '%dG' % (int(slurm_config['ntasks-per-node'] * 1.8))


##################################
# job script building
##################################

def build_initial_report(config):
    """
    report basic summary of the config files
    """
    print("\n[bold magenta]Summary of the config file:[/bold magenta]")
    print("Creating jobs for the following tasks: " + ','.join([" [italic green]"+job+"[/italic green]" for job in args.job_list]) + '\n')
    print("Input fastq file:")
    print("\tR1: " + config['sample']['fastq_file_R1'])
    print("\tR2: " + config['sample']['fastq_file_R2'])
    print("Output folder:")
    print("\t" + config['global']['out_result_folder'])
    print("Slurm script folder:")
    print("\t" + config['global']['out_job_script_folder'] + '\n')


def build_job_fastqc(config):
    """
    Create slurm script for running Fastqc on fastq files
    """
    config_section = 'job_fastqc'
    print("Job: fastqc")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    sample_prefix = config['sample']['out_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.fastqc.sbatch.sh')
    if os.path.exists(file_sbatch_script):
        if config['global']['force_overwrite'] == False:
            print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
            return True 
        elif config['global']['force_overwrite'] == True:
            print(f'[red]Warning: {file_sbatch_script} already exists, running the script will overwrite result files[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'fastqc'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_fastqc'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# fastqc")
    # check if job script exist
    job_script = os.path.join(config['global']['script_folder'], config_job['job_script'])
    if not os.path.exists(job_script):
        print(f'[red]Error: {job_script} does not exist, please check[/red]')
        exit(1)
    # check if input fastq file exists
    input_R1 = config['sample']['fastq_file_R1']
    input_R2 = config['sample']['fastq_file_R2']
    for fastq_file in [input_R1, input_R2]:
        if not os.path.exists(fastq_file):
            print(f'[red]Error: {fastq_file} does not exist, please check[/red]')
            exit(1)
    # check if output folder exists
    output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # check if output files exist
    output_R1 =  os.path.join(output_folder, sample_prefix + '_FastQC_R1')
    output_R2 = os.path.join(output_folder, sample_prefix + '_FastQC_R2')
    for output_file in [output_R1, output_R2]:
        if not os.path.exists(output_file):
            os.makedirs(output_file)
        else:
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {output_file} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {output_file} already exists, running the script will overwrite result files[/red]')
    # add output files to the config 
    config['sample']['output_R1'] = output_R1  
    config['sample']['output_R2'] = output_R2
    sbatch_script.append(f"{job_script} -t {config_job['slurm']['ntasks-per-node']} -o {config['sample']['output_R1']} {input_R1}")
    sbatch_script.append(f"{job_script} -t {config_job['slurm']['ntasks-per-node']} -o {config['sample']['output_R2']} {input_R2}")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tOutput R1: " + output_R1)
    print("\tOutput R2: " + output_R2)
    print("") 


def build_job_script_demultiplex(config):
    """
    Create slurm script for demultiplexing fastq files
    """
    config_section = 'job_demultiplex'
    print("Job: demultiplex")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    sample_prefix = config['sample']['out_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.demultiplex.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, running the script will overwrite result files[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'demultiplex'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_demultiplex'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# helper function")
    sbatch_script.append("function trimSeqName { awk -v OFS='\\t' 'NR==1 { idx=1; for(c=0;c<4;idx++) if(substr($0,idx,1)==\":\") c++; } NR%4==1 { $0 = \"@\"substr($0,idx) } 1' ; }")
    sbatch_script.append("# demultiplexing")
    # check if job script exist
    job_script = os.path.join(config['global']['script_folder'], config_job['job_script'])
    if not os.path.exists(job_script):
        print(f'[red]Error: {job_script} does not exist, please check[/red]')
        exit(1)
    # check if barcode file exist
    barcode_file_bc1 = config_job['adaptor_bc1']
    barcode_file_bc2 = config_job['adaptor_bc2']
    barcode_file_l1 = config_job['adaptor_l1']
    barcode_file_l2 = config_job['adaptor_l2']
    if config_job['barcode_keep'] == "NA":
        barcode_keep = ""
    else:
        barcode_keep = config_job['barcode_keep']
    for filename in [barcode_file_bc1, barcode_file_bc2, barcode_file_l1, barcode_file_l2]:
        if not os.path.exists(filename):
            print(f'[red]Error: {filename} does not exist, please check[/red]')
            exit(1)
    # check if input fastq file exists
    input_R1 = config['sample']['fastq_file_R1']
    input_R2 = config['sample']['fastq_file_R2']
    for fastq_file in [input_R1, input_R2]:
        if not os.path.exists(fastq_file):
            print(f'[red]Error: {fastq_file} does not exist, please check[/red]')
            exit(1)
    # check if output folder exists
    output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # check if output files exist
    output_R1 =  os.path.join(output_folder, sample_prefix + '_demulti_R1.fastq')
    output_R2 = os.path.join(output_folder, sample_prefix + '_demulti_R2.fastq')
    for output_file in [output_R1, output_R2]:
        if os.path.exists(output_file):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {output_file} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {output_file} already exists, running the script will overwrite result files[/red]')
    # add output files to the config 
    config['sample']['demulti_R1'] = output_R1  
    config['sample']['demulti_R2'] = output_R2
    sbatch_script.append(f"python {job_script} --barcode={barcode_file_bc2} --pos1=\"[]\" --pos2=\"slice(0,8)\"  --mode=1 --barcode={barcode_file_l2}  --pos1=\"[]\" --pos2=\"slice(8,23)\"  --mode=5 --barcode={barcode_file_bc1} --pos1=\"[]\" --pos2=\"slice(23,31)\" --mode=1  --infile1  <(gzip -cd \"{input_R1}\" | trimSeqName) --infile2  <(gzip -cd \"{input_R2}\" | trimSeqName) --outfile1 >(grep -a . > \"{output_R1}\") --outfile2 >(grep -a . > \"{output_R2}\") --barcode2keep=\"{barcode_keep}\"")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tOutput R1: " + output_R1)
    print("\tOutput R2: " + output_R2)
    print("")


def build_job_script_alignment(config):
    """
    Create slurm script for aligning fastq file
    """
    config_section = 'job_alignment'
    print("Job: alignment")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    sample_prefix = config['sample']['out_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if job script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.alignment.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # for alignment the minimum thread is 36
    if config_job['slurm']['ntasks-per-node'] < 36:
        print("[red]Warning: alignment thread is less than 36, reset to 36[/red]")
        config_job['slurm']['ntasks-per-node'] = 36
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'alignment'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_alignment'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# alignment")
    # check if bwa index files exist
    bwa_genome_index = config_job['genome_index']
    # .amb is text file, to record appearance of N (or other non-ATGC) in the ref fasta.
    # .ann is text file, to record ref sequences, name, length, etc.
    # .bwt is binary, the Burrows-Wheeler transformed sequence.
    # .pac is binary, packaged sequence (four base pairs encode one byte).
    # .sa is binary, suffix array index.
    for index_suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
        filename = bwa_genome_index + index_suffix
        if not os.path.exists(filename):
            print(f'[red]Error: Index file {filename} does not exist, please check[/red]')
            exit(1)
    # check if output folder exists
    output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # check if input files exist
    if 'demulti_R1' not in config['sample'] or 'demulti_R2' not in config['sample']:
        # first check if demultiplexed files exist
        input_R1 =  os.path.join(output_folder, sample_prefix + '_demulti_R1.fastq')
        input_R2 = os.path.join(output_folder, sample_prefix + '_demulti_R2.fastq')
    else:
        input_R1 = config['sample']['demulti_R1']
        input_R2 = config['sample']['demulti_R2']
    for filename in [input_R1, input_R2]:
        if not os.path.exists(filename) and config['global']['ignore_input_check'] == False:
            print(f'[red]Error: Input demultiplexed file {filename} does not exist, please check[/red]')
            exit(1)
    # check if output files exist
    output_file = os.path.join(output_folder, sample_prefix + '_alignment.bam')
    if os.path.exists(output_file):
        if config['global']['force_overwrite'] == False:
            print(f'[red]Warning: {output_file} already exists, skip job[/red]')
            return True 
        elif config['global']['force_overwrite'] == True:
            print(f'[red]Warning: {output_file} already exists, overwriting it[/red]')
        else:
            pass
    alignment_thread = config_job['slurm']['ntasks-per-node'] - 16
    samtools_thread = 16
    sbatch_script.append("bwa mem -SP5M -t%d \"%s\" <(awk 'NR%%2==0{print substr($0,13)} NR%%2==1{print}' \"%s\") <(awk 'NR%%2==0{print substr($0,36)} NR%%2==1{print}' \"%s\") | awk -v OFS='\\t' '$0~/^@/ { print; next } { split($1,bc,\"_\"); $1 = bc[1]; print $0 \"\\tBC:Z:\" bc[2] \",\" bc[4] }' | samtools view -hb@%d > \"%s\"" % (alignment_thread, bwa_genome_index, input_R1, input_R2, samtools_thread, output_file))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # add output files to the config 
    config['sample']['alignment'] = output_file
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput R1: " + input_R1)
    print("\tInput R2: " + input_R2)
    print("\tOutput file: " + output_file)
    print("")


def build_job_bam2pair(config):
    """
    convert bam to parsed pairs
    """
    config_section = 'job_bam2pair'
    print("Job: parse bam aligment to pairs")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    sample_prefix = config['sample']['out_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.bam2pair.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # check if job script exist
    job_script = os.path.join(config['global']['script_folder'], config_job['job_script'])
    if not os.path.exists(job_script):
        print(f'[red]Error: {job_script} does not exist, please check[/red]')
        exit(1)
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'bam2pair'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_bam2pair'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# bam2pair")
    # check if rescue_walk.py exist
    if 'script_rescue_walk' not in config_job or not os.path.exists(config_job['script_rescue_walk']):
        print(f'[red]Error: {config_job["script_rescue_walk"]} does not exist, please check[/red]')
        exit(1)
    file_rescue_walk = config_job['script_rescue_walk'] 
    # check if genome chromosome size file exist
    if 'genome_size' not in config_job or not os.path.exists(config_job['genome_size']):
        print(f'[red]Error: {config_job["genome_size"]} does not exist, please check[/red]')
        exit(1)
    # check if walk policy parameter exist
    if 'walk_policy' not in config_job:
        print(f'[red]Waning: walk_policy not set, use `complete` as default[/red]')
        config_job['walk_policy'] = 'complete'
    file_genome_size = config_job['genome_size']
    # check if output folder exists
    output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # check if input files exist
    if 'alignment' not in config['sample']:
        # first check if alignment files exist
        input_bam =  os.path.join(output_folder, sample_prefix + '_alignment.bam')
    else:
        input_bam = config['sample']['alignment']
    for filename in [input_bam]:
        if not os.path.exists(filename) and config['global']['ignore_input_check'] == False:
            print(f'[red]Error: Input alignment file {filename} does not exist, please check[/red]')
            exit(1)
    # check if output files exist
    walk_policy = config_job['walk_policy']
    if ',' in walk_policy: 
        walk_policy_list = walk_policy.split(',')
    else:
        walk_policy_list = [walk_policy]
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput file: " + input_bam)
    for walk_policy in walk_policy_list:
        output_file = os.path.join(output_folder, walk_policy, sample_prefix + '.pairs.gz')
        os.system("mkdir -p %s" % (os.path.join(output_folder, walk_policy)))
        if os.path.exists(output_file):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {output_file} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {output_file} already exists, overwriting it[/red]')
            else:
                pass
        #
        thread = config_job['slurm']['ntasks-per-node']
        genome_size = config_job['genome_size']
        #
        if 'complete' == walk_policy:
            sbatch_script.append(f"bash {job_script} -t {thread} -w {walk_policy} -b {input_bam} -r {file_rescue_walk} -c {file_genome_size} -p {output_file}")
        else:
            sbatch_script.append(f"bash {job_script} -t {thread} -w {walk_policy} -b {input_bam} -r {file_rescue_walk} -c {file_genome_size} -p {output_file}")
        print("\tOutput file: " + output_file)
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # add output files to the config 
    config['sample']['alignment'] = output_file
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("")


def build_job_script_merge_pair(config):
    """
    merge pairs 
    """
    config_section = 'job_merge_pair'
    print("Job: merge pairs from the same library")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    if 'output_prefix' not in config_job:
        print(f'[red]Error: output_prefix not set in merge_pair config, please check[/red]')
        exit(1)
    sample_prefix = config_job['output_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.merge_pair.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'merge_pair'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_merge_pair'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# merge pair")
    # check walk policy parameter exist
    if 'walk_policy' not in config_job:
        print(f'[red]Waning: walk_policy not set, use `complete` as default[/red]')
        config_job['walk_policy'] = 'complete'
    walk_policy_list = [] 
    if ',' in config_job['walk_policy']:
        walk_policy_list = config_job['walk_policy'].split(',')
    else:
        walk_policy_list = [config_job['walk_policy']]
    # for each walk policy
    output_file_list = []
    for walk_policy in walk_policy_list:
        # check if output folder exists
        output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix, walk_policy)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # check if input files exist
        if 'input_pair_list' not in config_job:
            print(f'[red]Error: input_pair_list not set in config file, please check[/red]')
            exit(1)
        else:
            input_pair_list = [filename.replace('<walk_policy>', walk_policy) for filename in config_job['input_pair_list']]
        for filename in input_pair_list:
            if not os.path.exists(filename):# and config['global']['ignore_input_check'] == False:
                print(f'[red]Error: Input pair file {filename} does not exist, please check[/red]')
                exit(1)
        # check if output files exist
        output_file = os.path.join(output_folder, sample_prefix + '.pairs.gz')
        if os.path.exists(output_file):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {output_file} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {output_file} already exists, overwriting it[/red]')
            else:
                pass
        #
        thread = config_job['slurm']['ntasks-per-node']
        #
        sbatch_script.append("pigz -cdp%d %s | pigz -6p%d > %s" % (2, ' '.join(input_pair_list), thread, output_file))
        # add output files to the config 
        config['sample']['merged_pairs'] = "%s:%s" % (walk_policy, output_file)
        output_file_list.append(output_file)
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput file: " + ','.join(input_pair_list))
    print("\tOutput file: " + ','.join(output_file_list))
    print("")


def build_job_script_dedup(config): 
    """
    merge pairs 
    """
    config_section = 'job_dedup'
    print("Job: remove PCR duplicates from pairs")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    if 'output_prefix' not in config_job:
        print(f'[red]Error: output_prefix not set in merge_pair config, please check[/red]')
        exit(1)
    sample_prefix = config_job['output_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.dedup_pair.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'dedup_pair'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_dedup_pair'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# dedup pair")
    # check walk policy parameter exist
    if 'walk_policy' not in config_job:
        print(f'[red]Waning: walk_policy not set, use `complete` as default[/red]')
        config_job['walk_policy'] = 'complete'
    walk_policy_list = [] 
    if ',' in config_job['walk_policy']:
        walk_policy_list = config_job['walk_policy'].split(',')
    else:
        walk_policy_list = [config_job['walk_policy']]
    output_file_dedup_list = []
    output_file_dedup_wo1k_list = []
    for walk_policy in walk_policy_list:
        # check if output folder exists
        output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix, walk_policy)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # check if input files exist
        if 'input_pair' not in config_job:
            print(f'[red]Error: input_pair not set in config file, please check[/red]')
            exit(1)
        else:
            input_pair = config_job['input_pair'].replace('<walk_policy>', walk_policy)
        for filename in [input_pair]:
            if not os.path.exists(filename) and config['global']['ignore_input_check'] == False:
                print(f'[red]Error: Input pair file {filename} does not exist, please check[/red]')
                exit(1)
        # check if output files exist
        output_file_prefix = os.path.join(output_folder, sample_prefix)
        output_file_dedup = output_file_prefix + '_nodup.pairs.gz' 
        output_file_dedup_wo1k = output_file_prefix + '_nodup_wo1k.pairs.gz'
        for filename in [output_file_dedup, output_file_dedup_wo1k]:
            if os.path.exists(filename):
                if config['global']['force_overwrite'] == False:
                    print(f'[red]Warning: {filename} already exists, skip job[/red]')
                    return True 
                elif config['global']['force_overwrite'] == True:
                    print(f'[red]Warning: {filename} already exists, overwriting it[/red]')
                else:
                    pass
                pass
        # check dedup_CPP program
        if 'script_dedup_CPP' not in config_job or not os.path.exists(config_job['script_dedup_CPP']):
            print(f'[red]Error: {config_job["script_dedup_CPP"]} does not exist, please check[/red]')
            exit(1)
        script_dedup = config_job['script_dedup_CPP']
        # check if job script exist
        job_script = os.path.join(config['global']['script_folder'], config_job['job_script'])
        if not os.path.exists(job_script):
            print(f'[red]Error: {job_script} does not exist, please check[/red]')
            exit(1)
        #
        thread = config_job['slurm']['ntasks-per-node']
        #
        sbatch_script.append("bash %s -t %d -s %s -p %s -o %s" % (job_script, thread, script_dedup, input_pair, output_file_prefix))
        config['sample']['merged_dedup_pair'] = "%s:%s" % (walk_policy, output_file_dedup)
        config['sample']['merged_dedup_pair_wo1k'] = "%s:%s" % (walk_policy, output_file_dedup_wo1k)
        output_file_dedup_list.append(output_file_dedup)
        output_file_dedup_wo1k_list.append(output_file_dedup_wo1k)
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # add output files to the config 
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput file: " + input_pair)
    print("\tOutput file (dedeup): " + ','.join(output_file_dedup_list))
    print("\tOutput file (dedeup wo1k): " + ','.join(output_file_dedup_wo1k_list))
    print("")


def build_job_script_qc_bam_file(config):
    """
    Conduct QC analysis of bam file
    """
    config_section = 'job_qc_bam'
    print("Job: perform quality control check for bam file")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    if 'output_prefix' not in config_job:
        print(f'[red]Error: output_prefix not set in merge_pair config, please check[/red]')
        exit(1)
    sample_prefix = config_job['output_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.qc_bam.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'bam2pair'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_bam2pair'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# bam2pair")
    # check if rescue_walk.py exist
    if 'script_rescue_walk' not in config_job or not os.path.exists(config_job['script_rescue_walk']):
        print(f'[red]Error: {config_job["script_rescue_walk"]} does not exist, please check[/red]')
        exit(1)
    file_rescue_walk = config_job['script_rescue_walk'] 
    # check if genome chromosome size file exist
    if 'genome_size' not in config_job or not os.path.exists(config_job['genome_size']):
        print(f'[red]Error: {config_job["genome_size"]} does not exist, please check[/red]')
        exit(1)
    # check if walk policy parameter exist
    if 'walk_policy' not in config_job:
        print(f'[red]Waning: walk_policy not set, use `complete` as default[/red]')
        config_job['walk_policy'] = 'complete'
    file_genome_size = config_job['genome_size']
    # check if output folder exists
    output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # check if input files exist
    if 'alignment' not in config['sample']:
        # first check if alignment files exist
        input_bam =  os.path.join(output_folder, sample_prefix + '_alignment.bam')
    else:
        input_bam = config['sample']['alignment']
    for filename in [input_bam]:
        if not os.path.exists(filename) and config['global']['ignore_input_check'] == False:
            print(f'[red]Error: Input alignment file {filename} does not exist, please check[/red]')
            exit(1)
    # check if output files exist
    walk_policy = config_job['walk_policy']
    if ',' in walk_policy: 
        walk_policy_list = walk_policy.split(',')
    else:
        walk_policy_list = [walk_policy]
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput file: " + input_bam)
    for walk_policy in walk_policy_list:
        output_file = os.path.join(output_folder, walk_policy, sample_prefix + '.pairs.gz')
        os.system("mkdir -p %s" % (os.path.join(output_folder, walk_policy)))
        if os.path.exists(output_file):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {output_file} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {output_file} already exists, overwriting it[/red]')
            else:
                pass
        #
        thread = config_job['slurm']['ntasks-per-node']
        genome_size = config_job['genome_size']
        #
        if 'complete' == walk_policy:
            sbatch_script.append(f"bash {job_script} -t {thread} -w {walk_policy} -b {input_bam} -r {file_rescue_walk} -c {file_genome_size} -p {output_file}")
        else:
            sbatch_script.append(f"bash {job_script} -t {thread} -w {walk_policy} -b {input_bam} -r {file_rescue_walk} -c {file_genome_size} -p {output_file}")
        print("\tOutput file: " + output_file)
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("")  
    sbatch_script.append("")
    # add output files to the config 
    config['sample']['alignment'] = output_file
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("")


def build_job_script_qc_pair_file(config):
    """
    Conduct QC analysis of pair file
    """
    config_section = 'job_qc_pair'
    print("Job: perform quality control check for pairs")
    # check if slurm parameters valid
    option_check_slurm(config, config_section)
    # config job
    config_job = config[config_section]
    # output prefix
    if 'output_prefix' not in config_job:
        print(f'[red]Error: output_prefix not set in merge_pair config, please check[/red]')
        exit(1)
    sample_prefix = config_job['output_prefix']
    # set the location of the sbatch script
    sample_folder = os.path.join(config['global']['out_job_script_folder'], sample_prefix)
    if not os.path.exists(sample_folder):
        os.makedirs(sample_folder)
    # check if sbatch script exist
    file_sbatch_script = os.path.join(sample_folder, sample_prefix  + '.qc_pair.sbatch.sh')
    if os.path.exists(file_sbatch_script):
            if config['global']['force_overwrite'] == False:
                print(f'[red]Warning: {file_sbatch_script} already exists, skip job[/red]')
                return True 
            elif config['global']['force_overwrite'] == True:
                print(f'[red]Warning: {file_sbatch_script} already exists, overwriting it[/red]') 
    # sbatch script start
    sbatch_script = ["#!/bin/bash"]
    sbatch_script.append("")
    sbatch_script.append("#SBATCH -p {}".format(config_job['slurm']['queue']))
    sbatch_script.append("#SBATCH --time={}".format(config_job['slurm']['time']))
    sbatch_script.append("#SBATCH --ntasks-per-node={}".format(config_job['slurm']['ntasks-per-node']))
    sbatch_script.append("#SBATCH --mem={}".format(config_job['slurm']['mem']))
    # set job name and log file
    if 'name' not in config_job:
        config_job['name'] = 'qc_pair'
    sbatch_script.append("#SBATCH --job-name={}".format(config_job['name']))
    if 'log' not in config_job:
        config_job['log'] = 'log_qc_pair'
    # set the location of the log file
    file_log = os.path.join(sample_folder, sample_prefix + '.' + config_job['log'])
    sbatch_script.append("#SBATCH --output={}-%J".format(file_log))
    sbatch_script.append("")
    if 'conda' in config and 'env_name' in config['conda']:
        sbatch_script.append("source activate {}".format(config['conda']['env_name']))
    sbatch_script.append("cd {}".format(config['global']['working_folder']))
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("# run the script")
    # actual script to run
    sbatch_script.append("# QC pair")
    # check walk policy parameter exist
    if 'walk_policy' not in config_job:
        print(f'[red]Waning: walk_policy not set, use `complete` as default[/red]')
        config_job['walk_policy'] = 'complete'
    walk_policy_list = [] 
    if ',' in config_job['walk_policy']:
        walk_policy_list = config_job['walk_policy'].split(',')
    else:
        walk_policy_list = [config_job['walk_policy']]
    input_pair_list = []
    output_file_list = [] 
    output_file_dedup_list = []
    for walk_policy in walk_policy_list:
        # check if input files exist
        if 'input_pair_prefix' not in config_job:
            print(f'[red]Error: input_pair not set in config file, please check[/red]')
            exit(1)
        input_pair_prefix = config_job['input_pair_prefix'].replace('<walk_policy>', walk_policy)
        input_pair = input_pair_prefix + '.pairs.gz'
        input_pair_dedup = input_pair_prefix + '_nodup.pairs.gz'
        input_pair_wo1k = input_pair_prefix + '_nodup_wo1k.pairs.gz'
        input_pair_list.append(input_pair)
        for filename in [input_pair, input_pair_dedup, input_pair_wo1k]:
            if not os.path.exists(filename) and config['global']['ignore_input_check'] == False:
                print(f'[red]Error: Input pair file {filename} does not exist, please check[/red]')
                exit(1)
        # check if output folder exists
        output_folder = os.path.join(config['global']['out_result_folder'], sample_prefix, walk_policy)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # check if sample label exist
        if 'sample_label' not in config_job:
            config_job['sample_label'] = sample_prefix
        # check if output files exist
        output_file_prefix = os.path.join(output_folder, sample_prefix)
        output_prefix_pair = output_file_prefix
        output_prefix_pair_dedup = output_file_prefix + '_nodup'
        for file_prefix in [output_prefix_pair, output_prefix_pair_dedup]:
            for filename in [file_prefix + '.qc.pair_type.tsv', file_prefix + '.qc.pair_range.tsv']:
                if os.path.exists(filename):
                    if config['global']['force_overwrite'] == False:
                        print(f'[red]Warning: {filename} already exists, skip job[/red]')
                        return True 
                    elif config['global']['force_overwrite'] == True:
                        print(f'[red]Warning: {filename} already exists, overwriting it[/red]')
                    else:
                        pass
        # check if job script exist
        job_script = os.path.join(config['global']['script_folder'], config_job['job_script'])
        if not os.path.exists(job_script):
            print(f'[red]Error: {job_script} does not exist, please check[/red]')
            exit(1)
        #
        thread = config_job['slurm']['ntasks-per-node']
        #
        # check pair type on the original pair file
        sbatch_script.append("python %s --pair_file %s --mode %s --label %s --output %s" % (job_script, input_pair, 'pair_type pair_range', config_job['sample_label'], output_prefix_pair))
        # chck pair range on dedup pair file
        sbatch_script.append("python %s --pair_file %s --mode %s --label %s --output %s" % (job_script, input_pair_dedup, 'pair_type pair_range dup_level', config_job['sample_label'], output_prefix_pair_dedup))
        output_file_list.append(output_prefix_pair)
        output_file_dedup_list.append(output_prefix_pair_dedup)
    sbatch_script.append("")
    sbatch_script.append("")
    sbatch_script.append("") 
    sbatch_script.append("")
    # write script
    with open(file_sbatch_script, 'w') as f:
        f.write('\n'.join(sbatch_script))
    print("\tSlurm script: " + file_sbatch_script)
    print("\tInput file: " + ','.join(input_pair_list))
    print("\tOutput file (original pair): " + ','.join(output_file_list))
    print("\tOutput file (dedeup wo1k): " + ','.join(output_file_dedup_list))
    print("")


def build_final_report(config):
    return True 


##################################
# Main Function
##################################


def main():
    global args
    args = parse_arg()
    # load the YAML file
    with open(args.config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    # check global parameters
    option_check_global(config)
    # create job script
    if args.job_list is None:
        args.job_list = ['fastqc', 'demultiplex', 'alignment', 'bam2pair', 'merge_pair', 'dedup']
    else:
        pass
    if 'fastqc' in args.job_list:
        option_check_sample(config)
    if 'demultiplex' in args.job_list:
        option_check_sample(config)
    # report basic summary of the config files
    if 'demultiplex' in args.job_list:
        build_initial_report(config)
    # job 0:  fastqc
    if 'fastqc' in args.job_list:
        build_job_fastqc(config)
    # job 1: demultiplex
    if 'demultiplex' in args.job_list:
        build_job_script_demultiplex(config)
    # job 2: alignment
    if 'alignment' in args.job_list:
        build_job_script_alignment(config)
    # job 3: parse bam to pairs
    if 'bam2pair' in args.job_list:
        build_job_bam2pair(config)
    # job 4: merge pairs
    if 'merge_pair' in args.job_list:
        build_job_script_merge_pair(config)
    # job 5: dedup
    if 'dedup' in args.job_list:
        build_job_script_dedup(config)
    # job 6: qc pair
    if 'qc_pair' in args.job_list:
        build_job_script_qc_pair_file(config)
    # job 7: qc bam
    if 'qc_bam' in args.job_list:
        build_job_script_qc_bam_file(config)
    exit(1)
    

if __name__=="__main__":
    main()

