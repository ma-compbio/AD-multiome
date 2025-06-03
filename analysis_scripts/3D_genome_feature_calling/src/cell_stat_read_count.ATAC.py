#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 13 Mar 2025 01:37:02 AM

import os,sys,argparse
import pandas as pd
from tqdm import tqdm
import numpy as np
import gzip


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument("--folder", dest="folder", type=str, default=None, help="The folder containing the cell read pairs files")
    p.add_argument("--patient_id", dest="patient_id", type=str, default=None, help="The patient id")
    p.add_argument('--cell_white_list',type=str,dest="cell_white_list",help="A file contains the the white list of cells by cell id")
    p.add_argument('--cell_white_list_format',type=str,dest="cell_white_list_format",default="csv",choices=['cell_label','cell_anno'],help="The format of the cell white list file")
    p.add_argument('--ATAC_folder',type=str,dest="ATAC_folder",help="folder containing the ATAC peak file, the peak file should match the cell type in the cell white list in the format of <cell_type>.bed")
    p.add_argument('--ATAC_format',type=str,dest="ATAC_format",default="bed",choices=['xiong', 'ruochi'],help="The format of the ATAC peak file")
    p.add_argument("--output_prefix", dest="output_prefix", type=str, default=None, help="The output file prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def build_MUSIC_bin():
    """
    To derive the contact frequency versus genomic distances heatmap for all frontal cortex cells we first generated 149 genomic bins spanning from 5,000 bp to 150 megabase pairs, with the nth bin spanning the genomic region 
    """
    bin_list = []
    for nn in range(151)[1:]:
        start = 2 ** (np.log2(5000) + (nn-1) * (np.log2(150000000) - np.log2(5000)) / 150)
        end = 2 **(np.log2(5000) + nn * (np.log2(150000000) - np.log2(5000)) / 150)
        bin_list.append((int(start), int(end)))
    # function to determine the bin index given the distance
    def find_bin_index(distance):
        for idx, (start, end) in enumerate(bin_list):
            if distance >= start and distance < end:
                return idx
    def find_bin_index_fast(distance):
        return int((np.log2(distance) - np.log2(5000)) / (np.log2(150000000) - np.log2(5000)) * 150)
    return bin_list, find_bin_index_fast


def build_Mattew_bin():
    """
    To count the number of cis (intra-chromosomal) contacts in each cell and bulk Hi-C data, we divided the contacts into 143 logarithmic bins, the first of which was for contacts that were separated by less than 1 kb. Each subsequent bin covered an exponent step of 0.125, using base 2.
    """
    bin_list = []
    for nn in range(144):
        start = 1000 * 2 ** (nn * 0.125)
        end = 1000 * 2 ** ((nn+1) * 0.125)
        bin_list.append((int(start), int(end)))
    # function to determine the bin index given the distance
    def find_bin_index(distance):
        for idx, (start, end) in enumerate(bin_list):
            if distance >= start and distance < end:
                return idx
    def find_bin_index_fast(distance):
        return int(np.log2(distance/1000) / 0.125)
    return bin_list, find_bin_index_fast


def parse_white_list(filename, format = 'cell_label'):
    """
    "","Idents.obj."
    "hPFC-lib-1_A1,A7","Micro"
    "hPFC-lib-1_A1,C5","OPC"
    """
    if format == 'cell_label':
        cell_id2group = pd.read_csv(filename, sep = ',', header = 0)
        cell_id2group.columns = ['cell_id', 'cell_type'] 
        # convert to dictionary with cell_id as key and cell_type as value
        cell_id2group = dict(zip(cell_id2group['cell_id'], cell_id2group['cell_type']))
    elif format == 'cell_anno':
        cell_id2anno = pd.read_csv(filename, sep = '\t', header = 0)
        # check if major column names are in the file
        if 'cell_id' not in cell_id2anno.columns or 'major' not in cell_id2anno.columns:
            print("The cell annotation file should contain cell_id and major columns")
            exit(1)
        cell_id2group = dict(zip(cell_id2anno['cell_id'], cell_id2anno['major']))
    else:
        print("Unknown format for cell white list file")
        exit(1)
    return cell_id2group


def parse_atac_peak(filename, bin_size = 1000):
    table = {}
    with open(filename, 'r') as fin:
        for line in tqdm(fin):
            row = line.strip().split()
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            if table.get(chrom) is None:
                table[chrom] = {}
            bin_idx = start // bin_size
            bin_idx_end = end // bin_size
            if bin_idx_end == bin_idx:
                bin_idx_end += 1
            for idx in range(bin_idx, bin_idx_end):
                table[chrom][idx] = True
    return table


def main():
    global args
    args = parse_arg()
    # get cell list 
    # get the white list of cells
    if args.cell_white_list is not None:
        cell_white_list = parse_white_list(args.cell_white_list, args.cell_white_list_format)
        print("Total number of selected cells: ", len(cell_white_list))
    else:
        cell_white_list = None
        print("No cell white list is provided, use all the cells in the folder")
    # list all the files in the folder
    file_list = []
    for filename in os.listdir(args.folder):
        if filename.endswith(".pair") or filename.endswith(".pair.gz"):
            cell_id = os.path.basename(filename).replace('.pair.gz', '').replace('.pair', '')
            if cell_white_list is not None:
                if cell_id in cell_white_list:
                    file_list.append(os.path.join(args.folder, filename))
            else:
                file_list.append(os.path.join(args.folder, filename))
    print("Total number of filtered files: ", len(file_list))
    file_list = sorted(file_list)
    # process all ATAC peak file in the ATAC folder
    ATAC_peak_size = 2500
    peak_anno = {}
    if args.ATAC_format == 'xiong':
        for filename in os.listdir(args.ATAC_folder):
            if filename.endswith(".bed"):
                cell_type = os.path.basename(filename).replace('.bed', '')
                peak_anno[cell_type] = parse_atac_peak(os.path.join(args.ATAC_folder, filename), bin_size = ATAC_peak_size)
        print("ATAC peak annotation has been loaded")
        print(f"Total number of cell types: {len(peak_anno)}")
    elif args.ATAC_format == 'ruochi':
        cell_type_name_convension = {'Ast': 'Astro', 'Exc': 'Exc', 'Inh': 'Inh', 'Mic': 'Micro', 'Oli': 'Oligo', 'Opc': 'OPC', 'Vas': 'VLMC'}
        for filename in os.listdir(args.ATAC_folder):
            if filename.endswith(".bed") and filename not in  ['tss6_atac_peaks.bed', 'atac_peaks.bed']:
                cell_type = os.path.basename(filename).replace('_cleaned_atac_peaks.bed', '').replace('tss6_', '')
                if 'CT_' in cell_type:
                    disease_group = 'CT'
                    cell_type = cell_type.replace('CT_', '')
                    cell_type = cell_type_name_convension[cell_type]
                else:
                    disease_group = 'AD'
                    cell_type = cell_type_name_convension[cell_type]
                if peak_anno.get(disease_group, None) is None:
                    peak_anno[disease_group] = {}
                peak_anno[disease_group][cell_type] = parse_atac_peak(os.path.join(args.ATAC_folder, filename), bin_size = ATAC_peak_size)
        peak_anno['AD']['Endo'] = peak_anno['AD']['VLMC']
        peak_anno['CT']['Endo'] = peak_anno['CT']['VLMC']
        if args.patient_id in ['12365619', '20780035', '50403446', '50410319', '66754397', '76733461', '78353027', '78452313', '84417209', '85171938']: # AD patients
            peak_anno = peak_anno['AD']
            print(f"The patient {args.patient_id} is in the AD group, use the AD ATAC peak annotation")
        else:
            peak_anno = peak_anno['CT']
            print(f"The patient {args.patient_id} is in the CT group, use the CT ATAC peak annotation")
    else:
        print(f"Unknown ATAC peak format: {args.ATAC_format}")
        exit(1)
    # determine the distance bin
    bin_list_MUSIC, find_bin_index_MUSIC = build_MUSIC_bin()
    bin_list_Mattew, find_bin_index_Mattew = build_Mattew_bin()
    #print(bin_list_MUSIC)
    #print(bin_list_Mattew)
    # read the files and count the number of reads
    out_folder = os.path.dirname(args.output_prefix)
    print(out_folder)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    fout_cell_distance = open(args.output_prefix + f".pair_distance.ATAC.{args.ATAC_format}.txt", 'w')
    print("Start to process the files")
    print("cell_id\tsample_id\tdistance_group\tbin_idx\tbin_start\tbin_end\tATAC_group\tcount\tcumu_count", file = fout_cell_distance)
    for filename in tqdm(file_list):
        cell_id = filename.split("/")[-1].replace('.pair.gz', '').replace(".pair", "")  
        cell_type = cell_white_list[cell_id]
        result = {}
        for peak_group in ['peak-peak', 'peak-none', 'none-none']:
            result[peak_group] = {'MUSIC': np.zeros(151), 'Mattew': np.zeros(144), 'read_count': 0, 'inter_chrom':0}
        if '.gz' in filename:
            fin = gzip.open(filename, 'rt')
        else:
            fin = open(filename, 'r')
        for line in fin:
            row = line.strip().split()
            chrom_1 = row[0]
            pos_1 = int(row[1])
            chrom_2 = row[2]
            pos_2 = int(row[3])
            # skip chrM
            if chrom_1 == 'chrM' or chrom_2 == 'chrM':
                continue
            result[peak_group]['read_count'] += 1
            bin_1 = pos_1 // ATAC_peak_size
            bin_2 = pos_2 // ATAC_peak_size
            # overlap with ATAC peak
            if peak_anno[cell_type].get(chrom_1, None) is not None and peak_anno[cell_type].get(chrom_2, None) is not None:
                overlap_peak_1 = peak_anno[cell_type][chrom_1].get(bin_1, None) 
                overlap_peak_2 = peak_anno[cell_type][chrom_2].get(bin_2, None)
            if overlap_peak_1 is not None and overlap_peak_2 is not None:
                peak_group = 'peak-peak'
            elif overlap_peak_1 is not None or overlap_peak_2 is not None:
                peak_group = 'peak-none'
            else:
                peak_group = 'none-none'
            # check if the two reads are in the same chromosome
            if chrom_1 != chrom_2:
                result[peak_group]['inter_chrom'] += 1
            else:
                distance = abs(pos_2 - pos_1) 
                if distance > 1000:
                    bin_idx_Mattew = find_bin_index_Mattew(distance)
                    result[peak_group]['Mattew'][bin_idx_Mattew] += 1
                    if distance > 5000 and distance < 150000000:
                        bin_idx_MUSIC = find_bin_index_MUSIC(distance)
                        result[peak_group]['MUSIC'][bin_idx_MUSIC] += 1
                    elif distance >= 150000000:
                        result[peak_group]['MUSIC'][150] += 1
        # output the result
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_MUSIC):
            start = bin_pair[0]
            end = bin_pair[1]
            count_peak2peak  = result['peak-peak']['MUSIC'][idx]
            count_peak2none = result['peak-none']['MUSIC'][idx]
            count_none2none = result['none-none']['MUSIC'][idx]
            count = count_peak2peak + count_peak2none + count_none2none
            cumu_count += count
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tPeak-Peak\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_peak2peak, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tPTCF-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_peak2none, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMUSIC\t%d\t%d\t%d\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_none2none, cumu_count), file = fout_cell_distance)
        count_peak2peak  = result['peak-peak']['MUSIC'][150]
        count_peak2none = result['peak-none']['MUSIC'][150]
        count_none2none = result['none-none']['MUSIC'][150] 
        count = count_peak2peak + count_peak2none + count_none2none
        cumu_count += count
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tPeak-Peak\t%d\t%d" % (cell_id, args.patient_id, count_peak2peak, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tPeak-Other\t%d\t%d" % (cell_id, args.patient_id, count_peak2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\t150\t150000000\t300000000\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)
        count_peak2peak  = result['peak-peak']['inter_chrom']
        count_peak2none = result['peak-none']['inter_chrom']
        count_none2none = result['none-none']['inter_chrom'] 
        count = count_peak2peak + count_peak2none + count_none2none
        cumu_count += count        
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tPeak-Peak\t%d\t%d" % (cell_id, args.patient_id, count_peak2peak, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tPeak-Other\t%d\t%d" % (cell_id, args.patient_id, count_peak2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMUSIC\tinter\t-1\t-1\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)
        cumu_count = 0
        for idx, bin_pair in enumerate(bin_list_Mattew):
            start = bin_pair[0]
            end = bin_pair[1]
            count_peak2peak  = result['peak-peak']['Mattew'][idx]
            count_peak2none = result['peak-none']['Mattew'][idx]
            count_none2none = result['none-none']['Mattew'][idx]
            count = count_peak2peak + count_peak2none + count_none2none
            cumu_count += count
            print("%s\t%s\tMattew\t%d\t%d\t%d\tPeak-Peak\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_peak2peak, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMattew\t%d\t%d\t%d\tPeak-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_peak2none, cumu_count), file = fout_cell_distance)
            print("%s\t%s\tMattew\t%d\t%d\t%d\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, idx, start, end, count_none2none, cumu_count), file = fout_cell_distance)
        count_peak2peak  = result['peak-peak']['inter_chrom']
        count_peak2none = result['peak-none']['inter_chrom']
        count_none2none = result['none-none']['inter_chrom'] 
        count = count_peak2peak + count_peak2none + count_none2none
        cumu_count += count        
        print("%s\t%s\tMattew\tinter\t-1\t-1\tPeak-Peak\t%d\t%d" % (cell_id, args.patient_id, count_peak2peak, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMattew\tinter\t-1\t-1\tPeak-Other\t%d\t%d" % (cell_id, args.patient_id, count_peak2none, cumu_count), file = fout_cell_distance)
        print("%s\t%s\tMattew\tinter\t-1\t-1\tOther-Other\t%d\t%d" % (cell_id, args.patient_id, count_none2none, cumu_count), file = fout_cell_distance)

    

if __name__=="__main__":
    main()

