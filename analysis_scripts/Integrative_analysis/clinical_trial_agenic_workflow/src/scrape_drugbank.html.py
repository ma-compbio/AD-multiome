#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import json
from bs4 import BeautifulSoup
import pandas as pd
from tqdm import tqdm


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_html',type=str,dest="input_html",help="input html file",default=None,required=True)
    p.add_argument('--output',type=str,dest="output",help="output file",default=None,required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def main():
    global args
    args = parse_arg()
    # read the HTML file and parse it with BeautifulSoup
    html_content_list = []
    with open(args.input_html, 'r', encoding='utf-8') as fin:
        block = ""
        for line in fin:
            if line.strip() == '':
                if block != "":
                    html_content_list.append(block)
                    block = ""
                continue
            elif line.strip().startswith('<table'):
                block = ""
            else:
                pass
            block += line 
        if block != "":
            html_content_list.append(block)
    print(f"Total {len(html_content_list)} HTML blocks found.")
    # loop over each HTML block and extract tables
    final_result = []
    for html_string in html_content_list:
        soup = BeautifulSoup(html_string, 'html.parser')
        # go through the html structure to find all tables
        table = soup.find('table', {'id': 'DataTables_Table_1'})
        # extract the header
        header_table = {}
        headers = [header.text.strip() for header in table.find_all('th')]
        for i, header in enumerate(headers):
            header_table[i] = header
        result_list = []
        # extract the rows
        for row in table.find_all('tr')[1:]:  # Skip the header row
            cells = row.find_all('td')
            if len(cells) == len(headers):  # Ensure the row has the correct number of cells
                # if one cell has ul list, extract the text from each li
                row = {}
                for row_idx, cell in enumerate(cells):
                    col_name = header_table[row_idx]
                    if col_name == 'Identifier':
                        row[col_name] = cell.text.strip()
                    elif col_name == 'Title':
                        row[col_name] = cell.text.strip()
                    elif col_name == 'Drug(s)':
                        # if ul in cell
                        row[col_name] = {}
                        if cell.find('ul'):
                            # find drug name and drug ID
                            for li in cell.find_all('li'):
                                drug_name = li.text.strip()
                                # check if there is a link
                                if li.find('a') is not None :
                                    drug_id = li.find('a')['href'].split('/')[-1]
                                else:
                                    drug_id = 'NA'
                                    print(f"Warning: No drug ID found for drug {drug_name} in trial {row['Identifier']}")
                                row[col_name][drug_name] = drug_id
                        else:
                            drug_name = cell.text.strip()
                            if drug_name != 'No drug interventions':
                                drug_id = cell.find('a')['href'].split('/')[-1]
                                row[col_name][drug_name] = drug_id
                            else:
                                row[col_name] = {}
                    elif col_name == 'Purpose':
                        row[col_name] = cell.text.strip()
                    elif col_name == 'Phase': 
                        row[col_name] = cell.text.strip()
                    else:
                        row[col_name] = cell.text.strip()
                    #
                result_list.append(row)
        # save the result
        final_result.extend(result_list)
    # save the final result
    trial2drug = {}
    for trial in final_result:
        nct_id = trial['Identifier']
        trial2drug[nct_id] = {}
        for col in ['Drug(s)', 'Title', 'Purpose', 'Phase']:
            trial2drug[nct_id][col] = trial[col]
    # convert dict to json file
    json.dump(trial2drug, open(args.output, 'w'), indent=4)
     
    
if __name__=="__main__":
    main()

