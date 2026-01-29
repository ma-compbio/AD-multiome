#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 03 Oct 2025 11:15:27 PM

import os,sys,argparse
import json
import html
import re
import pandas as pd
from tqdm import tqdm


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_json_folder', type=str, help='Input JSON folder', required=True)
    p.add_argument('--input_drug_target_file', type=str, help='Input drug target JSON file', required=False)
    p.add_argument('--input_drug_target_openai_folder', type=str, help='Input drug target OpenAI JSON folder', required=False)
    p.add_argument('--show_all', action='store_true', help='Show all JSON files, not just the latest one in each folder')
    p.add_argument('--output_folder', type=str, help='Output folder', required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


def load_trial_target_json(input_json_folder_target):
    """
    load trial target json files
    """
    result = {}
    for file in os.listdir(input_json_folder_target):
        if file.endswith('.json'):
            json_path = os.path.join(input_json_folder_target, file)
            data = json.load(open(json_path))
            trial_id = data.get('nct_id', 'N/A')
            if trial_id != 'N/A':
                result[trial_id] = data
    return result


def main():
    global args
    args = parse_arg()
    # create output folder
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    # create an index.html file to link all json files and show a summary of clinical trials
    json_files = []
    for folder in sorted(os.listdir(args.input_json_folder)):
        # check if folder is a directory
        folder_path = os.path.join(args.input_json_folder, folder)
        if os.path.isdir(folder_path):
            tmp_files = []
            for file in os.listdir(folder_path):
                if file.endswith('.json'):
                    json_files.append(os.path.join(folder, file))
            # sort tmp files by the run id
            tmp_files = sorted(tmp_files, key=lambda x: int(x.split('.')[-2]))
            # get the last file
            if tmp_files and args.show_all is False:
                json_files.append(tmp_files[-1])
            else:
                json_files.extend(tmp_files)
    print(f"Found {len(json_files)} JSON files in {args.input_json_folder}")
    index_file = os.path.join(args.output_folder, 'index.html')
    # load drug target json file
    drug_gene_target_df = pd.read_csv(args.input_drug_target_file, sep = '\t', header = 0) if args.input_drug_target_file else None
    drug_gene_target_table = {}
    for i, row in drug_gene_target_df.iterrows():
        trial_id = row['nct_id']
        drug_target_group =  row['drug_target_group']
        gene_list = row['gene_symbols'] if row['gene_symbols'] != 'NA' else []
        confidence_level = row['confidence_level']
        drug_gene_target_table[trial_id] = {'drug_target_group': drug_target_group, 'gene_list': gene_list, 'confidence_level': confidence_level}
    # load openai json files
    trial2target_openai = load_trial_target_json(args.input_drug_target_openai_folder) if args.input_drug_target_openai_folder else {}
    create_index_html_pretty(json_files, index_file, drug_gene_target_table, trial2target_openai)


def load_trial_target_json(input_json_folder_target):
    """
    load explanation json files
    """
    result = {}
    for file in os.listdir(input_json_folder_target):
        if file.endswith('.json'):
            json_path = os.path.join(input_json_folder_target, file)
            data = json.load(open(json_path))
            trial_id = data.get('nct_id', 'N/A')
            # get explanation 
            explanation_drug_gene_target = data.get('explanation_target_gene', [])
            result[trial_id] = explanation_drug_gene_target
    return result


def create_index_html(json_files, index_file, drug_gene_target_table=None, trial2target_openai=None):
    """
    Create an index.html file to link all json files and show a summary of clinical trials
    """
    index_content = """
    <html>
    <head>
        <title>Clinical Trials Summary</title>
        <style>
            table {
                border-collapse: collapse;
                width: 100%;
            }
            th, td {
                border: 1px solid black;
                padding: 8px;
                text-align: left;
            }
            th {
                background-color: #f2f2f2;
            }
        </style>
    </head>
    <body>
        <h1>Clinical Trials Summary</h1>
        <table>
            <tr>
                <th>Trial ID</th>
                <th>Run ID</th>
                <th>Trial title</th>
                <th>Staus</th>
                <th>Phase</th>
                <th>Primary purpose</th>
                <th>Target category</th>
                <th>Drugs</th>
                <th>Agent type</th> 
                <th>Official link</th>
            </tr>
    """
    for json_file in json_files:
        json_path = os.path.join(args.input_json_folder, json_file)
        # get run_id
        run_id = json_file.split('.')[-2]
        # get information from json file
        data = json.load(open(json_path))
        trial_id = data.get('nct_id', 'N/A')
        if trial_id == 'NCT05269394':
            print(data)
        trial_title = data.get('title', 'N/A')
        trial_status = format_text_first_letter_cap(data.get('status', 'N/A'))
        trial_phase = format_text_first_letter_cap(data.get('phase', 'N/A'))
        trial_description = html.escape(clean_citation(data.get('description_brief', 'N/A')))
        trial_purpose = format_text_first_letter_cap(data.get('primary_purpose', 'N/A'))
        trial_drugs = ', '.join(data.get('drug', [])) if data.get('drug') else 'N/A'
        trial_target = data.get('target_category', 'N/A')
        trial_agent_type = data.get('agent_type', 'N/A') 
        trial_explanation_target = html.escape(clean_citation(' '.join(data.get('explanation_target', []))))
        trial_explanation_agent = html.escape(clean_citation(' '.join(data.get('explanation_agent', []))))
        # add drug gene target information
        if drug_gene_target_table is not None:
            trial_drug_target_gene_group = drug_gene_target_table['target_group']
            trial_drug_target_gene_list = str(drug_gene_target_table['gene_list'])
            trial_drug_target_gene_confidence = drug_gene_target_table['confidence_level']
            if trial2target_openai is not None and trial_id in trial2target_openai:
                explanation_drug_gene_target = html.escape(clean_citation(' '.join(trial2target_openai[trial_id])))
            else:
                explanation_drug_gene_target = None
        else:
            trial_drug_target_gene_group = 'N/A'
            trial_drug_target_gene_list = 'N/A'
            trial_drug_target_gene_confidence = 'N/A'
            explanation_drug_gene_target = None
        # copy json file to output folder
        os.makedirs(os.path.join(args.output_folder, os.path.dirname(json_file)), exist_ok=True)
        # only copy files that do not exist
        if not os.path.exists(os.path.join(args.output_folder, json_file)):
            os.system(f"cp {json_path} {os.path.join(args.output_folder, json_file)}")
        # make the json path relative to the index.html file
        html_json_path = os.path.relpath(os.path.join(args.output_folder, json_file), os.path.dirname(index_file))
        # create link to official clinical trial page https://clinicaltrials.gov/study/NCT03378245
        official_link = f"https://clinicaltrials.gov/study/{trial_id}" if trial_id != 'N/A' else '#'
        index_content += f"""
            <tr>
                <td style="max-width:80px;"><a href="{html_json_path}">{trial_id}</a></td>
                <td style="max-width:30px;">{run_id}</td>  
                <td style="max-width:200px; white-space:pre-wrap; word-wrap:break-word;" title="{trial_description}">{trial_title}</td>
                <td style="max-width:60px; white-space:pre-wrap; word-wrap:break-word;">{trial_status}</td>
                <td style="max-width:100px;">{trial_phase}</td>
                <td style="max-width:50px;">{trial_purpose}</td>
                <td style="max-width:160px;" title="{trial_explanation_target}">{trial_target}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;">{trial_drugs}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;" title="{trial_explanation_agent}">{trial_agent_type}</td>
                <td style="max-width:60px; white-space:pre-wrap; word-wrap:break-word;">{trial_drug_target_gene_group}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;" title="{explanation_drug_gene_target}">{trial_drug_target_gene_list} |({trial_drug_target_gene_confidence})</td>
                <td><a href="{official_link}" target="_blank">link</a></td>
            </tr>
        """
    index_content += """
        </table>
    </body>
    </html>
    """
    with open(index_file, 'w') as f:
        f.write(index_content)
    print(f"Index HTML file created at {index_file}")


def create_index_html_pretty(json_files, index_file, drug_gene_target_table=None, trial2target_openai=None):
    index_content = """
    <html>
    <head>
        <title>Clinical Trials Summary</title>
        <style>
            table {
                border-collapse: collapse;
                width: 100%;
            }
            th, td {
                border: 1px solid black;
                padding: 8px;
                text-align: left;
            }
            th {
                background-color: #f2f2f2;
            }
            #pagination {
                margin: 20px 0;
                text-align: center;
            }
            #pagination button {
                margin: 0 5px;
                padding: 8px 12px;
                border: 1px solid #ccc;
                background: #f9f9f9;
                cursor: pointer;
            }
            #pagination button:hover {
                background: #e9e9e9;
            }
            #pagination button:disabled {
                background: #ddd;
                cursor: not-allowed;
            }
            input#search {
                padding: 8px;
                width: 300px;
                margin-bottom: 10px;
            }
        </style>
    </head>
    <body>
        <h1>Clinical Trials Summary</h1>
        <div>
            <input type="text" id="search" placeholder="Search trials (filters across all columns)..." onkeyup="applyFilter()">
        </div>
        <table>
            <tr>
                <th>Trial ID</th>
                <th>Run ID</th>
                <th>Trial title</th>
                <th>Status</th>
                <th>Phase</th>
                <th>Primary purpose</th>
                <th>Target category</th>
                <th>Drugs</th>
                <th>Agent type</th> 
                <th>Drug target group</th>
                <th>Drug target genes |(confidence)</th>
                <th>Official link</th>
            </tr>
    """
    for json_file in tqdm(json_files):
        json_path = os.path.join(args.input_json_folder, json_file)
        # get run_id
        run_id = json_file.split('.')[-2]
        # get information from json file
        data = json.load(open(json_path))
        trial_id = data.get('nct_id', 'N/A')
        trial_title = data.get('title', 'N/A')
        trial_status = format_text_first_letter_cap(data.get('status', 'N/A'))
        trial_phase = format_text_first_letter_cap(data.get('phase', 'N/A'))
        trial_description = html.escape(clean_citation(data.get('description_brief', 'N/A')))
        trial_purpose = format_text_first_letter_cap(data.get('primary_purpose', 'N/A'))
        trial_drugs = ', '.join(data.get('drug', [])) if data.get('drug') else 'N/A'
        trial_target = data.get('target_category', 'N/A')
        trial_agent_type = data.get('agent_type', 'N/A') 
        trial_explanation_target = html.escape(clean_citation(' '.join(data.get('explanation_target', []))))
        trial_explanation_agent = html.escape(clean_citation(' '.join(data.get('explanation_agent', []))))
        # add drug gene target information
        if drug_gene_target_table is not None:
            trial_drug_target_gene_group = drug_gene_target_table[trial_id]['drug_target_group']
            trial_drug_target_gene_list = str(drug_gene_target_table[trial_id]['gene_list'])
            trial_drug_target_gene_confidence = drug_gene_target_table[trial_id]['confidence_level']
            if trial2target_openai is not None and trial_id in trial2target_openai:
                explanation_drug_gene_target = html.escape(clean_citation(' '.join(trial2target_openai[trial_id])))
            else:
                explanation_drug_gene_target = None
        else:
            trial_drug_target_gene_group = 'N/A'
            trial_drug_target_gene_list = 'N/A'
            trial_drug_target_gene_confidence = 'N/A'
            explanation_drug_gene_target = None
        # make json output folder
        os.makedirs(os.path.join(args.output_folder, 'data'), exist_ok=True)
        # copy json file to output folder
        os.makedirs(os.path.join(args.output_folder, 'data', os.path.dirname(json_file)), exist_ok=True)
        # only copy files that do not exist
        if not os.path.exists(os.path.join(args.output_folder, json_file)):
            os.system(f"cp {json_path} {os.path.join(args.output_folder, 'data', json_file)}")
        # make the json path relative to the index.html file
        html_json_path = os.path.relpath(os.path.join(args.output_folder, 'data', json_file), os.path.dirname(index_file))
        # create link to official clinical trial page https://clinicaltrials.gov/study/NCT03378245
        official_link = f"https://clinicaltrials.gov/study/{trial_id}" if trial_id != 'N/A' else '#'
        index_content += f"""
            <tr>
                <td style="max-width:80px;"><a href="{html_json_path}">{trial_id}</a></td>
                <td style="max-width:30px;">{run_id}</td>  
                <td style="max-width:200px; white-space:pre-wrap; word-wrap:break-word;" title="{trial_description}">{trial_title}</td>
                <td style="max-width:60px; white-space:pre-wrap; word-wrap:break-word;">{trial_status}</td>
                <td style="max-width:80px;">{trial_phase}</td>
                <td style="max-width:50px;">{trial_purpose}</td>
                <td style="max-width:160px;" title="{trial_explanation_target}">{trial_target}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;">{trial_drugs}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;" title="{trial_explanation_agent}">{trial_agent_type}</td>
                <td style="max-width:60px; white-space:pre-wrap; word-wrap:break-word;">{trial_drug_target_gene_group}</td>
                <td style="max-width:120px; white-space:pre-wrap; word-wrap:break-word;" title="{explanation_drug_gene_target}">{trial_drug_target_gene_list} |({trial_drug_target_gene_confidence})</td>
                <td><a href="{official_link}" target="_blank">link</a></td>
            </tr>
        """
    index_content += """
            </table>
            <div id="pagination"></div>
            <script>
                let currentPage = 1;
                const rowsPerPage = 10; // Adjust as needed

                function applyFilter() {
                    const search = document.getElementById('search').value.toLowerCase();
                    const table = document.querySelector('table');
                    const dataRows = Array.from(table.querySelectorAll('tr:nth-child(n+2)'));
                    let visibleRows = dataRows.filter(row => {
                        return Array.from(row.cells).some(cell => 
                            cell.textContent.toLowerCase().includes(search)
                        );
                    });
                    const totalPages = Math.ceil(visibleRows.length / rowsPerPage);
                    displayPage(visibleRows, currentPage, rowsPerPage);
                    updatePagination(totalPages);
                }

                function displayPage(visibleRows, page, perPage) {
                    const start = (page - 1) * perPage;
                    const end = start + perPage;
                    const table = document.querySelector('table');
                    const allDataRows = table.querySelectorAll('tr:nth-child(n+2)');
                    allDataRows.forEach(row => row.style.display = 'none');
                    for (let i = start; i < Math.min(end, visibleRows.length); i++) {
                        visibleRows[i].style.display = '';
                    }
                }

                function updatePagination(totalPages) {
                    let pag = document.getElementById('pagination');
                    pag.innerHTML = '';
                    if (totalPages <= 1) return;
                    // Previous button
                    if (currentPage > 1) {
                        let btn = document.createElement('button');
                        btn.textContent = 'Previous';
                        btn.onclick = () => { currentPage--; applyFilter(); };
                        pag.appendChild(btn);
                    }
                    // Page numbers
                    for (let i = 1; i <= totalPages; i++) {
                        let btn = document.createElement('button');
                        btn.textContent = i;
                        if (i === currentPage) btn.disabled = true;
                        btn.onclick = () => { currentPage = i; applyFilter(); };
                        pag.appendChild(btn);
                    }
                    // Next button
                    if (currentPage < totalPages) {
                        let btn = document.createElement('button');
                        btn.textContent = 'Next';
                        btn.onclick = () => { currentPage++; applyFilter(); };
                        pag.appendChild(btn);
                    }
                }

                // Initialize on page load
                window.onload = applyFilter;
            </script>
        </body>
        </html>
    """
    with open(index_file, 'w') as f:
        f.write(index_content)
    print(f"Index HTML file created at {index_file}")


def format_text_first_letter_cap(text):
    ''' 
    format text to have the first letter of each word capitalized
    '''
    # if text is a list, join it with comma
    cleaned_text = ""
    if text == 'N/A' or text is None or text == ['N/A'] or text == 'NA' or text == ['NA']:
        return 'N/A'
    if isinstance(text, list):
        for sentence in text:
            cleaned_text += ' '.join([word.capitalize() for word in sentence.split(' ')]) + ', '
    elif isinstance(text, str):
        cleaned_text = ''.join([word.capitalize() for word in text.split(' ')])
    else:
        cleaned_text = str(text)
    # remove the last comma and space
    cleaned_text = cleaned_text.rstrip(', ')
    return cleaned_text


def clean_citation(text):
    ''' 
    remove text like \ue200cite\ue202turn0search7\ue201
    '''
    cleaned_text = re.sub(r'\ue200cite\S*?\ue201', '', text)
    return cleaned_text


if __name__=="__main__":
    main()

