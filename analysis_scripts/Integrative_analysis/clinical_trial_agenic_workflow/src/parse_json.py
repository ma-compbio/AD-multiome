#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import json
from openai import OpenAI
from pydantic import BaseModel
from tqdm import tqdm
import random
import time


def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_folder',type=str,dest="input_folder",help='Input folder with JSON files',default=None,required=True)
    p.add_argument('--whitelist',type=str,dest="whitelist",help='A file with a list of NCT IDs to process (one per line)',default=None,required=False)
    p.add_argument('--overwrite',action='store_true',dest="overwrite",help='Overwrite existing output files',default=False,required=False)
    p.add_argument('--append',action='store_true',dest="append",help='Append to existing output files',default=False,required=False)
    p.add_argument('--output_folder',type=str,dest="output_folder",help="output folder")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


class Trial(object):
    def __init__(self, data):
        """
        data is a dictionary parsed from a JSON file
        """
        # verify the data structure
        self.__verify_data_structure(data)
        self.nct_id = self.__get_nct_id(data) 
        self.title = self.__get_title(data) 
        self.status, self.last_update_time = self.__get_status(data) 
        self.description_brief = self.__get_description_brief(data)
        self.description_detailed = self.__get_description_detailed(data)
        self.phase = self.__get_phase(data)
        self.study_type = self.__get_study_type(data)
        self.primary_purpose = self.__get_primary_purpose(data)
        # results sections  
        self.results = self.__get_results(data)
        # openai section
        self.target_category = 'N/A'
        self.drug = []
        self.placebo = []
        self.explanation_target = []
        self.agent_type = 'N/A' # agent type for small molecule or biologic based on CADRO
        self.explanation_agent = []
    def __verify_data_structure(self, data):
        if 'protocolSection' not in data:
            print(f"Warning: 'protocolSection' not found in data for trial {data.get('nctId', 'N/A')}")
            if 'identificationModule' not in data.get('protocolSection', {}):
                print(f"Warning: 'identificationModule' not found in protocolSection for trial {data.get('nctId', 'N/A')}")
            if 'statusModule' not in data.get('protocolSection', {}):
                print(f"Warning: 'statusModule' not found in protocolSection for trial {data.get('nctId', 'N/A')}")
            if 'descriptionModule' not in data.get('protocolSection', {}):
                print(f"Warning: 'descriptionModule' not found in protocolSection for trial {data.get('nctId', 'N/A')}")
            if 'designModule' not in data.get('protocolSection', {}):
                print(f"Warning: 'designModule' not found in protocolSection for trial {data.get('nctId', 'N/A')}")
    def __get_nct_id(self, data):
        if 'protocolSection' in data:
            if 'identificationModule' in data['protocolSection']:
                if 'nctId' in data['protocolSection']['identificationModule']:
                    nctid = data['protocolSection']['identificationModule'].get('nctId', 'N/A')
                    return nctid
        return 'N/A'
    def __get_title(self, data):
        if 'protocolSection' in data:
            if 'identificationModule' in data['protocolSection']:
                if 'officialTitle' in data['protocolSection']['identificationModule']:
                    title = data['protocolSection']['identificationModule'].get('officialTitle', 'N/A')
                    return title
                elif 'briefTitle' in data['protocolSection']['identificationModule']:
                    title = data['protocolSection']['identificationModule'].get('briefTitle', 'N/A')
                    return title
        return 'N/A'
    def __get_status(self, data):
        if 'protocolSection' in data:
            if 'statusModule' in data['protocolSection']:
                if 'overallStatus' in data['protocolSection']['statusModule']:
                    status = data['protocolSection']['statusModule'].get('overallStatus', 'N/A')
                    last_update_time = data['protocolSection']['statusModule'].get('lastUpdateSubmitDate', 'N/A')
                    return status, last_update_time

        return 'N/A', 'N/A'
    def __get_description_brief(self, data):
        if 'protocolSection' in data:
            if 'descriptionModule' in data['protocolSection']:
                if 'briefSummary' in data['protocolSection']['descriptionModule']:
                    description = data['protocolSection']['descriptionModule'].get('briefSummary', 'N/A')
                    return description
        return 'N/A'
    def __get_description_detailed(self, data):
        if 'protocolSection' in data:
            if 'descriptionModule' in data['protocolSection']:
                if 'detailedDescription' in data['protocolSection']['descriptionModule']:
                    description = data['protocolSection']['descriptionModule'].get('detailedDescription', 'N/A')
                    return description  
        return 'N/A'
    def __get_phase(self, data):
        if 'protocolSection' in data:
            if 'designModule' in data['protocolSection']:
                if 'phases' in data['protocolSection']['designModule']:
                    phase = data['protocolSection']['designModule'].get('phases', 'N/A')
                    return phase
        return 'N/A'
    def __get_study_type(self, data):
        if 'protocolSection' in data:
            if 'designModule' in data['protocolSection']:
                if 'studyType' in data['protocolSection']['designModule']:
                    study_type = data['protocolSection']['designModule'].get('studyType', 'N/A')
                    return study_type
        return 'N/A'
    def __get_primary_purpose(self, data):
        if 'protocolSection' in data:
            if 'designModule' in data['protocolSection']:
                if 'designInfo' in data['protocolSection']['designModule']:
                    if 'primaryPurpose' in data['protocolSection']['designModule']['designInfo']:
                        primary_purpose = data['protocolSection']['designModule']['designInfo'].get('primaryPurpose', 'N/A')
                        return primary_purpose
        return 'N/A'
    def __get_results(self, data):
        if data['hasResults']:
            if 'resultsSection' in data:
                return data['resultsSection']
    def report_trial(self):
        """
        summarize the trial information into a tsv line
        """
        return f"{self.nct_id}\t{self.title}\t{self.status}\t{self.last_update_time}\t{self.phase}\t{self.study_type}\t{self.primary_purpose}"
    def dump_json_file(self, output_folder, update=False, append=False):
        """
        dump the trial information into a json file
        """
        # check output folder exists
        trial_output_folder = os.path.join(output_folder, self.nct_id)
        if not os.path.exists(trial_output_folder):
            os.makedirs(trial_output_folder)
        # check if there already exists json files in the output folder
        existing_output_file_list = []
        for fname in os.listdir(trial_output_folder):
            if fname.endswith('.json'):
                existing_output_file_list.append(os.path.join(fname))
                if not append:
                    print(f"Warning: JSON file {fname} already exists in {trial_output_folder}. Use --overwrite to overwrite.")
                    return
                else:
                    print(f"Info: JSON file {fname} already exists in {trial_output_folder}. Adding more records.")
        # determine the new output file name
        output_file = os.path.join(trial_output_folder, f"{self.nct_id}.{len(existing_output_file_list)+1}.json")
        # check if output file exists
        if os.path.exists(output_file):
            if update:
                print(f"Info: Output file {output_file} already exists. Overwrite it.")
                with open(output_file, 'w') as f:
                    json.dump(self.__dict__, f, indent=4)
            else:
                print(f"Warning: Output file {output_file} already exists. Use --overwrite to overwrite.")
        else:
            with open(output_file, 'w') as f:
                json.dump(self.__dict__, f, indent=4)


class TargetCategoryExtraction(BaseModel):
    target_catetory: str
    drug: list[str]
    placebo: list[str]
    explanation: list[str]


class AgentTypeExtraction(BaseModel):
    agent_type: str
    explanation: list[str]


def main():
    global args
    args = parse_arg()
    # list all json file in the input folder
    json_file_list = [f for f in os.listdir(args.input_folder) if f.endswith('.json')]
    json_file_list.sort()
    print(f"Found {len(json_file_list)} JSON files in the input folder.")
    trial_list = []
    for json_file in tqdm(json_file_list):
        json_path = os.path.join(args.input_folder, json_file)
        with open(json_path, 'r') as f:
            data = json.load(f)
            trial = Trial(data)
            trial_list.append(trial)
    print(f"Parsed {len(trial_list)} trials.")
    # use openai to summarize the trials
    client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    #openai_get_target_from_summary(client, trial_list[1372])
    #openai_get_target_from_summary_react(client, trial_list[1372])
    #exit(1)
    output_json_folder = os.path.join(args.output_folder, 'json_output')
    if not os.path.exists(output_json_folder):
        os.makedirs(output_json_folder)
    #for trial in random.sample(trial_list, 1):
    for trial in trial_list:
        print(trial.report_trial() + '\n')
        # check if the output file already exists
        trial_output_folder = os.path.join(output_json_folder, trial.nct_id)
        if os.path.exists(trial_output_folder):
            existing_output_file_list = [f for f in os.listdir(trial_output_folder) if f.endswith('.json')]
            if len(existing_output_file_list) > 0 and not args.overwrite and not args.append:
                print(f"Warning: Output files already exist for trial {trial.nct_id}. Use --overwrite or --append to process it.")
                continue
        # use prompt to summarize the trail description to get the target category
        try:
            openai_get_target_from_summary_react(client, trial)
        except Exception as e:
            print(f"Error processing trial {trial.nct_id}: {e}")
            continue
        # based on the target category, combine the trail information to further get the agent type
        if trial.target_category != 'N/A':
            print(f"Trial {trial.nct_id} has target category {trial.target_category}, proceeding to determine agent type.")
            try:
                openai_get_agent_type_react(client, trial)
            except Exception as e:
                print(f"Error processing trial {trial.nct_id} for agent type: {e}")
                continue
        # save the trial information into the output file
        trial.dump_json_file(output_json_folder, update = args.overwrite, append = args.append)
        # sleep for a while to avoid rate limit
        sleep_time = random.randint(5, 10)
        print(f"Sleeping for {sleep_time} seconds to avoid rate limit.")
        time.sleep(sleep_time)
    print("All Done.")


def openai_get_target_from_summary(client, trial):
    """
    The purpose of this function is to use openai to read the trial description and summarize the target of trial into four categories.
    """
    # print the trial information
    print(trial.report_trial())
    # set the schema for the output in json format
    schema = {
        "type": "object",
        "properties": {
            "target_category": {
                "type": "string",
                "description": "The target category of the clinical trial. If no target category can be inferred, respond with 'N/A'."
            }
        },
        "required": ["target_category"]
    }
    # use openai to summarize the trial
    print(">>>>>")
    print(f"trial description: {trial.description_brief}")
    print("")
    print("")
    print(f"trail description detailed: {trial.description_detailed}")
    print("<<<<<")
    #response = client.chat.completions.create(
    response = client.responses.parse(
        #model = "gpt-5-mini",
        model = "gpt-4o",
        tools = [{"type": "web_search"}],
        input = [
            {"role": "system", 
             "content": (
                 "You are an expert in clinical trial design and analysis. "
                 "Summarize the target category of the clinical trial for Alzheimer's disease from the provided description. You should fetch the drug target from the description if possible. "
                 "The four categories are: "
                 "1) disease-targeted biologic; "
                 "2) disease-targeted small molecule; "
                 "3) cognitive enhancer; "
                 "4) neuropsychiatric symptom improvement. "
                 "Return the result as a JSON object with a 'category' key. "
                 "If no target category can be inferred or the trial cannot be classified into any group (eg, PET scan for prediction), return 'N/A'."
             )},
              {"role": "user", "content": f"Extract the target category according to the following clinical trial description: {trial.description_brief} and title: {trial.title}. If the drug is not known, use web search to find the drug information. Also output the web search results as explanation."}
        ],
        #temperature=0.1,
        text_format = TargetCategoryExtraction,
        #response_format={"type": "json_schema", "schema": schema}
    )
    # print the response
    print(response.output_parsed)
    print("")
    # parse the response into the TargetCategoryExtraction model
    result_parse = response.output_parsed
    # get the "parsed" field from the response
    trial.target_category = result_parse.target_catetory
    trial.drug = result_parse.drug
    trial.placebo = result_parse.placebo
    trial.explanation_target = result_parse.explanation
  

def openai_get_target_from_summary_react(client, trial):
    """
    The purpose of this function is to use openai to read the trial description and summarize the target of trial into four categories.
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
                 "You are an expert in clinical trial design and analysis. Your task is to summarize the target category of a clinical trial for Alzheimer's disease based on the provided description, identifying the drug target if possible. The four possible categories are: "
                 "1) disease-targeted biologic; "
                 "2) disease-targeted small molecule; "
                 "3) cognitive enhancer; "
                 "4) neuropsychiatric symptom improvement. "
                 "Return the result in the format specified by `text_format`. If no target category can be inferred or the trial does not fit any category (e.g., PET scan for prediction), return 'N/A'. "
                 "Solve this using the ReAct (Reasoning and Acting) approach, structured as follows: "
                 "**Reason**: Analyze the description to identify the trial’s intervention, mechanism of action, and intended effect (e.g., disease modification, cognitive improvement, or symptom management). Define the categories: "
                 "- **Disease-targeted biologic**: Biologics (e.g., monoclonal antibodies, vaccines) targeting Alzheimer’s pathology (e.g., amyloid, tau). "
                 "- **Disease-targeted small molecule**: Small molecule drugs (e.g., inhibitors, modulators) targeting Alzheimer’s pathology. "
                 "- **Cognitive enhancer**: Drugs improving cognitive function (e.g., memory, attention) without targeting pathology. "
                 "- **Neuropsychiatric symptom improvement**: Interventions alleviating behavioral or psychiatric symptoms (e.g., agitation, depression). "
                 "If the trial involves diagnostics (e.g., PET scans) or other non-therapeutic interventions, classify as 'N/A'. "
                 "**Act**: Extract key details (e.g., drug name, type, or intervention) from the description and match to a category or 'N/A'. Note any ambiguity or missing information. If the drug is not known, use web search to find the drug information. "
                 "**Reflect**: Verify the classification aligns with the description and category definitions. Reconsider for overlooked details or alternative interpretations if uncertain. Confirm 'N/A' if no category fits. "
                 "**Output**: Return the result in the `text_format` with the category or 'N/A', accompanied by a brief explanation of the reasoning process for transparency. "
             )},
              {"role": "user", "content": f"Extract the target category according to the following clinical trial description: {trial.description_brief} and title: {trial.title}. If the drug is not known, use web search to find the drug information. Also output the web search results as explanation."}
        ],
        #temperature=0.1,
        text_format = TargetCategoryExtraction,
        #response_format={"type": "json_schema", "schema": schema}
    )
    # print the response
    #print("")
    #print(">>>>>")
    #print(response.output_parsed)
    #print("<<<<<")
    #print("")
    # parse the response into the TargetCategoryExtraction model
    result_parse = response.output_parsed
    # get the "parsed" field from the response
    trial.target_category = result_parse.target_catetory
    trial.drug = result_parse.drug
    trial.placebo = result_parse.placebo
    trial.explanation_target = result_parse.explanation
  

def openai_get_agent_type_react(client, trial):
    """
    The purpose of this function is to use openai to read the trial description and summarize the agent type based on CADRO classification.
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
                 "You are an expert in Alzheimer's disease. Your task is to determine the drug target type for a clinical trial based on the provided description, using the Common Alzheimer's and Related Dementias Research Ontology (CADRO) classification. The CADRO categories are: "
                 "A) Amyloid beta; B) Tau; C) ApoE, Lipids and Lipoprotein Receptors; D) Neurotransmitter Receptors; E) Neurogenesis; F) Inflammation; G) Oxidative Stress; H) Cell Death; I) Proteostasis/Proteinopathies; J) Metabolism and Bioenergetics; K) Vasculature; L) Growth Factors and Hormones; M) Synaptic Plasticity/Neuroprotection; N) Gut-Brain Axis; O) Circadian Rhythm; P) Environmental Factors; Q) Epigenetic Regulators; R) Multi-target; S) Unknown Target; T) Other. "
                 "Return the result in the format specified by `text_format`. If no target can be inferred or the trial does not fit any CADRO category (e.g., PET scan for prediction), return 'T) Other'. "
                 "Solve this using the ReAct (Reasoning and Acting) approach, structured as follows: "
                 "**Reason**: Analyze the description to identify the trial’s intervention, mechanism of action, or biological focus (e.g., targeting amyloid plaques, inflammation, or synaptic function). Compare these details against the CADRO categories to determine the most specific match. Consider whether the trial involves multiple targets (R) or a non-therapeutic focus (T). "
                 "**Act**: Extract relevant details (e.g., drug name, biological target, or intervention type) from the description and assign the appropriate CADRO category. If the description is ambiguous, note missing information. If multiple targets are indicated, consider 'R) Multi-target'. If the intervention is non-therapeutic (e.g., diagnostics), assign 'T) Other'. "
                 "**Reflect**: Verify the classification aligns with the description and CADRO definitions. Reassess for overlooked details or alternative interpretations. Confirm 'T) Other' if no therapeutic target fits or if information is insufficient. "
                 "**Output**: Return the result in the `text_format` with the CADRO category letter and name (e.g., 'A) Amyloid beta') or 'T) Other', accompanied by a brief explanation of the reasoning process for transparency. "
                 "Example: "
                 "- Description: 'A trial testing a monoclonal antibody targeting amyloid-beta plaques.' "
                 "- Reason: The trial uses a monoclonal antibody targeting amyloid-beta, a key Alzheimer’s pathology. "
                 "- Act: Assign 'A) Amyloid beta' based on the explicit amyloid-beta target. "
                 "- Reflect: Confirm the category fits, as the focus is solely on amyloid-beta, not multiple targets or non-therapeutic interventions. "
                 "- Output (assuming `text_format` is JSON): ```json\n{\"category\": \"A) Amyloid beta\"}\n``` "
                 "Explanation: The monoclonal antibody targets amyloid-beta plaques, aligning with CADRO category A. "
             )},
              {"role": "user", "content": f"Classifity drug target type according to the following clinical trial description: title: {trial.title}, trial category: {trial.target_category}, trial drug: {trial.drug}, trail target explanation: {trial.explanation_target}. If the drug or its target gene or pathway is not known, use web search to find the drug information. Also output the web search results as explanation."}
        ],
        #temperature=0.1,
        text_format = AgentTypeExtraction,
        #response_format={"type": "json_schema", "schema": schema}
    )
    # print the response
    print(response.output_parsed)
    print("")
    # parse the response into the TargetCategoryExtraction model
    result_parse = response.output_parsed
    # get the "parsed" field from the response
    trial.agent_type = result_parse.agent_type
    trial.explanation_agent = result_parse.explanation


def openai_web_search_for_drug(client, query):
    """
    The purpose of this function is to use openai to search the web for the query.
    """
    response = client.responses.create(
        model = "gpt-5",
        tools = [{"type": "web_search"}],
        input = f"Provide a short summary (one or two sentence) for the following drug : {query}."
    )
    print(response.output_text)


if __name__=="__main__":
    main()

