
# LLM-based agent framework linking AD clinical trials, mechanisms, and target genes.

This README describes the workflow to use OpenAI's GPT models to analyze clinical trials related to Alzheimer's Disease (AD) and integrate this information with drug data from DrugBank. The workflow includes web scraping, information extraction, summarization, integration, website building, and analysis.
Using the generated structured data, researchers can gain insights into clinical trials and their associated drugs for AD.
I also vibe code a website to visualize the clinical trial information.
The links to the repository and website are provided below:
Repository: [https://github.com/zocean/ad-clinical-trials-explorer](https://github.com/zocean/ad-clinical-trials-explorer)
Demo website: [https://zocean.github.io/ad-clinical-trials-explorer/](https://zocean.github.io/ad-clinical-trials-explorer/)

``` bash
##################
# Use OpenAI to anotate and analyze clinical trials, and integrate with Drugbank data for Alzheimer's Disease
##################
#https://clinicaltrials.gov/search?cond=Alzheimer%27s%20Disease&term=Alzheimer%20Disease&aggFilters=status:com%20act%20rec,studyType:int

# put your OPENAI API key here
OPENAI_API_KEY=""
export OPENAI_API_KEY
# step 0: web scrapting Drugbank to get the clinical and drug table
python src/scrape_drugbank.html.py --input_html script_custom/clinical_trial/data/drug_bank_AD/web_scrape.html --output script_custom/clinical_trial/results/drugbank_clinial2drug.json
# step 1: use Agent to extract information from clinical trial json files
python src/parse_json.py --input_folder script_custom/clinical_trial/data/AD_json --output script_custom/clinical_trial/results/AD_clinical_trial_annotation.tsv --output_folder script_custom/clinical_trial/results
# step 2: summarize the clinical trial information
python src/summarize_clinical_trial.py --input_json_folder script_custom/clinical_trial/results/json_output --output_summary script_custom/clinical_trial/results/AD_clinical_trial_annotation.tsv
# step 3: parse the full Drugbank xml file to get the drug information table
python src/parse_drugbank.py --input_xml script_custom/clinical_trial/data/drugbank_all_full_database.5.1.13.xml --output script_custom/clinical_trial/data/drugbank_all_full_database.5.1.13.csv
# step 4: integrate the clinical trial and drug information
python src/integrate_clinical_drug.py --input_clinical_trial script_custom/clinical_trial/results/AD_clinical_trial_annotation.tsv --input_drugbank script_custom/clinical_trial/data/drugbank_all_full_database.5.1.13.csv --input_trial2drug script_custom/clinical_trial/results/drugbank_clinial2drug.json --local_database_prefix script_custom/clinical_trial/results/AD_clinical_drug_local_db --output script_custom/clinical_trial/results/AD_clinical_drug_integration --input_openai_folder script_custom/clinical_trial/results/json_output --output_openai_folder script_custom/clinical_trial/results/json_output_drug
# step 5: build a website to visualize the clinical trial information
python src/build_website.v2.py --input_json_folder script_custom/clinical_trial/results/json_output --input_drug_target_file script_custom/clinical_trial/results/AD_clinical_drug_integration.tsv --input_drug_target_openai_folder script_custom/clinical_trial/results/json_output_drug --output_folder script_custom/clinical_trial/results/website
# step 6: analyze the clinical trial integrated results
python src/analyze_clinical_trial.prep.py --input_prefix script_custom/clinical_trial/results/AD_clinical_drug_integration --output_prefix script_custom/clinical_trial/results/AD_clinical_drug_gene_analysis
python src/analyze_clinical_trial.enrichment.py --input_drugbank_db script_custom/clinical_trial/results/AD_clinical_drug_local_db.json --input_trial_summary script_custom/clinical_trial/results/AD_clinical_drug_gene_analysis.summary.tsv --config_deg script_custom/clinical_trial/config/config_deg_MAST.json --config_trial script_custom/clinical_trial/config/config_trial_gene_enrichment.json --output_prefix script_custom/clinical_trial/results/AD_clinical_drug_enrichment.analysis
python src/analyze_clinical_trial.enrichment.geneset.py --input_drugbank_db script_custom/clinical_trial/results/AD_clinical_drug_local_db.json --input_trial_summary script_custom/clinical_trial/results/AD_clinical_drug_gene_analysis.summary.tsv --config_geneset script_custom/clinical_trial/config/config_genset.json --config_trial script_custom/clinical_trial/config/config_trial_gene_enrichment.json --output_prefix script_custom/clinical_trial/results/AD_clinical_drug_enrichment.analysis_geneset
```
