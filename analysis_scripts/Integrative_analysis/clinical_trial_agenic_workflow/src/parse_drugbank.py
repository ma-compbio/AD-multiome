#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 02 Feb 2023 11:27:48 PM

import os,sys,argparse
import xml.etree.ElementTree as ET
import time
from tqdm import tqdm
import pandas as pd

"""
This script use example script from https://gist.github.com/rosherbal/56461421c69a8a7da775336c95fa62e0
"""

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input_xml', type=str, help='input xml file', required=True)
    p.add_argument('--output', type=str, help='output file', required=True)
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()


# Classes #
class Drug:
    """
    docstring for Drug.
    """
    def __init__(self, features):

        self.id = features['id']
        self.name = features['name']
        self.synonyms = features['synm']
        self.kingdom = features['kgd']
        self.superclass = features['sclass']
        self.interaction = features['itrc']
        self.external_id = features['ext_id']
        self.pathways = features['pathways']
        self.target = []
        self.target_source = None # DrugBank-direct, DrugBank-match, or OpenAI

    def getDrugfeatures(self):
        drug_dict = {"dg_id":self.id,
                    "dg_name":self.name,
                    "dg_synm":self.synonyms,
                    "dg_kingdom":self.kingdom,
                    "dg_superclass":self.superclass,
                    "dg_interactions":self.interaction,
                    "dg_ext_id":self.external_id,
                    "dg_pathways":self.pathways}
        return drug_dict

    def addTarget(self, feature_target):
        self.target.append(feature_target)


def main():
    global args
    args = parse_arg()
    xtree = ET.parse(args.input_xml)
    xroot = xtree.getroot()
    drugs = list(xroot)
    print("Total drugs in the database:", len(drugs))

    organism = 'Humans'

    drug_targets = []
    for i in tqdm(range(len(drugs))):
        drug = drugs[i]
        idDB = drug[0].text # Drug Bank ID

        for idx,feature in enumerate(drug):
            if 'name' in str(feature): # drug name
                drug_name = drug[idx].text

            if 'synonyms' in str(feature): # drug's synonyms
                drug_synm = ';'.join([synm.text \
                                        for synm in list(drug[idx])])

            if 'classification' in str(feature): #type of drug
                drug_class_kingdom = list(drug[idx])[2].text
                drug_class_superclass = list(drug[idx])[3].text

            if 'drug-interactions' in str(feature): #interaction other drugs
                drug_interaction = ';'.join([di[0].text
                                            for di in list(drug[idx])])

            if 'external-identifiers' in str(feature): #other drug's IDs
                aux = [ext_id[0].text + ":" + ext_id[1].text \
                                            for ext_id in list(drug[idx])]
                drug_external_id = ';'.join(aux)

            if 'pathways' in str(feature): #related pathways
                drug_pathway = ';'.join([pathway[1].text \
                                        for pathway in list(drug[idx])])

            if 'targets' in str(feature): #if polypeptide, drug's targets
                targets = list(drug[idx])

        # get all drug-related information in a dictionary
        drug_dict = {"id":idDB,
                    "name":drug_name,
                    "synm":drug_synm,
                    "kgd":drug_class_kingdom,
                    "sclass":drug_class_superclass,
                    "itrc":drug_interaction,
                    "ext_id":drug_external_id,
                    "pathways":drug_pathway}
        drug = Drug(drug_dict)

        # get information of polypeptide targets
        if len(targets) > 0:
            for target in targets:
                idx_pep = None
                # get indexes
                for idx,feature in enumerate(target): # check features of targets
                    if 'organism' in str(feature):
                        idx_org = idx
                    if 'name' in str(feature):
                        idx_name = idx
                    if 'actions' in str(feature):
                        idx_act = idx
                    if 'polypeptide' in str(feature):
                        idx_pep = idx

                # Get information for polypeptide
                if target[idx_org].text == organism:
                    target_name = target[idx_name].text
                    actions = ';'.join([action.text
                                        for action in list(target[idx_act])])
                    # Get information for polypeptide
                    if idx_pep is not None: #if there is polypeptide's info...
                        for idx,feature in enumerate(target[idx_pep]):
                            if 'gene-name' in str(feature):
                                gene_name = target[idx_pep][idx].text
                            if 'cellular-location' in str(feature):
                                cell_loc = target[idx_pep][idx].text
                            if 'external-identifiers' in str(feature):
                                for ext_id in list(target[idx_pep][idx]):
                                    if ext_id[0].text == "UniProtKB":
                                        uniprot = ext_id[1].text
                    else:
                        gene_name = None
                        action = None
                        cell_loc = None
                        uniprot = None

                    row = {
                            "dg_id":drug.id,
                            "dg_name":drug.name,
                            "dg_synm":drug.synonyms,
                            "dg_kingdom":drug.kingdom,
                            "dg_superclass":drug.superclass,
                            #"dg_interactions":drug.interaction,
                            #"dg_ext_id":drug.external_id,
                            "dg_pathways":drug.pathways,
                            "target_name":target_name,
                            "target_uniprot":uniprot,
                            "target_gene_name":gene_name,
                            "action":actions,
                            "cell_loc":cell_loc,
                            }

                    drug_targets.append(row)
        else:
            row = {
                    "dg_id":drug.id,
                    "dg_name":drug.name,
                    "dg_synm":drug.synonyms,
                    "dg_kingdom":drug.kingdom,
                    "dg_superclass":drug.superclass,
                    #"dg_interactions":drug.interaction,
                    #"dg_ext_id":drug.external_id,
                    "dg_pathways":drug.pathways,
                    "target_name":None,
                    "target_uniprot":None,
                    "target_gene_name":None,
                    "action":None,
                    "cell_loc":None,
                    }
            drug_targets.append(row)


    dt = pd.DataFrame.from_dict(drug_targets, orient='columns')
    dt.shape
    dt.to_csv(args.output, sep = ',', index = False)
    
if __name__=="__main__":
    main()

