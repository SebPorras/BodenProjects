from tkinter import N
from unicodedata import name
import xmlschema
import requests
import regex as re
import pickle

# define xml schema

def uniprot_retrieve(ids: list):
    
    entries = dict()

    for id in ids:
        schema = xmlschema.XMLSchema('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd')
        url    = f'https://www.uniprot.org/uniprotkb/{id}.xml'
        r      = requests.get(url)
        responseBody = r.text
        entry_dict = schema.to_dict(responseBody)
        entries[id] = entry_dict

    return entries

    # to bring in lineage information ( or anything else ...follow the xml nested structure)
    #content['organism']['lineage']['taxon']


#this returns a dictionary with the following keys
#'@dataset', '@created', '@modified', '@version', 'accession',
#'name', 'protein', 'gene', 'organism', 'reference', 'comment',
#'dbReference', 'proteinExistence', 'keyword', 'feature', 'evidence',
#'sequence'


#features 
#'@evidence': [1], 'location': {'position': {'@position': 227, '@status': 'certain'}}, 
#'ligand': {'name': 'Mg(2+)', 'dbReference': {'@type': 'ChEBI', '@id': 'CHEBI:18420'}, 'label': '2'}}



test = uniprot_retrieve(['P29107'])

#can access the entry with syntax below 

def record_annotations(entry_xml: dict):
    
    annot_dict = {}

    for entry in entry_xml.keys():

        annotations = {"cofactor_binding": {}, "domains": {}}

        info = entry_xml[entry]['entry'][0]

       
        features = info['feature']
        comments = info['comment']

        for x in comments:
            if x['@type'] == 'cofactor':
                cofactor = x['cofactor'][0]
                
        for x in features:
            if x['@type'] == 'domain': 

                start = x['location']['begin']['@position']
                end = x['location']['end']['@position']
                
                annotations['domains'][x['@description']] = {'start': start, 'end': end}

               

            elif x['@type'] == 'binding site':
                
                print(cofactor)

                print(x['@type'])
                
                print('Location')
                print(x['location'])

                print('ligand')
                print(x['ligand'])


record_annotations(test)
    