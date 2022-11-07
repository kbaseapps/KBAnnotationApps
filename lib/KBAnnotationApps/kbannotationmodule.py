from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json, requests
import aiohttp, asyncio
from requests.exceptions import ConnectionError, HTTPError, RequestException
from python_graphql_client import GraphqlClient
import pandas as pd
import hashlib
import re
from kbbasemodules.basemodule import BaseModule

logger = logging.getLogger(__name__)

meta_data_query ="""{
    entries(entry_ids:["%s"]) {
        rcsb_id
        exptl {
            method
        }
        rcsb_primary_citation {
            title
            rcsb_authors
            journal_abbrev
            year
            pdbx_database_id_PubMed
        }
        polymer_entities {
            rcsb_id
            entity_poly {
                pdbx_seq_one_letter_code
                pdbx_strand_id
            }
            rcsb_entity_source_organism {
                ncbi_taxonomy_id
                ncbi_scientific_name
            }
            rcsb_polymer_entity_container_identifiers {
                reference_sequence_identifiers {
                    database_accession
                    database_name
                }
            }
            rcsb_polymer_entity {
                rcsb_ec_lineage {
                    id
                }
            }
            uniprots {
              rcsb_uniprot_protein {
                name {
                  value
                }
                ec {
                  number
                  provenance_code
                }
              }
            }
        }
        nonpolymer_entities {
          nonpolymer_comp {
            rcsb_chem_comp_synonyms {
              name
              provenance_source
            }
            rcsb_chem_comp_descriptor {
              InChI
              InChIKey
              SMILES
            }
          }
        }
    }
}
"""

class KBAnnotationModule(BaseModule):
    def __init__(self,name,ws_client,anno_api_client,working_dir,config):
        BaseModule.__init__(self,name,ws_client,working_dir,config)
        self.anno_api_client = anno_api_client
        self.graphqlClient = GraphqlClient(endpoint='https://data.rcsb.org/graphql')
        self.genome_info_hash = {}
        self.inchK_cpd_jsonObj = {}
        json_file_path = os.path.join(os.path.dirname(__file__), 'inchikey_cpd.json')
        with open(json_file_path) as DATA:
            self.inchK_cpd_jsonObj = json.load(DATA)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
    
    #Main function for annotating genomes with PDB
    def PDBAnnotation(self,params):
        self.initialize_call("PDBAnnotation",params,True)
        self.validate_args(params,["workspace","genome_ref"],{
            "suffix":".pdb",
            "similarity_threshold_type":"evalue",
            "similarity_threshold":0.00001,
            "return_data":False,
            "save_annotated_genome":True
        })
        sequence_list = self.genome_to_proteins(params["genome_ref"])
        query_table = self.query_rcsb_with_proteins(sequence_list,params["similarity_threshold_type"],params["similarity_threshold"])
        if params["save_annotated_genome"]:
            self.add_annotations_to_genome(params["genome_ref"],params["suffix"],query_table)
        output = self.build_report(query_table,params["genome_ref"])
        if params["return_data"]:
            output["data"] = query_table
        return output
    
    #Utility functions
    def genome_to_proteins(self,ref):
        output = self.get_object(ref,self.ws_id)
        self.genome_info_hash[ref] = output["info"]
        sequence_list = []
        for ftr in output["data"]["features"]:
            if "protein_translation" in ftr:
                sequence_list.append([ftr["id"],ftr["protein_translation"]])
        return sequence_list
    
    def query_rcsb_with_sequence(self,sequence,theshold_type,threshold_value,maxhits=500):
        #Creating query datastructure
        query_input = {
          "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
              "sequence_type": "protein",
              "value": sequence
            }
          },
          "return_type": "polymer_entity",
          "request_options": {
            "results_verbosity": "verbose",
            "paginate": {
              "start": 0,
              "rows": maxhits
            },
            "results_content_type": [
              "experimental"
            ],
            "sort": [
              {
                "sort_by": "score",
                "direction": "desc"
              }
            ],
            "scoring_strategy": "combined"
          }
        }
        #Setting similarity thresholds based on parameters
        if theshold_type == "evalue":
            query_input["query"]["parameters"]["evalue_cutoff"] = threshold_value
        else:
            query_input["query"]["parameters"]["identity_cutoff"] = threshold_value
        #Querying RCSB for IDs that match the input sequence
        results = {}
        try:
            reqH = requests.post('https://search.rcsb.org/rcsbsearch/v2/query',json=query_input)
            reqH.raise_for_status()
            results = json.loads(reqH.text)
        except (HTTPError, ConnectionError, RequestException) as e:
            logging.info(" _queryRCSB ERROR ".center(30, "-"))
            logging.info(f'Querying RCSB db with {query_input} had an Error: {e}')
            return {}
        except Exception as e:
            print(e)
            return {}
        formatted_results = {}
        if "result_set" in results:
            for item in results["result_set"]:
                if "identifier" in item:
                    if "services" in item and len(item["services"]) >= 1:
                        service = item["services"][0]
                        if "nodes" in service and len(service["nodes"]) >= 1:
                            node = service["nodes"][0]
                            if "match_context" in node and len(node["match_context"]) >= 1:
                                matchdata = node["match_context"][0]
                                formatted_results[item["identifier"]] = matchdata
                                #Example match data:
                                    # "sequence_identity": 0.521,
                                    # "evalue": 6.592e-58,
                                    # "bitscore": 211,
                                    # "alignment_length": 184,
                                    # "mismatches": 88,
                                    # "gaps_opened": 0,
                                    # "query_beg": 74,
                                    # "query_end": 257,
                                    # "subject_beg": 28,
                                    # "subject_end": 211,
                                    # "query_length": 426,
                                    # "subject_length": 382,
                                    # "query_aligned_seq": "DSFYMRKCVELAKRAIGCTSPNPMVGCVIVKDGDIVGQGFHPKAGQPHAEVFALRDAGELAENATAYVSLEPCNHYGRTPPCTEALIKAKVRRVVIGMVDPNPIVFSSGISRLKDAGIDVTVSVEEELCKKMNEGFIHRMLTGKPFLALRYSMSVNGCLLDKIGQGASDSGGYYSKLLQEYDAI",
                                    # "subject_aligned_seq": "DQYWMQQAIELAKRGLYSTKPNPNVGCVIVKDDQLIGEGFHPKAGQPHAEVFALRQAGEQAQGATAYVTLEPCAHYGRTPPCAEALVKAQVKKVVVACPDPNPLVAGKGVQILKNAGIEVEIGICEDLAAKLNQGFLKAMSTGMPYVRLKVASSLDGRTAMASGESKWITGSAARQDVQHWRAI"
        return formatted_results
    
    def query_rcsb_metadata_by_id(self,idlist,bundle_size=100):
        count = 0
        current_bundle = []
        metadata_hash = {}
        for id in idlist:
            count += 1
            current_bundle.append(id)
            if len(current_bundle) >= bundle_size or count >= len(idlist):
                query_string = meta_data_query % '", "'.join(current_bundle)
                raw_data = {}
                try:
                    #graphql_ret = self.__graphqlClient.execute(query=queryString)
                    # asyncio.run() is available ONLY for python 3.7 or newer
                    # graphql_ret = asyncio.run(self.__graphqlClient.execute_async(query=queryString))
                    evt_loop = asyncio.get_event_loop()
                    raw_data = evt_loop.run_until_complete(
                                  self.graphqlClient.execute_async(query=query_string))
                    if 'errors' in raw_data:
                        raise ConnectionError(graphql_ret['errors'][0]['message'])
                except (aiohttp.ClientConnectionError, asyncio.TimeoutError):
                    logging.info(" _queryGraphql ERROR ".center(30, "-"))
                    logging.info(f'Connecting to RCSB GraphQL host at https://data.rcsb.org/graphql errored.')
                    return {'data': None}
                except (HTTPError, ConnectionError, RequestException) as e:
                    # not raising error to allow continue with other chunks
                    logging.info(f'Querying RCSB GraphQL had a Connection Error:************\n {e}.\n'
                                 'Or database connection request had no response!')
                    return {'data': None}
                except (RuntimeError, TypeError, KeyError, ValueError) as e:
                    err_msg = f'Querying RCSB errored with message: {e.message} and data: {e.data}'
                    raise ValueError(err_msg)
                if raw_data.get('data', None) and raw_data['data'].get('entries', None):
                    # short-naming the long rcsb data attribute strings
                    for entry in raw_data["data"]["entries"]:
                        id = entry['rcsb_id']
                        #Setting primary metadata for structure
                        metadata_hash[id] = {
                            "reference":None,
                            "methods":[],
                            "compounds":[],
                            "proteins":{}
                        }
                        if 'rcsb_primary_citation' in entry and entry["rcsb_primary_citation"]:
                            pubmedid = entry["rcsb_primary_citation"].get('pdbx_database_id_PubMed', '')
                            title = entry["rcsb_primary_citation"].get('title', '')
                            journal = entry["rcsb_primary_citation"].get('journal_abbrev', '')
                            year = entry["rcsb_primary_citation"].get('year', 'NA')
                            authorlist = entry["rcsb_primary_citation"].get('rcsb_authors', [])
                            author = ""
                            if authorlist and len(authorlist) > 0:
                                author = authorlist[0]
                            metadata_hash[id]["reference"] = [pubmedid,title,journal,author,str(year)]
                        for exp in entry.get('exptl', []):
                            if exp.get('method', ''):
                                metadata_hash[id]["methods"].append(exp['method'])
                        #Adding subproteins from this structure complex
                        for pe in entry['polymer_entities']:
                            protein_data = {
                                "strands":pe['entity_poly']['pdbx_strand_id'].split(','),
                                'taxonomy':[],
                                'uniprotID':[],
                                'ec_numbers':[],
                                'uniprot_name':[],
                                'uniprot_ec':[]
                            }
                            if (pe.get('rcsb_entity_source_organism', None)):
                                for srco in pe['rcsb_entity_source_organism']:
                                    protein_data['taxonomy'].append((srco.get('ncbi_taxonomy_id', ''), srco.get('ncbi_scientific_name', '')))
                            if (pe.get('rcsb_polymer_entity_container_identifiers', None) and
                                    pe['rcsb_polymer_entity_container_identifiers'].get('reference_sequence_identifiers', None)):
                                protein_data['ref_sequence_ids'] = pe['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers']
                                for rsid in protein_data['ref_sequence_ids']:
                                    protein_data['uniprotID'].append(rsid.get('database_accession', ''))
                            if (pe.get('rcsb_polymer_entity', None) and pe['rcsb_polymer_entity'].get('rcsb_ec_lineage', None)):
                                for idic in pe['rcsb_polymer_entity']['rcsb_ec_lineage']:
                                    if idic['id'] and (idic['id'] not in protein_data['ec_numbers']):
                                        protein_data['ec_numbers'].append(idic['id'])
                            if pe.get('uniprots', None):
                                for unp in pe['uniprots']:
                                    if unp.get('rcsb_uniprot_protein', None):
                                        uniprot_prot = unp['rcsb_uniprot_protein']
                                        if uniprot_prot.get('name', None):
                                            protein_data['uniprot_name'].append(
                                                uniprot_prot['name'].get('value', ''))
                                        if uniprot_prot.get('ec', None):
                                            protein_data['uniprot_ec'].extend(uniprot_prot['ec'])
                            metadata_hash[id]["proteins"][pe['rcsb_id']] = protein_data
                        #Adding compound metadata
                        if entry.get('nonpolymer_entities', None):
                            for npe in entry.get('nonpolymer_entities', []):
                                if npe.get('nonpolymer_comp', None):
                                    compound_data = {
                                        'pdb_ref_name':''
                                    }
                                    if npe['nonpolymer_comp'].get('rcsb_chem_comp_synonyms', []):
                                        for synm in npe['nonpolymer_comp'].get('rcsb_chem_comp_synonyms', []):
                                            # assuming at least one 'PDB Reference Data' is always available
                                            if synm['provenance_source'] == 'PDB Reference Data':
                                                compound_data['pdb_ref_name'] = synm['name']
                                                break
                                    if npe['nonpolymer_comp'].get('rcsb_chem_comp_descriptor', {}):
                                        for desType in ('InChI', 'InChIKey', 'SMILES'):
                                            if npe['nonpolymer_comp']['rcsb_chem_comp_descriptor'].get(desType, None):
                                                dt_val = npe['nonpolymer_comp']['rcsb_chem_comp_descriptor'][desType]
                                                compound_data[desType] = dt_val
                                        inchiK = compound_data.get('InChIKey', '')
                                        ik = inchiK[:len(inchiK) - 2]
                                        for k in self.inchK_cpd_jsonObj.keys():
                                            if ik in k:
                                                compound_data["mscpd"] = self.inchK_cpd_jsonObj[k]
                                                break                                    
                                    metadata_hash[entry['rcsb_id']]['compounds'].append(compound_data)
                current_bundle = []
                print("Done querying metadata for ",count," of ",len(idlist))
        return metadata_hash
    
    def query_rcsb_with_proteins(self,sequence_list,cutoff_type="evalue",threhold=0.00001,bundle_size=100):
        output_table = {"id":[],"rcsbid":[],"method":[],"strand":[],"similarity":[],"taxonomy":[],"name":[],"components":[],"rcsbec":[],"uniprotec":[],"uniprotID":[],"references":[]}
        #First we get all the hits in RCSB one at a time (maybe bulk this later?)
        all_hits = {}
        distinct_ids = {}
        count = 0
        for item in sequence_list:
            count += 1
            hits = self.query_rcsb_with_sequence(item[1],cutoff_type,threhold)
            for id, data in hits.items():
                if item[0] not in all_hits:
                    all_hits[item[0]] = {}
                all_hits[item[0]][id] = data
                array = id.split("_")
                distinct_ids[array[0]] = 1
            if count % 100 == 0:
                print("Done querying sequences for ",count," of ",len(sequence_list))
        print(len(all_hits)," distinct PDB IDs hit across ",len(sequence_list)," input features!")
        #Now we bulk query for metadata about all of our hits
        self.print_json_debug_file("MetadataList.json",list(distinct_ids.keys()))
        metadata_hash = self.query_rcsb_metadata_by_id(list(distinct_ids.keys()))
        #Now we build our results table
        for item in sequence_list:
            if item[0] in all_hits:
                retained_hits = []
                ecnumbers = {}
                for hit in all_hits[item[0]]:
                    array = hit.split("_")
                    struct_row = metadata_hash[array[0]]
                    prot_row = metadata_hash[array[0]]["proteins"][hit]
                    if prot_row["ec_numbers"] or prot_row["uniprot_ec"]:
                        ec_hash = {}
                        if prot_row["ec_numbers"]:
                            longest_ec = ""
                            longest_ec_size = 0
                            for ec in prot_row["ec_numbers"]:
                                array = ec.split(".")
                                if len(array) == 4:
                                    ec_hash[ec] = "rcsbec"
                                    longest_ec = None
                                if longest_ec and len(array) > longest_ec_size:
                                    longest_ec_size = len(array)
                                    longest_ec = ec
                            if longest_ec:
                                ec_hash[longest_ec] = "rcsbec"
                        if prot_row["uniprot_ec"]:
                            for ec in prot_row["uniprot_ec"]:
                                ec_hash[ec["number"]] = "uniprot_ec"
                        for ec in ec_hash:
                            if ec not in ecnumbers or all_hits[item[0]][ecnumbers[ec]]["evalue"] > all_hits[item[0]][hit]["evalue"]:
                                ecnumbers[ec] = hit
                    else:
                        retained_hits.append(hit)
                for ec in ecnumbers:
                    if ecnumbers[ec] not in retained_hits:
                        retained_hits.append(ecnumbers[ec])
                for hit in retained_hits:
                    output_table["id"].append(item[0])
                    output_table["rcsbid"].append(hit)
                    if len(prot_row["uniprot_name"]) > 0:
                        output_table["name"].append(prot_row["uniprot_name"][0])
                    else:
                        output_table["name"].append("")
                    if len(struct_row["methods"]) > 0:
                        output_table["method"].append(struct_row["methods"][0])
                    else:
                        output_table["method"].append("")
                    output_table["strand"].append(prot_row["strands"])
                    output_table["similarity"].append([all_hits[item[0]][hit]["evalue"],all_hits[item[0]][hit]["sequence_identity"]])
                    output_table["taxonomy"].append(prot_row["taxonomy"])
                    output_table["components"].append(struct_row["compounds"])
                    output_table["rcsbec"].append(prot_row["ec_numbers"])
                    output_table["uniprotec"].append(prot_row["uniprot_ec"])
                    output_table["uniprotID"].append(prot_row["uniprotID"])
                    output_table["references"].append(struct_row["reference"])
        return output_table

    def add_annotations_to_genome(self,ref,suffix,pdb_query_output):
        ontology_inputs = {"RCSBID":{}}
        self.print_json_debug_file("RCSBQueryResults.json",pdb_query_output)
        for count, geneid in enumerate(pdb_query_output["id"]):
            if geneid not in ontology_inputs["RCSBID"]:
                ontology_inputs["RCSBID"][geneid] = []
            ontology_inputs["RCSBID"][geneid].append({
                "term":"RCSBID:"+pdb_query_output["rcsbid"][count],
                "name":pdb_query_output["name"][count]+";"+pdb_query_output["method"][count]+";"+",".join(pdb_query_output["strand"][count])+";Evalue="+str(pdb_query_output["similarity"][count][0])
            })
            if pdb_query_output["taxonomy"][count]:
                if "TAXID" not in ontology_inputs:
                    ontology_inputs["TAXID"] = {}
                if geneid not in ontology_inputs["TAXID"]:
                    ontology_inputs["TAXID"][geneid] = []
                ontology_inputs["TAXID"][geneid].append({
                    "term":"TAXID:"+str(pdb_query_output["taxonomy"][count][0][0]),
                    "name":pdb_query_output["taxonomy"][count][0][1]+";"+pdb_query_output["rcsbid"][count]
                })
            if pdb_query_output["uniprotID"][count]:
                if "UNIPROT" not in ontology_inputs:
                    ontology_inputs["UNIPROT"] = {}
                if geneid not in ontology_inputs["UNIPROT"]:
                    ontology_inputs["UNIPROT"][geneid] = []
                for item in pdb_query_output["uniprotID"][count]:
                    ontology_inputs["UNIPROT"][geneid].append({
                        "term":"UNIPROT:"+item,
                        "name":pdb_query_output["name"][count]+";"+pdb_query_output["rcsbid"][count]
                    })
            if pdb_query_output["rcsbec"][count] or pdb_query_output["uniprotec"][count]:
                if "EC" not in ontology_inputs:
                    ontology_inputs["EC"] = {}
                if geneid not in ontology_inputs["EC"]:
                    ontology_inputs["EC"][geneid] = []
                ec_hash = {}
                if pdb_query_output["rcsbec"][count]:
                    for item in pdb_query_output["rcsbec"][count]:
                        array = item.split(".")
                        if len(array) == 4:
                            ec_hash[item] = ";"+pdb_query_output["rcsbid"][count]
                if pdb_query_output["uniprotec"][count]:
                    for item in pdb_query_output["uniprotec"][count]:
                        if item["number"] in ec_hash:
                            ontology_inputs["EC"][geneid].append({
                                "term":"EC:"+item["number"],
                                "suffix":ec_hash[item["number"]]+";"+pdb_query_output["uniprotID"][count][0]
                            })
                        else:
                            ontology_inputs["EC"][geneid].append({
                                "term":"EC:"+item["number"],
                                "suffix":";"+pdb_query_output["uniprotID"][count][0]
                            })
            if pdb_query_output["uniprotec"][count]:
                if "EC" not in ontology_inputs:
                    ontology_inputs["EC"] = {}
                if geneid not in ontology_inputs["EC"]:
                    ontology_inputs["EC"][geneid] = []
                
            if pdb_query_output["references"][count]:
                if pdb_query_output["references"][count][0] != None:
                    if "PUBMED" not in ontology_inputs:
                        ontology_inputs["PUBMED"] = {}
                    if geneid not in ontology_inputs["PUBMED"]:
                        ontology_inputs["PUBMED"][geneid] = []
                    ontology_inputs["PUBMED"][geneid].append({
                        "term":"PUBMED:"+str(pdb_query_output["references"][count][0]),
                        "name":pdb_query_output["rcsbid"][count]+":"+pdb_query_output["references"][count][3]+" et al."+pdb_query_output["references"][count][1]+"."+pdb_query_output["references"][count][2]+"("+pdb_query_output["references"][count][2]+")"
                    })
                else:
                    if "REF" not in ontology_inputs:
                        ontology_inputs["REF"] = {}
                    if geneid not in ontology_inputs["REF"]:
                        ontology_inputs["REF"][geneid] = []
                    ontology_inputs["REF"][geneid].append({
                        "term":"REF:"+pdb_query_output["references"][count][1],
                        "name":pdb_query_output["rcsbid"][count]+":"+pdb_query_output["references"][count][3]+" et al."+pdb_query_output["references"][count][2]+"("+pdb_query_output["references"][count][2]+")"
                    })
            if pdb_query_output["components"][count]:
                for item in pdb_query_output["components"][count]:
                    if "mscpd" in item:
                        m = re.search('(cpd\d+)\s\((.+)\)', item["mscpd"])
                        if "MSCPD" not in ontology_inputs:
                            ontology_inputs["MSCPD"] = {}
                        if geneid not in ontology_inputs["MSCPD"]:
                            ontology_inputs["MSCPD"][geneid] = []
                        ontology_inputs["MSCPD"][geneid].append({
                            "term":"MSCPD:"+m[1],
                            "name":m[2]+";cocrystalized in "+pdb_query_output["rcsbid"][count]
                        })
                    elif "InChI" in item:
                        if "InChI" not in ontology_inputs:
                            ontology_inputs["InChI"] = {}
                        if geneid not in ontology_inputs["InChI"]:
                            ontology_inputs["InChI"][geneid] = []
                        ontology_inputs["InChI"][geneid].append({
                            "term":"InChI:"+item["InChI"],
                            "name":m[1]+";cocrystalized in "+pdb_query_output["rcsbid"][count]
                        })
        anno_api_input = {
            "input_ref":ref,
            "output_name":self.genome_info_hash[ref][1]+suffix,
            "output_workspace":self.ws_id,
            "overwrite_matching":1,
            "save":1,
            "provenance":self.provenance(),
            "events":[]
        }
        for ontology in ontology_inputs.keys():
            anno_api_input["events"].append({
                "ontology_id":ontology,
                "method":self.name+"."+self.method,
                "method_version":self.config["version"],
                "timestamp":self.timestamp,
                "ontology_terms":ontology_inputs[ontology]
            })
        anno_api_output = self.anno_api_client.add_annotation_ontology_events(anno_api_input)
        self.obj_created.append({"ref":anno_api_output["output_ref"],"description":"Saving PDB annotation for "+self.genome_info_hash[ref][1]})
        return anno_api_output
            
    def build_report(self,data,ref):        
        table = pd.DataFrame(data)
        for index, row in table.iterrows():
            row["id"] = '<a href="https://appdev.kbase.us/#dataview/'+ref+'?sub=Feature&subid='+row["id"]+'" target="_blank">'+row["id"]+'</a>'
            array = row["rcsbid"].split("_")
            row["rcsbid"] = '<a href="https://www.rcsb.org/3d-view/'+array[0]+'" target="_blank">'+row["rcsbid"]+'</a>'
            row["strand"] = ", ".join(row["strand"])
            newdata = ""
            inchidata = ""
            for item in row["components"]:
                if "mscpd" in item:
                    m = re.search('(cpd\d+)\s\((.+)\)', item["mscpd"])
                    if len(newdata) > 0:
                        newdata += "<br>"
                    newdata += '<a href=" https://modelseed.org/biochem/compounds/'+m[1]+'" target="_blank">'+m[2]+'</a>'
                elif "InChI" in item:
                    if len(inchidata) > 0:
                        inchidata += "<br>"
                    inchidata += item["InChI"]  
            if len(newdata) > 0 and len(inchidata) > 0:
                row["components"] = newdata+"<br>"+inchidata
            elif len(newdata) == 0 and len(inchidata) > 0:
                row["components"] = inchidata
            elif len(newdata) > 0 and len(inchidata) == 0:
                row["components"] = newdata
            else:
                row["components"] = ""
            newsim = "Identity: "+str(row["similarity"][1])+"<br>Evalue: "+str(row["similarity"][0])
            row["similarity"] = newsim
            newec = ""
            echash = {}
            for ec in row["rcsbec"]:
                array = ec.split(".")
                if len(array) == 4:
                    echash[ec] = ["RCSB"]
            for item in row["uniprotec"]:
                if "number" in item:
                    if item["number"] not in echash:
                        echash[item["number"]] = []
                    echash[item["number"]].append("UniProt")
            for ec in echash:
                if len(newec) > 0:
                    newec += "<br>"
                newec += '<a href="https://www.kegg.jp/entry/'+ec+'" target="_blank">'+ec+"("+"/".join(echash[ec])+')</a>'
            row["rcsbec"] = newec
            newuniprot = ""
            for id in row["uniprotID"]:
                if len(newuniprot) > 0:
                    newuniprot += "<br>"
                newuniprot += '<a href="https://www.uniprot.org/uniprotkb/'+id+'/entry" target="_blank">'+id+'</a>'
            row["uniprotID"] = newuniprot
            refdata = ""
            if row["references"]:
                if row["references"][0]:
                    refdata = '<a href="https://pubmed.ncbi.nlm.nih.gov/'+str(row["references"][0])+'/" target="_blank">'+str(row["references"][0])+'</a>'
                elif row["references"][4] != "None":
                    refdata = row["references"][1]+". "+row["references"][2]+" ("+row["references"][4]+")"
                else:
                    refdata = row["references"][1]+". "+row["references"][2]
                row["references"] = refdata
            taxonomy = ""
            for item in row["taxonomy"]:
                if len(taxonomy) > 0:
                    taxonomy += "<br>"
                taxonomy += '<a href="https://narrative.kbase.us/#dataview/12570/'+str(item[0])+'/" target="_blank">'+item[1]+'</a>'
            row["taxonomy"] = taxonomy
        #columns=column_list
        table = table.drop(columns=['uniprotec'])
        html_data = f"""
    <html>
    <header>
        <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
    </header>
    <body>
    {table.to_html(escape=False,notebook=False,table_id="table",index=False,justify="left")}
    <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready( function () {{
            $('#table').DataTable({{
                // paging: false,    
                // scrollY: 400,
            }});
        }});
    </script>
    </body>
    </html>
    """
        report_name = str(uuid.uuid4())
        html_report_folder = os.path.join(self.working_dir, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        with open(os.path.join(html_report_folder, 'index.html'), 'w') as f:
            f.write(html_data)
        return {
            'file_path':os.path.join(html_report_folder, 'index.html'),
            'report_params':{
                'objects_created': self.obj_created,
                'workspace_name': self.ws_name,
                'html_links': [{
                    'name' : 'index.html',
                    'shock_id': None
                }],
                'direct_html_link_index': 0,
                'html_window_height': 700,
                'report_object_name': report_name
            }
        }