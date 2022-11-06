from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import pandas as pd
import hashlib
import re
from kbbasemodules.basemodule import BaseModule

logger = logging.getLogger(__name__)

class KBAnnotationModule(BaseModule):
    def __init__(self,name,ws_client,anno_api_client,pdb_query_client,working_dir,config):
        BaseModule.__init__(self,name,ws_client,working_dir,config)
        self.pdb_query_client = pdb_query_client
        self.anno_api_client = anno_api_client
        self.genome_info_hash = {}
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
            "bundle_size":20,
            "save_annotated_genome":False
        })
        sequence_list = self.genome_to_proteins(params["genome_ref"])
        query_table = self.query_rcsb_with_proteins(sequence_list,params["similarity_threshold_type"],params["similarity_threshold"],params["bundle_size"])
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
    
    def query_rcsb_with_proteins(self,sequence_list,cutoff_type="evalue",threhold=0.00001,bundle_size=10):
        output_table = {"id":[],"rcsbid":[],"method":[],"strand":[],"similarity":[],"taxonomy":[],"name":[],"components":[],"rcsbec":[],"uniprotec":[],"uniprotID":[],"references":[]}
        current_bundle = []
        count = 0
        for item in sequence_list:
            count += 1
            current_bundle.append(item)
            if len(current_bundle) >= bundle_size or item == sequence_list[-1]:
                query_rcsb_input = {
                    "sequence_strings":[],
                    "workspace_name":self.ws_name
                }
                if cutoff_type == "evalue":
                    query_rcsb_input["evalue_cutoff"] = threhold
                else:
                    query_rcsb_input["identity_cutoff"] = threhold
                for item in current_bundle:
                    query_rcsb_input["sequence_strings"].append(item[1])
                pdb_query_output = self.pdb_query_client.query_rcsb_annotations(query_rcsb_input)
                self.print_json_debug_file("QueryResults.json",pdb_query_output)
                for item in current_bundle:
                    if item[1] in pdb_query_output:
                        for row in pdb_query_output[item[1]]:
                            output_table["id"].append(item[0])
                            output_table["rcsbid"].append(row["rcsbid"])
                            if len(row["name"]) > 0:
                                output_table["name"].append(row["name"][0])
                            else:
                                output_table["name"].append(None)
                            if len(row["method"]) > 0:
                                output_table["method"].append(row["method"][0])
                            else:
                                output_table["method"].append(None)
                            output_table["strand"].append(row["pdbx_strand_id"])
                            output_table["similarity"].append([row["evalue"],row["identity"]])
                            output_table["taxonomy"].append(row["taxonomy"])
                            output_table["components"].append(row["components"])
                            output_table["rcsbec"].append(row["rcsbec"])
                            output_table["uniprotec"].append(row["uniprotec"])
                            output_table["uniprotID"].append(row["uniprotID"])
                            output_table["references"].append(row["references"])
                current_bundle = []
                print("Done querying ",count," of ",len(sequence_list))
        return output_table

    def add_annotations_to_genome(self,ref,suffix,pdb_query_output):
        ontology_inputs = {"RCSBID":{}}
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
                    if re.search('(cpd\d+)\s\((.+)\)', item):
                        m = re.search('(cpd\d+)\s\((.+)\)', item)
                        if "MSCPD" not in ontology_inputs:
                            ontology_inputs["MSCPD"] = {}
                        if geneid not in ontology_inputs["MSCPD"]:
                            ontology_inputs["MSCPD"][geneid] = []
                        ontology_inputs["MSCPD"][geneid].append({
                            "term":"MSCPD:"+m[1],
                            "name":m[2]+";cocrystalized in "+pdb_query_output["rcsbid"][count]
                        })
                    elif re.search('^InChI', item):
                        if "InChI" not in ontology_inputs:
                            ontology_inputs["InChI"] = {}
                        if geneid not in ontology_inputs["InChI"]:
                            ontology_inputs["InChI"][geneid] = []
                        ontology_inputs["InChI"][geneid].append({
                            "term":"InChI:"+m[1],
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
                if re.search('(cpd\d+)\s\((.+)\)', item):
                    m = re.search('(cpd\d+)\s\((.+)\)', item)
                    if len(newdata) > 0:
                        newdata += "<br>"
                    newdata += '<a href=" https://modelseed.org/biochem/compounds/'+m[1]+'" target="_blank">'+m[2]+'</a>'
                elif re.search('^InChI', item):
                    if len(inchidata) > 0:
                        inchidata += "<br>"
                    inchidata += item  
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
            if row["references"][0]:
                refdata = '<a href="https://pubmed.ncbi.nlm.nih.gov/'+row["references"][0]+'/" target="_blank">'+row["references"][0]+'</a>'
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