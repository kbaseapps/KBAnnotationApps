from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import pandas as pd
import hashlib
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
        self.validate_args(params,["workspace","genome_refs"],{
            "suffix":".pdb",
            "similarity_threshold_type":"evalue",
            "similarity_threshold":0.00001,
            "return_data":False
        })
        tables = {}
        for ref in params["genome_refs"]:
            sequence_list = self.genome_to_proteins(ref)
            query_table = self.query_rcsb_with_proteins (sequence_list,params["similarity_threshold_type"],params["similarity_threshold"])
            self.add_annotations_to_genome(ref,params["suffix"],query_table)
            tables[self.genome_info_hash[ref][1]] = pd.DataFrame(query_table)
        output = self.build_report(tables)
        if params["return_data"]:
            output["data"] = data
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
    
    def query_rcsb_with_proteins(self,sequence_list,cutoff_type="evalue",threhold=0.00001):
        output_table = {"id":[],"rcsbid":[],"method":[],"strand":[],"similarity":[],"taxonomy":[],"name":[],"components":[],"rcsbec":[],"uniprotec":[],"uniprotID":[],"references":[]}
        seq_hash = {}
        query_rcsb_input = {
            "sequence_strings":[],
            "workspace_name":self.ws_name
        }
        if cutoff_type == "evalue":
            query_rcsb_input["evalue_cutoff"] = threhold
        else:
            query_rcsb_input["identity_cutoff"] = threhold
        for item in sequence_list:
            query_rcsb_input["sequence_strings"].append(item[1])
            seq_hash[item[0]] = hashlib.md5(item[1].encode()).hexdigest()
        pdb_query_output = self.pdb_query_client.query_rcsb_structures(query_rcsb_input)
        for item in sequence_list:
            row = {
                "rcsbid": "1A49_1",
                "name": ["Carbonic anhydrase 2"],
                "pdbx_sequence": "SKSHSEAGSAFIQTQQLHAAMADTFLEHMCRLDIDSAPITARNTGIICTIGPASRSVETLKEMIKSGMNVARMNFSHGTHEYHAETIKNVRTATESFASDPILYRPVAVALDTKGPEIRTGLIKGSGTAEVELKKGATLKITLDNAYMEKCDENILWLDYKNICKVVDVGSKVYVDDGLISLQVKQKGPDFLVTEVENGGFLGSKKGVNLPGAAVDLPAVSEKDIQDLKFGVEQDVDMVFASFIRKAADVHEVRKILGEKGKNIKIISKIENHEGVRRFDEILEASDGIMVARGDLGIEIPAEKVFLAQKMIIGRCNRAGKPVICATQMLESMIKKPRPTRAEGSDVANAVLDGADCIMLSGETAKGDYPLEAVRMQHLIAREAEAAMFHRKLFEELARSSSHSTDLMEAMAMGSVEASYKCLAAALIVLTESGRSAHQVARYRPRAPIIAVTRNHQTARQAHLYRGIFPVVCKDPVQEAWAEDVDLRVNLAMNVGKARGFFKKGDVVIVLTGWRPGSGFTNTMRVVPVP",
                "rcsbec": ["2", "2.7", "2.7.1", "2.7.1.40", "2.7.11", "2.7.11.1", "2.7.10", "2.7.10.2"],
                "uniprotec": [{"number": "4.2.1.1", "provenance_code": "up:Protein.up:recommendedName"}, {"number": "4.2.1.69", "provenance_code": "up:Protein.up:alternativeName"}],
                "identity": 1.0,
                "uniprotID": ["P11974"],
                "pdbx_strand_id": ["A", "B", "C", "D", "E", "F", "G", "H"],
                "taxonomy": [9986, "Oryctolagus cuniculus"],
                "evalue": 0.001,
                "method": ["X-RAY DIFFRACTION"],
                "components": {"InChIKey": ["NPYPAHLBTDXSSS-UHFFFAOYSA-N", "MUBZPKHOEPUJKR-UHFFFAOYSA-L", "JLVVSXFLKOJNIY-UHFFFAOYSA-N", "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"]},
                "references": [1910042, "Structure of the bis(Mg2+)-ATP-oxalate complex of the rabbit muscle pyruvate kinase at 2.1 A resolution: ATP binding over a barrel.", "Biochemistry", "Larsen, T.M.", "1998"]
            }
            #if seq_hash[item[0]] in pdb_query_output:
            #for row in pdb_query_output[seq_hash[item[0]]]:
            output_table["id"].append(item[0])
            output_table["rcsbid"].append(row["rcsbid"])
            output_table["name"].append(row["name"][0])
            output_table["method"].append(row["method"][0])
            output_table["strand"].append(row["pdbx_strand_id"])
            output_table["similarity"].append(row["evalue"])
            output_table["taxonomy"].append(row["taxonomy"])
            output_table["components"].append(row["components"]["InChIKey"])
            output_table["rcsbec"].append(row["rcsbec"])
            output_table["uniprotec"].append(row["uniprotec"])
            output_table["uniprotID"].append(row["uniprotID"])
            output_table["references"].append(row["references"])
        return output_table

    def add_annotations_to_genome(self,ref,suffix,pdb_query_output):
        ontology_inputs = {"RCSBID":{}}
        for count, geneid in enumerate(pdb_query_output["id"]):
            if geneid not in ontology_inputs["RCSBID"]:
                ontology_inputs["RCSBID"][geneid] = []
            ontology_inputs["RCSBID"][geneid].append({
                "term":"RCSBID:"+pdb_query_output["rcsbid"][count],
                "name":pdb_query_output["name"][count]+";"+pdb_query_output["method"][count]+";"+",".join(pdb_query_output["strand"][count])+";Evalue="+str(pdb_query_output["similarity"][count])
            })
            if pdb_query_output["taxonomy"][count]:
                if "TAXID" not in ontology_inputs:
                    ontology_inputs["TAXID"] = {}
                if id not in ontology_inputs["TAXID"]:
                    ontology_inputs["TAXID"][geneid] = []
                ontology_inputs["TAXID"][geneid].append({
                    "term":"TAXID:"+str(pdb_query_output["taxonomy"][count][0]),
                    "name":pdb_query_output["taxonomy"][count][1]+";"+pdb_query_output["rcsbid"][count]
                })
            if pdb_query_output["uniprotID"][count]:
                if "UNIPROT" not in ontology_inputs:
                    ontology_inputs["UNIPROT"] = {}
                if id not in ontology_inputs["UNIPROT"]:
                    ontology_inputs["UNIPROT"][geneid] = []
                for item in pdb_query_output["uniprotID"][count]:
                    ontology_inputs["UNIPROT"][geneid].append({
                        "term":"UNIPROT:"+item,
                        "name":pdb_query_output["name"][count]+";"+pdb_query_output["rcsbid"][count]
                    })
            if pdb_query_output["rcsbec"][count]:
                if "EC" not in ontology_inputs:
                    ontology_inputs["EC"] = {}
                if id not in ontology_inputs["EC"]:
                    ontology_inputs["EC"][geneid] = []
                for item in pdb_query_output["rcsbec"][count]:
                    array = item.split(".")
                    if len(array) == 4:
                        ontology_inputs["EC"][geneid].append({
                            "term":"EC:"+item,
                            "suffix":";"+pdb_query_output["rcsbid"][count]
                        })
            if pdb_query_output["uniprotec"][count]:
                if "EC" not in ontology_inputs:
                    ontology_inputs["EC"] = {}
                if id not in ontology_inputs["EC"]:
                    ontology_inputs["EC"][geneid] = []
                for item in pdb_query_output["uniprotec"][count]:
                    ontology_inputs["EC"][geneid].append({
                        "term":"EC:"+item["number"],
                        "suffix":";"+pdb_query_output["uniprotID"][count][0]
                    })
            if pdb_query_output["references"][count]:
                if "PUBMED" not in ontology_inputs:
                    ontology_inputs["PUBMED"] = {}
                if id not in ontology_inputs["PUBMED"]:
                    ontology_inputs["PUBMED"][geneid] = []
                ontology_inputs["PUBMED"][geneid].append({
                    "term":"PUBMED:"+str(pdb_query_output["references"][count][0]),
                    "name":pdb_query_output["rcsbid"][count]+":"+pdb_query_output["references"][count][3]+" et al."+pdb_query_output["references"][count][1]+"."+pdb_query_output["references"][count][2]+"("+pdb_query_output["references"][count][2]+")"
                })
            if pdb_query_output["components"][count]:
                if "InChIKey" not in ontology_inputs:
                    ontology_inputs["InChIKey"] = {}
                if id not in ontology_inputs["InChIKey"]:
                    ontology_inputs["InChIKey"][geneid] = []
                for item in pdb_query_output["components"][count]:
                    ontology_inputs["InChIKey"][geneid].append({
                        "term":"InChIKey:"+item,
                        "name":item+";cocrystalized in "+pdb_query_output["rcsbid"][count]
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
            
    def build_report(self,tables):
        genomes = tables.keys()
        table = tables[genomes[0]]
        #columns=column_list
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
            'data':table,
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