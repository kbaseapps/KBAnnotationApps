# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
sys.path.append("/deps/KBBaseModules/")
import sys
import json
from os.path import exists
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.ProteinStructureUtilsClient import ProteinStructureUtils
from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api
from KBAnnotationApps.kbannotationmodule import KBAnnotationModule
#END_HEADER


class KBAnnotationApps:
    '''
    Module Name:
    KBAnnotationApps

    Module Description:
    A KBase module: KBAnnotationApps
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def save_report_to_kbase(self,output):
        report_shock_id = self.dfu.file_to_shock({'file_path': output["file_path"],'pack': 'zip'})['shock_id']
        output["report_params"]["html_links"][0]["shock_id"] = report_shock_id
        repout = self.kbreport.create_extended_report(output["report_params"])
        return {"report_name":output["report_params"]["report_object_name"],"report_ref":repout["ref"],'workspace_name':self.api.ws_name}
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config["version"] = self.VERSION
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.wsclient = Workspace(self.config["workspace-url"], token=self.token)
        self.kbreport = KBaseReport(self.callback_url,token=self.token)
        self.struct_utils = ProteinStructureUtils(self.callback_url,token=self.token)
        self.anno_api = cb_annotation_ontology_api(self.callback_url,token=self.token)
        self.dfu = DataFileUtil(self.callback_url,token=self.token)
        self.api = KBAnnotationModule("KBAnnotationApps",self.wsclient,self.anno_api,self.struct_utils,config['scratch'],self.config)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass



    def PDBAnnotation(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "PDBAnnotationParams" -> structure:
           parameter "genomes" of list of type "Genome_ref" (Reference to a
           Genome object in the workspace @id ws KBaseGenomes.Genome
           KBaseGenomeAnnotations.GenomeAnnotation), parameter "suffix" of
           String
        :returns: instance of type "PDBAnnotationResults" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN PDBAnnotation
        api_output = self.api.PDBAnnotation(params)            
        output = self.save_report_to_kbase(api_output,)
        self.api.transfer_outputs(output,api_output,["data"])
        #END PDBAnnotation

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method PDBAnnotation return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
