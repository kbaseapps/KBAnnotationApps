# -*- coding: utf-8 -*-
import os
import time
import unittest
import json
from configparser import ConfigParser

from KBAnnotationApps.KBAnnotationAppsImpl import KBAnnotationApps
from KBAnnotationApps.KBAnnotationAppsServer import MethodContext
from KBAnnotationApps.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class KBAnnotationAppsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('KBAnnotationApps'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'KBAnnotationApps',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = KBAnnotationApps(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        api = self.serviceImpl.api
        api.initialize_call("PDBAnnotation",{"genome_refs":["68716/ecoli-5"],"workspace":68716},True)
        
        #Testing retrieval of genome proteins
        output = api.genome_to_proteins("68716/ecoli-5")
        self.assertEqual(len(output),5)
        self.assertEqual(output[0][0],"b0001")
        
        #Testing query function
        output = api.query_rcsb_with_proteins(output)
            
        #Testing saving of annotation genome
        output = api.add_annotations_to_genome("68716/ecoli-5",".pdb",output)
        #self.assertEqual(output["ref"],"68716/ecoli-5.pdb")
        
        output = api.PDBAnnotation({
            "workspace":68716,
            "genome_ref":"68716/ecoli-5",
            "suffix":".pdb",
            "similarity_threshold_type":"evalue",
            "similarity_threshold":0.00001,
            "return_data":False
        })