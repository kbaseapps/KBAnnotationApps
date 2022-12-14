# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except ImportError:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport


class ProteinStructureUtils(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='dev',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def batch_import_pdbs_from_metafile(self, params, context=None):
        """
        batch_import_pdbs_from_metafile: import a batch of ProteinStructures from PDB files
        :param params: instance of type "BatchPDBImportParams" (Input/Output
           of the batch_import_pdbs_from_metafile structures_name:
           Proteinstructures object name workspace_name: workspace name for
           object to be saved to metadata_staging_file_path: path to a
           spreadsheet file that lists the metadata of PDB files and their
           KBase metadata) -> structure: parameter
           "metadata_staging_file_path" of String, parameter
           "structures_name" of String, parameter "workspace_name" of type
           "workspace_name" (workspace name of the object)
        :returns: instance of type "BatchPDBImportOutput" -> structure:
           parameter "structures_ref" of String, parameter "report_name" of
           String, parameter "report_ref" of String
        """
        return self._client.run_job('ProteinStructureUtils.batch_import_pdbs_from_metafile',
                                    [params], self._service_ver, context)

    def import_rcsb_structures(self, params, context=None):
        """
        :param params: instance of type "ImportRCSBParams" (Input/output of
           the import_rcsb_structures function rcsb_infos: a list of
           RCSBInfoStruct's structures_name: Proteinstructures object name
           workspace_name: workspace name for object to be saved to) ->
           structure: parameter "rcsb_infos" of list of type "RCSBInfoStruct"
           (The information required by the importing app rcsb_id: rcsb
           structure id extension: file extension for the structure ('pdb' or
           'cif') narrative_id: a KBase narrative id genome_name: a KBase
           genome name in the respective narrative of narrative_id
           feature_id: a KBase feature id in the respective narrative of
           narrative_id is_model: a value of 0 or 1 to indicate the structure
           is exprimental or computational) -> structure: parameter "rcsb_id"
           of String, parameter "extension" of String, parameter
           "narrative_id" of String, parameter "genome_name" of String,
           parameter "feature_id" of String, parameter "is_model" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "structures_name" of String, parameter "workspace_name"
           of type "workspace_name" (workspace name of the object)
        :returns: instance of type "ImportRCSBStructOutput" -> structure:
           parameter "structures_ref" of String, parameter "report_name" of
           String, parameter "report_ref" of String
        """
        return self._client.run_job('ProteinStructureUtils.import_rcsb_structures',
                                    [params], self._service_ver, context)

    def export_pdb_structures(self, params, context=None):
        """
        :param params: instance of type "ExportParams" (Input/output of the
           export_pdb_structures function input_ref: generics object
           reference) -> structure: parameter "input_ref" of type "obj_ref"
           (An X/Y/Z style reference @id ws)
        :returns: instance of type "ExportStructOutput" -> structure:
           parameter "shock_ids" of list of String
        """
        return self._client.run_job('ProteinStructureUtils.export_pdb_structures',
                                    [params], self._service_ver, context)

    def query_rcsb_structures(self, params, context=None):
        """
        :param params: instance of type "QueryRCSBStructsParams"
           (Input/output of the query_rcsb_structures function
           sequence_strings: a list of protein sequences uniprot_ids: a list
           of uniprot ids ec_numbers: a list of ec numbers inchis: a list of
           InChI strings smiles: a list of SMILES strings evalue_cutoff:
           threshold of homology search identity_cutoff: threshold for
           sequence identity match workspace_name: workspace name for objects
           to be saved to @optional sequence_strings uniprot_ids ec_numbers
           inchis smiles evalue_cutoff identity_cutoff) -> structure:
           parameter "sequence_strings" of list of String, parameter
           "uniprot_ids" of list of String, parameter "ec_numbers" of list of
           String, parameter "inchis" of list of String, parameter "smiles"
           of list of String, parameter "evalue_cutoff" of Double, parameter
           "identity_cutoff" of Double, parameter "logical_and" of type
           "boolean" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "workspace_name" of type "workspace_name" (workspace
           name of the object)
        :returns: instance of type "QueryRCSBStructsOutput" -> structure:
           parameter "rcsb_ids" of list of String, parameter "rcsb_scores" of
           unspecified object, parameter "report_name" of String, parameter
           "report_ref" of String
        """
        return self._client.run_job('ProteinStructureUtils.query_rcsb_structures',
                                    [params], self._service_ver, context)

    def query_rcsb_annotations(self, params, context=None):
        """
        :param params: instance of type "QueryRCSBAnnotationsParams"
           (Input/output of the query_rcsb_annotations function
           sequence_strings: a list of protein sequences evalue_cutoff:
           threshold of homology search identity_cutoff: threshold for
           sequence identity match workspace_name: workspace name for objects
           to be saved to @optional evalue_cutoff identity_cutoff) ->
           structure: parameter "sequence_strings" of list of String,
           parameter "evalue_cutoff" of Double, parameter "identity_cutoff"
           of Double, parameter "workspace_name" of type "workspace_name"
           (workspace name of the object)
        :returns: instance of type "QueryRCSBAnnotationsOutput" -> structure:
           parameter "rcsb_hits" of unspecified object
        """
        return self._client.run_job('ProteinStructureUtils.query_rcsb_annotations',
                                    [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.run_job('ProteinStructureUtils.status',
                                    [], self._service_ver, context)
