/*
A KBase module: KBAnnotationApps
*/

module KBAnnotationApps {
    /*
    Reference to a Genome object in the workspace
    @id ws KBaseGenomes.Genome KBaseGenomeAnnotations.GenomeAnnotation
    */
    typedef string Genome_ref;
    
    typedef structure {
        list<Genome_ref> genome_refs;
        string workspace;
        string similarity_threshold_type;
        float similarity_threshold;
        string suffix;
    } PDBAnnotationParams;
    
    typedef structure {
        string report_name;
        string report_ref;
    } PDBAnnotationResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef PDBAnnotation(PDBAnnotationParams params) returns (PDBAnnotationResults output) authentication required;

};
