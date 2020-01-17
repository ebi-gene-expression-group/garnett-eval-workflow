#!/usr/bin/env nextflow 

// build query and reference CDS object 
REF_10X_DIR = Channel.fromPath(params.ref_10x_dir)
QUERY_10X_DIR = Channel.fromPath(params.query_10x_dir)
process build_CDS_objects{
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input: 
        file(ref_10x_dir) from REF_10X_DIR
        file(query_10x_dir) from QUERY_10X_DIR

    output:
        file("ref_cds.rds") into REF_CDS
        file("query_cds.rds") into QUERY_CDS
    
    """
    parse_expr_data.R\
            --ref-10x-dir ${ref_10x_dir}\
            --query-10x-dir ${query_10x_dir}\
            --ref-output-cds ref_cds.rds\
            --query-output-cds query_cds.rds
    """
}

// transform markers from SCXA format into Garnett
SCXA_MARKER_GENES = Channel.fromPath(params.marker_genes)
process transform_markers{
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(scxa_markers) from SCXA_MARKER_GENES

    output:
        file("garnett_markers.txt") into GARNETT_MARKERS
        file("markers_list.rds") into MARKERS_LIST

    """
    transform_marker_file.R\
            --input-marker-file ${cxa_markers}\
            --marker-list markers_list.rds\
            --garnett-marker-file garnett_markers.txt
    """
}


// check supplied markers 
//value channels for re-use
REF_CDS = REF_CDS.first()
process check_markers{ 
    publishDir "data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(marker_genes) from GARNETT_MARKERS
        file(ref_cds) from REF_CDS

    output:
        file("marker_genes_checked.txt") into CHECKED_MARKER_GENES
        file("marker_plot.png") into MARKER_PLOT

    """
    garnett_check_markers.R\
            --cds-object ${ref_cds}\
            --marker-file-path ${marker_genes}\
            -d ${params.database}\
            --cds-gene-id-type ${params.ref_cds_gene_id_type}\
            --marker-output-path marker_genes_checked.txt\
            --plot-output-path marker_plot.png

    """
}

process update_markers {

    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(marker_list) from MARKERS_LIST
        file(marker_summary) from CHECKED_MARKER_GENES

    output: 
        file("garnett_markers_upd.txt") into GARNETT_MARKERS_UPD

    """
    update_marker_file.R\
            --marker-list-obj ${marker_list}\
            --marker-check-file ${marker_summary}\
            --updated-marker-file garnett_markers_upd.txt
    """

}
GARNETT_MARKERS_UPD = GARNETT_MARKERS_UPD.first()

process train_classifier{
    conda "${baseDir}/envs/garnett-cli.yaml" 

    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(ref_cds) from REF_CDS
        file(marker_genes) from GARNETT_MARKERS_UPD

    output:
        file("trained_classifier.rds") into TRAINED_CLASSIFIER 

    """
    garnett_train_classifier.R\
            --cds-object ${ref_cds}\
            --marker-file-path ${marker_genes}\
            -d ${params.database}\
            --cds-gene-id-type ${params.ref_cds_gene_id_type}\
            --marker-file-gene-id-type ${params.marker_gene_id_type}\
            --classifier-gene-id-type ${params.classifier_gene_type}\
            -n ${params.n_outgroups}\
            --output-path trained_classifier.rds 
    """
}

//get feature genes 
process get_feature_genes{

    publishDir "${baseDir}/data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml"

    // resource handling 
    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file("feature_genes.txt") into FEATURE_GENES 

    """
    garnett_get_feature_genes.R\
            --classifier-object ${classifier}\
            --database ${params.database}\
            --output-path feature_genes.txt
    """

}

// classify cells 
process classify_cells{
    publishDir "${baseDir}/data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.attempt < 4  ? 'retry' : 'ignore' }   
    maxRetries 4
    memory { 16.GB * task.attempt }

    input:
        file(query_cds) from QUERY_CDS
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file("cds_classified.rds") into CLASSIFIED_CELL_TYPES


    """
    garnett_classify_cells.R\
            --cds-object ${query_cds}\
            --cds-gene-id-type ${params.query_cds_gene_id_type}\
            --classifier-object ${classifier}\
            --database ${params.database}\
            --cds-output-obj cds_classified.rds
    """
}

// get cell metadata for further processing 
process get_metadata{
    publishDir "${baseDir}/data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml"

    input:
        file(cds_classified) from CLASSIFIED_CELL_TYPES

    output:
        file("garnett_metadata.tsv") into CDS_METADATA

    """
    get_metadata.R ${garnett_cds} garnett_metadata.tsv 
    """
} 

//obtain output in standard format 
process get_output{
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/standardised_output.yaml"

    input:
        file(garnett_meta) from CDS_METADATA

    output:
        file("garnett_output.txt") into ANNOTATION_OUTPUT

    """
    get_output.R\
        --input-file ${garnett_meta}\
        --cell-id-field ${params.cell_id_field}\
        --predicted-cell-type-field ${params.predicted_cell_type_field}\
        --output-file-path garnett_output.txt
    """
}
