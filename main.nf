#!/usr/bin/env nextflow 

// build reference CDS object 
REF_10X_DIR = Channel.fromPath(params.ref_10x_dir)
process build_ref_CDS_object{
    conda "${baseDir}/envs/monocle3-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input: 
        file(ref_10x_dir) from REF_10X_DIR

    output:
        file("ref_cds.rds") into REF_CDS
    
    """
    monocle3 create ref_cds.rds\
                --expression-matrix ${ref_10x_dir}/matrix.mtx\
                --cell-metadata ${ref_10x_dir}/barcodes.tsv\
                --gene-annotation ${ref_10x_dir}/genes.tsv\
    """
}
REF_CDS = REF_CDS.first()

// build query CDS object
QUERY_10X_DIR = Channel.fromPath(params.query_10x_dir)
process build_query_CDS_object{
    conda "${baseDir}/envs/monocle3-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input: 
        file(query_10x_dir) from QUERY_10X_DIR

    output:
        file("query_cds.rds") into QUERY_CDS
    
    """
    monocle3 create query_cds.rds\
                --expression-matrix ${query_10x_dir}/matrix.mtx\
                --cell-metadata ${query_10x_dir}/barcodes.tsv\
                --gene-annotation ${query_10x_dir}/genes.tsv\
    """
}
QUERY_CDS = QUERY_CDS.first()

// transform markers from SCXA format into Garnett
SCXA_MARKER_GENES = Channel.fromPath(params.marker_genes).first()
process transform_markers{
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(scxa_markers) from SCXA_MARKER_GENES

    output:
        file("garnett_markers.txt") into GARNETT_MARKERS
        file("markers_list.rds") into MARKERS_LIST

    """
    transform_marker_file.R\
            --input-marker-file ${scxa_markers}\
            --pval-col ${params.pval_col}\
            --gene-names ${params.gene_names}\
            --groups-col ${params.groups_col}\
            --marker-list markers_list.rds\
            --garnett-marker-file garnett_markers.txt
    """
}


// check supplied markers 
process check_markers{ 
    publishDir "data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    input:
        file(marker_genes) from GARNETT_MARKERS
        file(ref_cds) from REF_CDS

    output:
        file("marker_genes_checked.txt") into CHECKED_MARKER_GENES
        

    """
    garnett_check_markers.R\
            --cds-object ${ref_cds}\
            --marker-file-path ${marker_genes}\
            --marker-file-gene-id-type ${params.marker_gene_id_type}\
            -d ${params.database}\
            --cds-gene-id-type ${params.ref_cds_gene_id_type}\
            --marker-output-path marker_genes_checked.txt\
    """
}


process update_markers {

    conda "${baseDir}/envs/garnett-cli.yaml" 
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
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

process train_classifier{
    conda "${baseDir}/envs/garnett-cli.yaml" 

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
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

// classify cells 
process classify_cells{
    publishDir "${baseDir}/data/output_dir", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
    maxRetries 10
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

//obtain output in standard format 
process get_output{
    publishDir "${params.results_dir}", mode: 'copy'
    conda "${baseDir}/envs/garnett-cli.yaml"

    input:
        file(classified_cells) from CLASSIFIED_CELL_TYPES
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file("garnett_output.txt") into ANNOTATION_OUTPUT

    """
    garnett_get_std_output.R\
            --input-object ${classified_cells}\
            --predicted-cell-type-field ${params.predicted_cell_type_field}\
	    --classifier ${classifier}\
            --output-file-path garnett_output.txt\
    """
}
