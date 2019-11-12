#!/usr/bin/env nextflow 

// parse 10x directory into CDS object 
input_mtx = params.matrix
input_bracodes = params.barcodes
input_genes = params.genes
INPUT_MTX = Channel.fromPath(input_mtx)
INPUT_BARCODES = Channel.fromPath(input_bracodes)
INPUT_GENES = Channel.fromPath(input_genes)

process build_CDS{
    conda "envs/garnett-cli.yaml"

    //TODO: add resourse handling

    input: 
        file(matrix) from INPUT_MTX
        file(barcodes) from INPUT_BARCODES
        file(genes) from INPUT_GENES

    output:
        file("garnett_cds.rds") into GARNETT_CDS
    
    """
    parse_expr_data.R\
            --expression-matrix ${matrix}\
            --phenotype-data ${barcodes}\
            --feature-data ${genes}\
            --output-file garnett_cds.rds
    """
}

// check supplied markers 
marker_genes = params.marker_genes
MARKER_GENES = Channel.fromPath(marker_genes)
//value channels for re-use
GARNETT_CDS = GARNETT_CDS.first()
MARKER_GENES = MARKER_GENES.first()

process check_markers{ 
    publishDir "data/output_dir", mode: 'copy'
    
    conda "envs/garnett-cli.yaml" 

    input:
        file(marker_genes) from MARKER_GENES
        file(garnett_cds) from GARNETT_CDS

    output:
        file("marker_genes_checked.txt") into CHECKED_MARKER_GENES
        file("marker_plot.png") into MARKER_PLOT

    """
    garnett_check_markers.R\
            --cds-object ${garnett_cds}\
            --marker-file-path ${marker_genes}\
            -d ${params.database}\
            --cds-gene-id-type ${params.cds_gene_id_type}\
            --marker-output-path marker_genes_checked.txt\
            --plot-output-path marker_plot.png

    """
}


//TODO: at this point, the pipeline should pause and let the user inspect the markers (?)

process train_classifier{
    conda "envs/garnett-cli.yaml" 

    input:
        file(garnett_cds) from GARNETT_CDS
        file(marker_genes) from MARKER_GENES

    output:
        file("trained_classifier.rds") into TRAINED_CLASSIFIER 

    """
    garnett_train_classifier.R\
            --cds-object ${garnett_cds}\
            --marker-file-path ${marker_genes}\
            -d ${params.database}\
            --cds-gene-id-type ${params.cds_gene_id_type}\
            --marker-file-gene-id-type ${params.marker_gene_id_type}\
            --classifier-gene-id-type ${params.classifier_gene_type}\
            -n ${params.n_outgroups}\
            --output-path trained_classifier.rds 

    """

}

//get feature genes 
//TRAINED_CLASSIFIER = TRAINED_CLASSIFIER.first()
process get_feature_genes{

    publishDir "data/output_dir", mode: 'copy'
    conda "envs/garnett-cli.yaml"

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
    publishDir "data/output_dir", mode: 'copy'
    conda "envs/garnett-cli.yaml"

    input:
        file(garnett_cds) from GARNETT_CDS
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file(garnett_cds) into CLASSIFIED_CELL_TYPES


    """
    garnett_classify_cells.R\
            --cds-object ${garnett_cds}\
            --cds-gene-id-type ${params.cds_gene_id_type}\
            --classifier-object ${classifier}\
            --database ${params.database}

    """
}

// get cell metadata for further processing 
process get_metadata{
    publishDir "data/output_dir", mode: 'copy'
    conda "envs/garnett-cli.yaml"

    input:
        file(garnett_cds) from CLASSIFIED_CELL_TYPES

    output:
        file("garnett_metadata.tsv") into CDS_METADATA

    """
    get_metadata.R ${garnett_cds} garnett_metadata.tsv 
    """

} 

//obtain output in standard format 
process get_output{
    publishDir "data/output_dir", mode: 'copy'
    conda "envs/standardised_output.yaml"

    input:
        file(garnett_meta) from CDS_METADATA

    output:
        file("garnett_cell_type_annotations.tsv") into ANNOTATION_OUTPUT

    // make a conda package for our own channel 
    """
    get_output.R\
        --input-file ${garnett_meta}\
        --cell-id-field ${params.cell_id_field}\
        --predicted-cell-type-field ${params.predicted_cell_type_field}\
        --output-file-path garnett_cell_type_annotations.tsv 
    """
}