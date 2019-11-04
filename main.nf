#!/usr/bin/env nextflow 

// 10x --> cds --> garnett stuff 


// parse 10x directory into CDS object 
input_mtx = params.matrix
input_bracodes = params.barcodes
input_genes = params.genes
INPUT_MTX = channel.fromPath(input_mtx)
INPUT_BARCODES = channel.fromPath(input_bracodes)
INPUT_GENES = channel.fromPath(input_genes)

process build_CDS{
    conda envs/moncle3.yaml

    //TODO: add resourse handling

    input: 
        file(10x_matrix) from INPUT_10X
        file(10x_barcodes) from INPUT_BARCODES
        file(10x_genes) from INPUT_GENES

    output:
        file(garnett_cds.rds) into GARNETT_CDS
    
    """
    monocle3 create\
            --expression-matrix ${10x_matrix}\
            --cell-metadata ${10x_barcodes}\
            --gene-annotation ${10x_genes}\
            garnett_cds.rds

    """
}


// check supplied markers 
marker_genes = params.marker_genes
MARKER_GENES = channel.fromPath(marker_genes)
//value channels for re-use
GARNETT_CDS = GARNETT_CDS.first()
MARKER_GENES = MARKER_GENES.first()

process check_markers{ 
    publishDir "data/output_dir", mode: 'copy'
    
    conda envs/garnett-cli.yaml 

    input:
        file(marker_genes) from MARKER_GENES
        file(garnett_cds) from GARNETT_CDS

    output:
        file(marker_genes_checked.txt) into CHECKED_MARKER_GENES
        file(marker_plot.png) into MARKER_PLOT

    """
    garnett_check_markers 
            --cds-object ${garnett_cds}\
            --marker-file-path ${marker_genes}
            --database ${params.database}\
            --cds-gene-id-type ${params.cds_gene_id_type}\
            --marker-file-gene-id-type ${marker_gene_id_type}\
            --marker-output-path marker_genes_checked.txt\
            --plot-output-path marker_plot.png

    """
}


//TODO: at this point, the pipeline should pause and let the user inspect the markers (?)

process train_classifier{
    conda envs/garnett-cli.yaml 

    input:
        file(garnett_cds) from GARNETT_CDS
        file(marker_genes) from MARKER_GENES

    output:
        file(trained_classifier.rds) into TRAINED_CLASSIFIER 

    """
    garnett_train_classifier.R\
            --cds-object ${garnett_cds}\
            --marker-file-path ${marker_genes}\
            --database params.database\
            --output-path trained_classifier.rds 

    """

}

//get feature genes 
TRAINED_CLASSIFIER = TRAINED_CLASSIFIER.first()
process get_feature_genes{

    publishDir "data/output_dir", mode: 'copy'
    conda envs/garnett-cli.yaml

    input:
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file(feature_genes.txt) into FEATURE_GENES 

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
    conda envs/garnett-cli.yaml

    input:
        file(garnett_cds) from GARNETT_CDS
        file(classifier) from TRAINED_CLASSIFIER

    output:
        file(garnett_cds) into CLASSIFIED_CELL_TYPES


    """
    garnett_classify_cells.R\
            --cds-object ${garnett_cds}
            --classifier-object ${classifier}
            --database ${params.database}

    """
}

















