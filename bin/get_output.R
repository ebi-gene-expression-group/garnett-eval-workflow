#!/usr/bin/env Rscript 

# Produce standardised output from a variety of tools 

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))
# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-i", "--input-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the input metadata file in .tsv format"
    ),

    make_option(
        c("-c", "--cell-id-field"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Column name of the cell id annotations. If not supplied, it is assumed
                that cell ids are represented by index'
    ),

    make_option(
        c("-p", "--predicted-cell-type-field"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Column name of the predicted cell type annotation'
    ), 

    make_option(
        c("-o", "--output-file-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the produced output file in .tsv format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c('input_file',
                                                'predicted_cell_type_field',
                                                'output_file_path'))


# parse the input table
table = read.table(opt$input_file, sep="\t")

# can't use NAs in workflows? 
if(opt$cell_id_field != "null"){
    cell_id = table[, opt$cell_id_field]
} else{
    cat("no index column provided; use row names\n")
    cell_id = row.names(table)
}

predicted_label = as.character(table[, opt$predicted_cell_type_field])
output_table = data.frame(cbind(cell_id, predicted_label))
write.csv(output_table, file = opt$output_file_path, sep="\t")
