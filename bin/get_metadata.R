#!/usr/bin/env Rscript 
suppressPackageStartupMessages(require(garnett))
cds_path = commandArgs(TRUE)[1]
metadata_path = commandArgs(TRUE)[2]

cds = readRDS(cds_path)
cell_meta = pData(cds)
write.table(cell_meta, file=metadata_path, sep="\t")
