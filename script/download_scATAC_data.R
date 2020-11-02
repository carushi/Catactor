# Downloads

library(GEOquery)

gse <- getGEO("GSE111586", GSEMatrix=TRUE)
filePaths = getGEOSuppFiles("GSE111586")

# SNARE-seq
dir <- "GSE126074"
peaks <- "GSE126074_P0_BrainCortex_SNAREseq_chromatin.peaks.tsv"
counts <- "GSE126074_P0_BrainCortex_SNAREseq_chromatin.counts.mtx"
barcodes <- "GSE126074_P0_BrainCortex_SNAREseq_chromatin.barcodes.tsv"


dir <- "GSE100033"
peaks <- "GSM266*.ygi.txt"
counts <- "GSM266*.cell.mat" # non-sparse
barcodes <- "GSM266*.xgi.txt"
