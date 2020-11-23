convert.en2eg <- function() {
    require(org.Mm.eg.db)
    require(biomaRt)
    xx <- mappedkeys(org.Mm.egENSEMBLTRANS2EG)
    print(length(xx))
    print(xx[1:5])
    mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast")
    mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
    print(listAttributes(mart))
    genes <- getBM(
      filters="ensembl_transcript_id",
      attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "entrezgene_id", "refseq_mrna"),
      values=xx, 
      mart=mart)
    colnames(genes) <- c("Ensembl", "En_gene", "Name", "Entrez", "Refseq")
    write.table(genes, file="mm.en.et.csv", sep=",", quote=F)
}

write.tfbs.region <- function(jasp_dir) {
    library(TFBSTools)
    suppressMessages(library(JASPAR2018))
    require(BSgenome.Mmusculus.UCSC.mm10)
    opts <- list()
    opts[["species"]] <- 10090
    opts[["type"]] <- "SELEX"
    opts[["all_versions"]] <- TRUE
    PFMatrixList <- getMatrixSet(JASPAR2018, opts)
    for (n in names(PFMatrixList)) {
        tf_name <- name(PFMatrixList[[n]])
        pwm <- toPWM(PFMatrixList[[n]])
        print(c(n, tf_name))
        next
        for (chr in seqnames(BSgenome.Mmusculus.UCSC.mm10)) {
            if (nchar(chr) > 5) next
            siteset <- searchSeq(pwm, BSgenome.Mmusculus.UCSC.mm10[[chr]], min.score="80%", mc.cores=10)
            write.table(writeGFF3(siteset), file=file.path(jasp_dir, paste0("MM_jaspar_", n, "_", chr, ".tsv")), sep="\t")
        }
    }
}
get.delim <- function(annot_file) {
    sep <- '\t'
    if (grepl(".csv", annot_file)) sep <- ','
    return(sep)
}

annotate.tfbs.region <- function(annot_file, jasp_dir) {
    #if (!file.exists(file.path(jasp_dir, "MM_jaspar_MA0004.1_chr1.tsv")))
        write.tfbs.region(jasp_dir)
    tfs <- unique(unlist(lapply(list.files(jasp_dir), function(x){return(unlist(strsplit(x, "_", fixed=TRUE))[3])})))
    bin_ann <- read.table(annot_file, header=T, sep=get.delim(annot_file))
    chrs <- unique(bin_ann[,"chr"])
    for (tf in tfs) {
        dict = list()
    #     for (chr in chrs) dict[chr] = 
    # subject <- promoters(txdb, upstream=dist, downstream=dist)
    # ov <- findOverlaps(query, subject, minoverlap=1L, type="any", select="all", ignore.strand=TRUE)

    #     bin_ann <- cbind(bin_ann, )
    #     names(bin_ann)[dim(bin_ann)[2]] <- tf
    }
    file_base = unlist(strsplit(annot_file, ".", fixed=TRUE))[1]
    annot_out_file <- paste0(dirname(annot_file), "/", file_base, "_tfbs.csv")
    write.table(bin_ann, annot_out_file, sep=get.delim(annot_out_file))
}

annotate.genic.region.single <- function(annot_file) {
    options(scipen=99)
    require(annotatr)
    require(GenomicRanges)
    require(GenomicFeatures)
    require(tools)
    bin_ann <- read.table(annot_file, header=T, sep=get.delim(annot_file), stringsAsFactor=F)
    if (any(colnames(bin_ann) == "chr")) {
        bin_chr_ann <- bin_ann
        bin_chr_ann['chr'] <- sapply(bin_chr_ann['chr'], function(x){return(paste0("chr", gsub("^chr", "", x)))})
        if (colnames(bin_chr_ann)[1] == "X")
            bin_chr_ann <- bin_chr_ann[,-c(1)]
        print(head(bin_chr_ann))
        colnames(bin_chr_ann)[colnames(bin_chr_ann) == "strand"] = "gene_strand"
        regions <- makeGRangesFromDataFrame(data.frame(apply(bin_chr_ann, c(1,2), as.character)), keep.extra.columns = TRUE, ignore.strand=TRUE)
        annots = c('mm10_basicgenes', 'mm10_genes_intergenic',
            'mm10_genes_intronexonboundaries')
        annotations = build_annotations(genome = 'mm10', annotations = annots)
        dm_annotated = annotate_regions(
            regions = regions,
            annotations = annotations,
            ignore.strand = TRUE,
            quiet = FALSE)
        annots_order = c(
            'mm10_genes_promoters',
            'mm10_genes_5UTRs',
            'mm10_genes_exons',
            'mm10_genes_intronexonboundaries',
            'mm10_genes_introns',
            'mm10_genes_3UTRs',
            'mm10_genes_1to5kb',
            'mm10_genes_intergenic')
        df_dm_annotated = data.frame(dm_annotated)
        stopifnot(all(sapply(unique(df_dm_annotated[,"annot.type"]), function(x){return(any(x == annots_order))})))
        df_dm_annotated <- cbind(df_dm_annotated, annot.index=sapply(df_dm_annotated[,"annot.type"], function(x){return(which(x == annots_order))}))
        df_dm_annotated <- df_dm_annotated[order(df_dm_annotated[,"annot.index"]),]
        df <- df_dm_annotated[!duplicated(df_dm_annotated[,c("seqnames", "start", "end")]),c("seqnames","start", "end","annot.gene_id", "annot.symbol", "annot.type", "annot.index")]
        df <- cbind(df, 'annot.proximal'=unlist(sapply(df[,"annot.index"], function(x){if(x == 1)return(1); return(0)})))
        df <- cbind(df, 'annot.distal'=unlist(sapply(df[,"annot.index"], function(x){if(any(x == c(2, 3, 4, 5, 6, 7)))return(1); return(0)})))
        stopifnot(dim(df)[1] <= dim(bin_ann)[1])
        df[,1] <- unlist(sapply(df[,1], function(x){return(substr(x, 4, nchar(as.character(x))))}))
        colnames(df)[1] <- "chr"
        m <- merge(bin_ann, df, by=c("chr", "start", "end"), all.x=TRUE)
        print(head(m))
        m <- m[!duplicated(m[,1:3]),]
        m <- m[order(m[,1], m[,2], m[,3]),]
        vec <- rep(NA, dim(m)[1]);
        vec[m[,"annot.proximal"] > 0 & !is.na(m[,"annot.proximal"])] <- 1:length(which(m[,"annot.proximal"] > 0 & !is.na(m[,"annot.proximal"])))
        m <- cbind(m, 'annot.proximal.uniq'=vec)
        vec <- rep(NA, dim(m)[1]);
        vec[m[,"annot.distal"] > 0 & !is.na(m[,"annot.distal"])] <- 1:length(which(m[,"annot.distal"] > 0 & !is.na(m[,"annot.distal"])))
        m <- cbind(m, 'annot.distal.uniq'=vec)
    } else {
        m <- bin_ann
    }
    colnames(m)[colnames(m) == "gene_strand"] = "strand"
    if (any("X" == colnames(m))) m <- m[,-which("X" == colnames(m))]
    if (any("name" == colnames(m))) m <- m[,-which("name" == colnames(m))]
    if (any("score" == colnames(m))) m <- m[,-which("score" == colnames(m))]
    file_base = annot_file
    annot_out_file <- paste0(dirname(annot_file), "/", file_path_sans_ext(basename(file_base)), "_annot.csv")
    print(c(annot_out_file, basename(file_base), file_path_sans_ext(basename(file_base))))
    write.table(m, annot_out_file, sep=get.delim(annot_out_file), quote=F, row.names=F)
}

assign.gene.names <- function(query=c(), dist=50000) {
    if (length(query) == 0) {
        df <- data.frame(chr=c("1", "2", "3"), start=c(1000, 10000, 100000), end=c(10000000, 11000, 101000))
        query <- makeGRangesFromDataFrame(df)
    }
    txdb <- makeTxDbFromEnsembl(organism="Mus musculus",
                    release=NA,
                    server="ensembldb.ensembl.org",
                    username="anonymous", password=NULL, port=0L,
                    tx_attrib=NULL)
    subject <- promoters(txdb, upstream=dist, downstream=dist)
    ov <- findOverlaps(query, subject, minoverlap=1L, type="any", select="all", ignore.strand=TRUE)
    #print(ov)
    ov <- data.frame(ov)
    print(paste(names(subject)[ov[which(ov[,1] == 1),2]], collapse=","))
    return(unlist(sapply(1:length(query), function(x){return(paste(names(subject)[ov[which(ov[,1] == x),2]], collapse=","))})))
}


args <- commandArgs(trailingOnly=TRUE)
print(args)
stopifnot(length(args) >= 2)
fname <- args[2]
if (args[1] == "tfbs") {
    annotate.tfbs.region(fname, args[3])
} else if (args[1] == "promoter") {
    for (fname in args[2:length(args)]) {
        annotate.genic.region.single(fname)
    }
}
