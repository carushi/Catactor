library(parallel)
library(ggplot2)
library(viridis)
library(GenomicFeatures)
library(dplyr)
library(multidplyr)

annots_order <<- c(
'mm10_genes_promoters',
'mm10_genes_5UTRs',
'mm10_genes_exons',
'mm10_genes_intronexonboundaries',
'mm10_genes_introns',
'mm10_genes_3UTRs',
'mm10_genes_1to5kb',
'mm10_genes_intergenic')
CORES <<- 10

analyze.GO.result <- function(go_fname, header) {
    result <- read.table(go_fname, sep="\t", header=T, stringsAsFactor=F)
    annot <- read.table("mm.go.description.csv", sep=",", header=T, stringsAsFactor=F)
    annot = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
    rownames(annot) <- gsub(":", ".", annot[,"GOID"])
    print(dim(result))
    m <- merge(result, annot, by.x=0, by.y=0, all.x=TRUE)
    print(dim(m))
    m[,"max_pvalue"] <- unlist(sapply(1:dim(m)[1], function(x){return(max(m[x, "pgreat.x"], m[x, "pgreat.x"]))}))
    m <- m[order(m[,"max_pvalue"], decreasing=FALSE),]
    print(head(m))
    write.table(m, paste0("GO_enrichment_annotated_", header, ".tsv"), quote=T, sep="\t")
    
}

compute.enrichment <- function(row_or_column, chr, batch, bin, size, rel=-1, sup=-1) {
    if (sup > 0) {
        header <- paste(row_or_column, chr, batch, bin, "size", size, "sup", sup, sep="_")
    } else if (rel > 0) {
        header <- paste(row_or_column, chr, batch, bin, "size", size, rel, sep="_")
    } else {
        return(0)
    }
    fname <- paste0("data_frame_", header, "_entrez.tsv")
    print(fname)
    if (!file.exists(fname))
        map.to.go.genes(header)
    else
        print('ok')
    compute.enrichment.go(header)
    compute.enrichment.chr(header)
    compute.networks(header)
    compute.enrichment.biomarker(header)
}

convert.to.coexp <- function(mat) {
    coexp_mat <- matrix(0, nrow=rownames(mat), ncol=rownames(mat))
    rownames(coexp_mat) <- rownames(mat); colnames(coexp_mat) <- rownames(mat);
    for (i in 1:dim(mat)[2]) {
        genes <- combn(where(mat[i,]), 2)
        for (j in 1:dim(genes)[2]) {
            coexp_mat[genes[1,j],genes[2,j]] <- 1
            coexp_mat[genes[2,j],genes[2,j]] <- 1
        }
    }
    for (i in 1:dim(coexp_mat)[1]) {
        coexp_mat[i,i] <- 0
    }
    return(coexp_mat)
}

serach.symbol.mat <- function(mat, lnames, rnames) {
    return(factor(unlist(sapply(1:length(lnames), function(x) {
        lgenes <- unlist(strsplit(lnames[i], ',', fixed=TRUE))
        rgenes <- unlist(strsplit(rnames[i], ',', fixed=TRUE))
        lindex <- sapply(lgenes, function(x){return(which(rownames(mat) == x))})
        rindex <- sapply(rgenes, function(x){return(which(rownames(mat) == x))})
        if (length(lindex) == 0 && length(rindex) == 0) {
            return("None")
        } else if (length(lindex) > 0 || length(rindex) > 0) {
            return("Either")
        } else {
            if (any(mat[lindex, rindex] == 1))
                return("Coacc")
            else
                return("Both")
        }
    })), levels=c("None", "Either", "Both", "Coacc")))
}

serach.symbol.vec <- function(vec, lnames, rnames) {
    return(factor(unlist(sapply(1:length(lnames), function(x) {
        lgenes <- unlist(strsplit(lnames[i], ',', fixed=TRUE))
        rgenes <- unlist(strsplit(rnames[i], ',', fixed=TRUE))
        lindex <- sapply(lgenes, function(x){return(which(vec == x))})
        rindex <- sapply(rgenes, function(x){return(which(vec == x))})
        if (length(lindex) == 0 && length(rindex) == 0) {
            return("None")
        } else if (length(lindex) > 0 || length(rindex) > 0) {
            return("Either")
        } else return("Both")
    })), levels=c("None", "Either", "Both")))
}

compute.enrichment.biomarkers <- function(header, value="count") {
    require(ggplot2)
    require(cowplot)
    require(viridis)
    fname <- paste0("data_frame_", header, "_entrez.tsv")
    data <- read.table(fname, sep="\t", header=T, stringsAsFactor=F)
    for (edge in c("all", "non_neighbor")) {
        if (edge != 'all') {
            data <- subset(data, ddist %in% c("inter", "<=10", ">10"))
        }
        fpath <- "/data/rkawaguc/scATAC/data/190425/gene_annotation_from_scRNA/"
        for (biomarker in c("tasic2016_gli_mat", "tasic2016_glu_mat", "tasic2016_gaba_mat", "tasic2016_mat", "markers_uniq")) {
            mat <- convert.to.coexp(read.table(paste0(fpath, biomarker, ".txt"), header=T, row.names=1))
            vec <- search.symbol.mat(mat, data[,"lnames"], data[,"rnames"])
            df <- data.frame(x=1:length(vec), y=data[,"count"], color=vec)
            g <- ggplot(data=df, aes(x=x, y=y, color=color))+geom_point()+scale_color_viridis(discrete=TRUE)+theme_cowplot()
            png(paste0("biomarker_overlap_", header, "_", edge, ".png"))
            plot(g)
            dev.off()
        }
        for (biomarker in c("random", "hvg")) {
            mat <- read.table(paste0(fpath, biomarker, ".txt"), header=F)
            vec <- search.symbol.vec(mat, data[,"lnames"], data[,"rnames"])
            df <- data.frame(x=1:length(vec), y=data[,"count"], color=vec)
            g <- ggplot(data=df, aes(x=x, y=y, color=color))+geom_point()+scale_color_viridis(discrete=TRUE)+theme_cowplot()
            png(paste0("biomarker_overlap_", header, "_", edge, ".png"))
            plot(g)
            dev.off()
        }
    }
}

compute.enrichment.chr <- function(header, value="count") {
    require(GenomicRanges)
    require(ggbio)
    #require(UpSetR)
    require(corrplot)
    chrom <- read.table("/data/rkawaguc/scATAC/data/190402/Annotation/mm10.chrom.sizes", sep="\t", header=F, stringsAsFactor=F)
    chrom <- subset(chrom, apply(chrom, c(1), function(x){return(nchar(x[1]) <= 5 && !any(x[1] == c("chrY", "chrM")))}))
    seqlength_mm10 <- chrom[,2]; names(seqlength_mm10) <- chrom[,1]
    chrom <- data.frame(chr=chrom[,1], start=rep(1, dim(chrom)[1]), end=chrom[,2])
    chrom <- makeGRangesFromDataFrame(chrom, keep.extra.columns=TRUE, ignore.strand=TRUE)
    seqlengths(chrom) <- seqlength_mm10
    fname <- paste0("data_frame_", header, "_entrez.tsv")
    data <- read.table(fname, sep="\t", header=T, stringsAsFactor=F)
    data[, "chr"] <- paste0("chr", data[, "chr"])
    data[, "chr.1"] <- paste0("chr", data[, "chr.1"])
    data[, "chr"] <- factor(data[, "chr"], levels=paste0("chr", c(1:19, "X", "Y", "M")))
    data[, "chr.1"] <- factor(data[, "chr.1"], levels=paste0("chr", c(1:19, "X", "Y", "M")))
    data <- data[1:20,]
    print(head(data))
    for (edge in c("all", "non_neighbor")) {
        if (edge == "all") next
        if (edge != 'all') {
            data <- subset(data, ddist %in% c("inter", "<=10", ">10"))
        }
        chr_set <- unique(unique(data[,"chr"]), unique(data[,"chr.1"]))
        print(length(chr_set))
        mat <- matrix(0, nrow=22, ncol=22)
        rownames(mat) <- paste0("chr", c(1:19, "X", "Y", "M"));
        colnames(mat) <- paste0("chr", c(1:19, "X", "Y", "M"));
        print(rownames(mat))
        print(colnames(mat))
        for (i in 1:dim(data)[1]) {
            x <- data[i, "chr"]; y <- data[i, "chr.1"]
            mat[x, y] <- mat[x,y]+1
            mat[y, x] <- mat[y,x]+1            
        }
        #print(mat)
        tmat <- mat/max(mat)
        png(paste0("chr_enrichment_corrplot_", header, "_", edge, ".png"))
        corrplot(tmat, method="circle")
        dev.off()
        require("scales")
        tmat <- do.call(rbind, lapply(1:dim(mat)[1], function(x) {return(rescale(mat[x,], c(-1, 1)))}))
        print(tmat)
        rownames(tmat) <- paste0("chr", c(1:19, "X", "Y", "M"));
        colnames(tmat) <- paste0("chr", c(1:19, "X", "Y", "M"));
        print(dim(tmat))
        png(paste0("chr_enrichment_corrplot_scaled_", header, "_", edge, ".png"))
        corrplot(tmat, method="circle")
        dev.off()
        lgenes <- GRanges(as.character(data[,"chr"]), IRanges(start=data[,"start"], width=data[,"end"]-data[,"start"]), mcols=data.frame(cov=data[,"count"]), seqlengths=seqlength_mm10)
        rgenes <- GRanges(as.character(data[,"chr.1"]), IRanges(start=data[,"start.1"], width=data[,"end.1"]-data[,"start.1"]), mcols=data.frame(cov=data[,"count"]))
        lgenes$to_gr <- rgenes
        print(chrom)
        print(lgenes)
        print(seqnames(lgenes))
        print(names(seqlength_mm10))
        #seqlengths(lgenes) <- seqlength_mm10
        #next
        print(lgenes)
        p <- ggbio(radius = 30) + circle(lgenes, geom = "link", linked.to = "to_gr", alpha=0.001) + 
         circle(chrom, geom = "ideo", fill = "gray70") +
         circle(chrom, geom = "scale", size = 2) +
         circle(chrom, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
        png(paste0("chr_enrichment_circ_", header, "_", edge, ".png"))
        p
        dev.off()
        ggsave(paste0("chr_enrichment_circ_gg_", header, "_", edge, ".png"))
    }
}

compute.networks <- function(header) {
    require(igraph)
    fname <- paste0("data_frame_", header, "_entrez.tsv")
    data <- read.table(fname, sep="\t", header=T, stringsAsFactor=F)
    for (edge in c("all", "non_neighbor")) {
        print(data[,"ddist"][1:10])
        if (edge != 'all') {
            data <- subset(data, ddist %in% c("inter", "<=10", ">10"))
        }
        print(dim(data))
        nodes <- unique(c(data[,"first"], data[,"second"]))
        quant_value <- rev(quantile(data[,"count"], seq(0.0, 0.75, 0.25)))
        print(quant_value)
        adj_mat <- matrix(0, nrow=length(nodes), ncol=length(nodes))
        nodes_degree <- NULL
        for (i in 1:dim(data)[1]) {
            l <- which(nodes == data[i,"first"])
            r <- which(nodes == data[i, "second"])
            adj_mat[l,r] <- adj_mat[l,r]+1
            adj_mat[r,l] <- adj_mat[r,l]+1
            if (any(data[i, "count"] == quant_value) && (i == dim(data)[1] || data[i+1, "count"] < data[i, "count"])) {
                q <- as.integer(100-which(data[i, "count"] == quant_value)/4*100)
                g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
                V(g)$color=colSums(adj_mat) 
                png(paste0("plot_network_", header, "_", edge, "_", q, ".png"))
                plot(g, vertex.size=5, vertex.label=NA, layout=layout_with_kk)
                dev.off()
                if (is.null(nodes_degree)) {
                    nodes_degree <- cbind(colSums(adj_mat), rep(q, dim(adj_mat)[2]))
                } else {
                    nodes_degree <- rbind(nodes_degree, cbind(colSums(adj_mat), rep(q, dim(adj_mat)[2])))
                }
                print(centr_eigen(g, directed=FALSE)$centralization)
                print(centr_degree(g)$centralization)
                #print(betweenness(g, directed=FALSE))
            }
        }
        colnames(nodes_degree) <- c("degree", "percent")
        nodes_degree <- data.frame(nodes_degree)
        print(head(nodes_degree))
        print(class(nodes_degree))
        nodes_degree[, "percent"] <- factor(nodes_degree[, "percent"], levels=c(0, 25, 50, 75))
        g <- ggplot(data.frame(nodes_degree), aes(x=degree, group=percent, color=percent))+geom_freqpoly()+scale_color_viridis(discrete=TRUE)+xlim(0, 500)
        png(paste0("hist_network_degree_", header, "_", edge, ".png"))
        plot(g)
        dev.off()
        #ggsave(g, "hist_network_degree_gg_", header, "_", edge, ".png")
    }
}

add.go.description <- function(all_mat) {
    require(GO.db)
    annot <- read.table("mm.go.description.csv", sep=",", header=T, stringsAsFactor=F)
    annot = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
    rownames(annot) <- gsub(":", ".", annot[,"GOID"])
    print(dim(all_mat))
    m <- merge(all_mat, annot, by.x=0, by.y=0, all.x=TRUE)
    print(dim(m))
    #m[,"max_pvalue"] <- unlist(sapply(1:dim(m)[1], function(x){return(max(m[x, "pgreat.x"], m[x, "pgreat.x"]))}))
    m <- m[order(m[,"max_pvalue"], decreasing=FALSE),]
    print(head(m))
    return(m)
}

compute.enrichment.go <- function(header) {
    fname <- paste0("data_frame_", header, "_entrez.tsv")
    data <- read.table(fname, sep="\t", header=T, stringsAsFactor=F)
    genes <- colSums(go_eg_mat)
    min_threshold <- 10
    max_threshold <- 200
    all_mat <- NULL
    for (ord in c("a", "b")) {
        if (ord == "a") {
            go_gene_hit <- as.vector(rbind(data[,"legenes"], data[,"regenes"]))
        } else
            go_gene_hit <- as.vector(rbind(data[,"regenes"], data[,"legenes"]))            
        go_gene_hit <- go_gene_hit[go_gene_hit != ""]
        print(go_gene_hit)
        mapped_index <- unlist(sapply(go_gene_hit, function(x){return(which(x == rownames(go_eg_mat)))}))
        print(mapped_index)
        print(dim(mapped_index))
        mapped_index <- data.frame(go_mat_index=mapped_index, rank=1:length(mapped_index))
        enrich_mat <- do.call(rbind, mclapply(which(genes >= min_threshold & genes <= max_threshold), function(x) {
            #GO-gene matrix coordinate according to the order of residual
                gene_index <- which(go_eg_mat[,x] > 0)
                x <- subset(mapped_index, (go_mat_index %in% gene_index))[,"rank"]
                y <- subset(mapped_index, !(go_mat_index %in% gene_index))[,"rank"]
                if (length(x) <= min_threshold)
                    return(c(NA, NA, length(x), length(gene_index)))
                wil_res_less <- wilcox.test(x, y, alternative="less")
                wil_res_greater <- wilcox.test(x, y, alternative="greater")
                return(c(wil_res_greater$p.value, wil_res_less$p.value, length(x), length(gene_index)))
            }, mc.cores=1))
        rownames(enrich_mat) <- colnames(go_eg_mat)[genes >= min_threshold & genes <= max_threshold]
        if (ord == "a") all_mat <- enrich_mat
        else {
            all_mat <- cbind(all_mat, enrich_mat)
        }
    }
    colnames(all_mat) <- c("pgreat.x", "pless.x", "obs.x", "total.x", "pgreat.y", "pless.y", "obs.y", "total.y")
    all_mat <- cbind(all_mat, max_pvalue=apply(all_mat, c(1), function(x){return(max(x["pgreat.x"], x["pgreat.y"]))}))
    add.go.description(all_mat)
    write.table(all_mat, file=paste0("GO_enrichment_", header, ".tsv"), sep="\t", quote=F)
}
# genes.to.name.list <- function(x, en_eg_mat, keywords=c("Name")) {
#     l <- list()
#     for (key in keywords) {
#         ens <- data.frame(Ensembl=unique(unlist(strsplit(x, ",", fixed=TRUE))))
#         m <- merge(ens, en_eg_mat, by="Ensembl")
#         if (dim(m)[1] == 0) return('')
#         l[key] <- paste0(m[,key], collapse=',')
#     }
#     return(l)
# }
map.to.go.genes <- function(header, size=2) {
    data <- read.table(paste0("data_frame_", header, ".tsv"), header=T, stringsAsFactor=F, sep="\t")
    data[,"ddist"] <- dist.to.discrete(log10(data[,"dist"]))
    #data <- subset(data, data[,"ddist"] == "inter" & data[,"ddist"] != ">10") # remove neighbors supposedly from same genes
    #data <- data[1:100,]
    for (h in c("l", "r", "t")[1:2]) {
        lgenes <- unlist(sapply(data[,paste0(h, "genes")], function(x) {genes.to.uniq.entrez(x, en_eg_mat)}))
        lnames <- unlist(sapply(data[,paste0(h, "genes")], function(x) {genes.to.name.list(x, en_eg_mat)}))
        print(c(length(lgenes), length(lnames), dim(data)[1]))
        stopifnot(dim(data)[1] == length(lgenes) && dim(data)[1] == length(lnames))
        l <- list(lgenes, lnames)
        names(l) <- paste0(h, c("genes", "names"))
        tdata <- t(do.call(rbind.data.frame, l))
        colnames(tdata) <- paste0(h, c("egenes", "names"))
        data <- data.frame(data, tdata)
   }
    print(head(data))
    write.table(data, paste0("data_frame_", header, "_entrez.tsv"), sep="\t", quote=F, row.names=FALSE)
}

#test.significance.gene.enrichment("gene", "all", "3C1", bin, 2, sup=200)
dist_factor <<- c("inter", "intra", "<5", "neighbor")

plot.ratio.dist.all <- function(df, header) {
    sup <- df[,"bin"]
    vec <- as.integer(log10(df[,"dist"]))
    vec <- dist.to.discrete(vec)
    print(unique(vec))
    output_data <- data.frame(dist_disc=vec, dist=df[,"dist"], support=sup, count=df[,"count"])
    print(head(output_data))
    g <- ggplot(output_data, aes(support, fill=dist_disc))+geom_histogram()+scale_fill_viridis(discrete=TRUE)
    png(paste0("ratio_hist_", header, ".png"))
    plot(g)
    dev.off()
    print(output_data %>% group_by(support) %>% top_n(5000, wt=count))
    top <- (output_data %>% group_by(support) %>% top_n(5000, wt=count))
    print(head(top))
    stopifnot(dim(top)[1] < dim(output_data)[1])
    g <- ggplot(top, aes(support, fill=dist_disc))+geom_histogram()+scale_fill_viridis(discrete=TRUE)
    png(paste0("ratio_hist_", header, "_5000.png"))
    plot(g)
    dev.off()
    output_data <- data.frame(dist=vec, support=sup, value=df[,1], count=1)
    output_data <- output_data %>% group_by(support) %>% mutate(rank=rank(-value, ties.method="average"))
    output_data[,"dist"] <- factor(output_data[,"dist"], levels=dist_factor)
    pdf = output_data %>% group_by(support, dist) %>% arrange(rank) %>% mutate(cs=cumsum(count))
    print(head(pdf))
    for (f in dist_factor) {
        print(f)
        temp = subset(pdf, dist %in% c(f))
        print(head(temp))
        g <- ggplot(temp, aes(x=rank, y=cs, color=support, group=support))+geom_line(size=2, alpha=0.8)+xlim(0, 50000)+theme_bw()+scale_color_distiller(palette="RdBu")
        png(paste0("ratio_cumsum_", header, "_", f, ".png"))
        plot(g)
        dev.off()
        print('??????')
        g <- ggplot(temp, aes(x=rank, y=cs, color=support, group=support))+geom_line(size=2, alpha=0.8)+ scale_color_distiller(palette="RdBu")+theme_bw()+xlim(0, 5000)+ylim(0, max(subset(temp, temp$rank <= 5000)$cs)+100)
        png(paste0("ratio_cumsum_lim5000_", header, "_", f, ".png"))
        plot(g)
        dev.off()
    }
}

plot.ratio.dist <- function(df, header) {
    if (!any(colnames(df) == "sup"))
        sup <- df[,"bin"]
    else
        sup <- df[,"sup"]
    vec <- as.integer(log10(df[,"dist"]))
    vec <- dist.to.discrete(vec)
    output_data <- data.frame(dist=vec, support=sup)
    mock <- data.frame(dist=c("inter", "neighbor"), support=c(2, 2))
    output_data <- rbind(output_data, mock) 
    print(head(output_data))
    g <- ggplot(output_data, aes(support, fill=dist))+geom_histogram()+scale_fill_viridis(discrete=TRUE)
    png(paste0("ratio_hist_", header, ".png"))
    plot(g)
    dev.off()
    g <- ggplot(output_data %>% group_by(support) %>% top_n(1000), aes(support, fill=dist))+geom_histogram()+scale_fill_viridis(discrete=TRUE)
    png(paste0("ratio_hist_", header, "_1000.png"))
    plot(g)
    dev.off()
    dist_factor = c("inter", "intra", "<5", "neighbor")
    output_data <- data.frame(dist=vec, support=sup, value=df[,1], count=1)
    output_data <- output_data %>% group_by(support) %>% mutate(rank=rank(-value, ties.method="average"))
    output_data[,"dist"] <- factor(output_data[,"dist"], levels=dist_factor)
    pdf = output_data %>% group_by(support, dist) %>% arrange(rank) %>% mutate(cs=cumsum(count))
    print(head(pdf))
    g <- ggplot(pdf, aes(x=rank, y=cs, color=dist))+geom_point()+geom_line()+scale_color_viridis(discrete=TRUE)+xlim(0, 50000)+theme_bw()
    png(paste0("ratio_cumsum_", header, ".png"))
    plot(g)
    dev.off()
    g <- ggplot(pdf, aes(x=rank, y=cs, color=dist))+geom_point()+geom_line()+scale_color_viridis(discrete=TRUE)+theme_bw()+xlim(0, 10000)+ylim(0, max(subset(pdf, pdf$rank <= 10500)$cs)+100)
    png(paste0("ratio_cumsum_lim5000_", header, ".png"))
    plot(g)
    dev.off()
}

map.to.genes <- function(df, header) {
    value_column <- "count"
    if (!any(colnames(df) == value_column))
        value_column <- "observed"
    left <- bin_ann[df[,"first"],]
    right <- bin_ann[df[,"second"],]
    print(head(left))
    print(head(right))
    lranges <- makeGRangesFromDataFrame(left[,c("chr", "start", "end")])
    lgenes <- assign.gene.names(lranges)
    rranges <- makeGRangesFromDataFrame(right[,c("chr", "start", "end")])
    rgenes <- assign.gene.names(rranges)
    stopifnot(dim(df)[1] == length(lgenes) && dim(df)[1] == length(rgenes))
    adf <- data.frame(df, left[,c("chr", "start", "end")], right[,c("chr", "start", "end")], lgenes=lgenes, rgenes=rgenes)
    adf <- adf[order(adf[,which(colnames
                                (adf) == value_column)], decreasing=TRUE),]
    write.table(adf,  paste0("data_frame_", header, ".tsv"), sep="\t", quote=F, row.names=FALSE)
}


map.all.bins.to.genes <- function(bin) {
    require(annotatr)
    ranges <- makeGRangesFromDataFrame(bin_ann[,c("chr", "start", "end")])
    genes <- assign.gene.names(ranges)
    adf <- data.frame(bin_ann, genes=genes)
    adf <- adf[order(adf$cov, decreasing=TRUE),]
    write.table(adf,  paste0("data_frame_whole_", bin, ".tsv"), sep="\t", quote=F, row.names=FALSE)
    
}

test.significance.gene.enrichment <- function(row_or_column, chr, batch, bin, size, rel=-1, sup=-1, plot=TRUE) {
    if (sup > 0) {
        header <- paste(row_or_column, chr, batch, bin, "size", size, "sup", sup, sep="_")
    } else if (rel > 0) {
        header <- paste(row_or_column, chr, batch, bin, "size", size, rel, sep="_")
    } else {
        return(0)
    }
    fname <- paste0(data_dir, "lcm_result/", "lcm_", header, ".txt")
    #df <- read.itemset.data(fname, sup)
    compute.pvalue.of.gene.enrichment(header)
    annotate.genic.region(header)
}
annotate.genic.region <- function(header) {
    require(annotatr)
    bin_ann <- read.table(paste0(data_dir, "global_ng_list/", "global_bin_ng_", bin, ".csv"), header=T, sep=",")
    regions <- makeGRangesFromDataFrame(cbind(bin_ann[,c("chr", "start", "end", "cov")], batch=rep(1, dim(bin_ann)[1])), keep.extra.columns = TRUE)
    df <- rbind(s[,c("chr", "start", "end", "observed")], data.frame(chr=s[,"chr.1"], start=s[,"start.1"], end=s[,"end.1"], observed=s[,"observed"]))
    df <- df[!duplicated(df[,c(1:3)]),]
    df[,"chr"] <- unlist(sapply(df[,"chr"], function(x){return(paste0("chr", x))}))
    sregions <- makeGRangesFromDataFrame(cbind(df, batch=rep(1, dim(df)[1])), keep.extra.columns = TRUE)
    # require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    # enhancer_file <- "/data/rkawaguc/scATAC/data/190711/enhancer/mouse_permissive_enhancers_phase_1_and_2_mm10.bed"
    # read_annotations(con = enhancer_file, genome = 'mm10', name = 'enh', format = 'bed')

    annots = c('mm10_basicgenes', 'mm10_genes_intergenic',
        'mm10_genes_intronexonboundaries')
    annotations = build_annotations(genome = 'mm10', annotations = annots)
    # annotate.genic.region <- function(header) {
        #annotations_2 = build_enhancer_annots(genome = c("mm10"))
        # Build the annotations (a single GRanges object)
        tfbs_codes = c('AH6120', 'AH6125')
        # Fetch ah_codes from AnnotationHub and create annotations annotatr understands
    #     build_ah_annots(genome = 'mm10', ah_codes = tfbs_codes, annotation_class = 'tfbs')
    #     ah_names = c('mm10_tfbs')
    #     print(annotatr_cache$get('mm10_tfbs'))
    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(
        regions = regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)
    s_annotated = annotate_regions(
        regions = sregions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)
    # A GRanges object is returned
    print(dm_annotated)
    # Coerce to a data.frame
    df_dm_annotated = data.frame(dm_annotated)
    # See the GRanges column of dm_annotaed expanded
    print(head(df_dm_annotated))
    # # Randomize the input regions
    # dm_random_regions = randomize_regions(
    #     regions = dm_regions,
    #     allow.overlaps = TRUE,
    #     per.chromosome = TRUE)

    # # Annotate the random regions using the same annotations as above
    # # These will be used in later functions
    # dm_random_annotated = annotate_regions(
    #     regions = dm_random_regions,
    #     annotations = annotations,
    #     ignore.strand = TRUE,
    #     quiet = TRUE)
    dm_annsum = summarize_annotations(
        annotated_regions = dm_annotated,
        quiet = TRUE)
    print(dm_annsum)
    s_annsum = summarize_annotations(
        annotated_regions = s_annotated,
        quiet = TRUE)
    print(s_annsum)
    #     dm_catsum = summarize_categorical(
    #         annotated_regions = dm_annotated,
    #         by = c('annot.type', 'DM_status'),
    #         quiet = TRUE)
    #     print(dm_catsum)
    annots_order = c(
        'mm10_genes_1to5kb',
        'mm10_genes_promoters',
        'mm10_genes_5UTRs',
        'mm10_genes_exons',
        'mm10_genes_intronexonboundaries',
        'mm10_genes_introns',
        'mm10_genes_3UTRs',
        'mm10_genes_intergenic')
    png(paste0("plot_annotation_count_whole_", header, ".png"))
        dm_vs_kg_annotations = plot_annotation(
            annotated_regions = dm_annotated,
            annotation_order = annots_order,
            plot_title = '# of Sites',
            x_label = 'knownGene Annotations',
            y_label = 'Count')
    dev.off()
    #     print(dm_vs_kg_annotations)
    #     annots_order = c(
    #         'hg19_custom_ezh2',
    #         'hg19_H3K4me3_Gm12878',
    #         'hg19_genes_1to5kb',
    #         'hg19_genes_promoters',
    #         'hg19_genes_5UTRs',
    #         'hg19_genes_exons',
    #         'hg19_genes_intronexonboundaries',
    #         'hg19_genes_introns',
    #         'hg19_genes_3UTRs',
    #         'hg19_genes_intergenic')
    png(paste0("plot_annotation_count_comp", header, ".png"))
    dm_vs_kg_annotations_wrandom = plot_annotation(
        annotated_regions = dm_annotated,
        annotated_random = s_annotated,
        annotation_order = annots_order,
        plot_title = 'Dist. of Sites',
        x_label = 'Annotations',
        y_label = 'Proportion')
    print(dm_vs_kg_annotations_wrandom)
    dev.off()
    print(head(dm_annotated))
    print(s_annotated)
    png(paste0("plot_annotation_density_comp", header, ".png"))
    dm_vs_cpg_cat_random = plot_categorical(
        annotated_regions = dm_annotated, annotated_random = s_annotated,
        x='batch', fill='annot.type',
        x_order = c(1), fill_order = annots_order, position='fill',
        plot_title = 'Annotation Proportions',
        legend_title = 'Annotations',
        x_label = 'batch',
        y_label = 'Proportion')
    print(dm_vs_cpg_cat_random)
    dev.off()
    if (!grepl("sup", header, fixed = TRUE))
        return()
    png(paste0("plot_coverage_density_comp", header, ".png"))
    dm_vs_regions_annot = plot_numerical(
        annotated_regions = dm_annotated,
        x = 'cov',
        facet = 'annot.type',
        facet_order = c('mm10_genes_1to5kb','mm10_genes_promoters',
            'mm10_genes_5UTRs','mm10_genes_3UTRs',
            'mm10_genes_intergenic'),
        bin_width = 5,
        plot_title = 'Group 0 Region Methylation In Genes',
        x_label = 'Group 0')
    dev.off()
#     # View a heatmap of regions occurring in pairs of annotations
#     annots_order = c(
#         'hg19_custom_ezh2',
#         'hg19_H3K4me3_Gm12878',
#         'hg19_genes_promoters',
#         'hg19_genes_5UTRs',
#         'hg19_genes_exons',
#         'hg19_genes_introns',
#         'hg19_genes_3UTRs',
#         'hg19_genes_intergenic')
#     dm_vs_coannotations = plot_coannotations(
#         annotated_regions = dm_annotated,
#         annotation_order = annots_order,
#         axes_label = 'Annotations',
#         plot_title = 'Regions in Pairs of Annotations')
#     print(dm_vs_coannotations)
#     dm_vs_regions_annot = plot_numerical(
#         annotated_regions = dm_annotated,
#         x = 'mu0',
#         facet = 'annot.type',
#         facet_order = c('hg19_genes_1to5kb','hg19_genes_promoters',
#             'hg19_genes_5UTRs','hg19_genes_3UTRs', 'hg19_custom_ezh2',
#             'hg19_genes_intergenic', 'hg19_cpg_islands'),
#         bin_width = 5,
#         plot_title = 'Group 0 Region Methylation In Genes',
#         x_label = 'Group 0')
#     print(dm_vs_regions_annot)
# }
}
# list of annotation hub
# library(AnnotationHub)
# ah <- AnnotationHub()
# ah <- subset(ah, species == "Mus musculus")
# hist <- display(ah)
test.wilcoxon.vector <- function(vec, target) {
    x <- which(vec == target)
    y <- which(vec != target)
    w <- wilcox.test(x, y, alternative="less")
    return(w)
}
compute.pvalue.of.gene.enrichment <- function(header) {
    #header <- paste0("gene_all_3C1_", bin, "_size_2_200")
    s <- read.table(paste0("data_frame_", header, ".tsv"), sep="\t", header=T, stringsAsFactor=F)
    a <- read.table(paste0("data_frame_whole_", bin, ".tsv"), sep="\t", header=T, stringsAsFactor=F)
    top_genes <- apply(s, c(1), function(x){
        if(nchar(x[["lgenes"]]) > 0 && nchar(x[["rgenes"]]) > 0) return("Both");
        if (nchar(x[["lgenes"]]) > 0 || nchar(x[["rgenes"]]) > 0) return("Either");
        return("None");
    })
    top_genes_any <- top_genes
    top_genes_any[top_genes_any == "Both"] <- "Either"
    df <- data.frame(gene=top_genes, any_gene=top_genes_any, count=1, rank=1:length(top_genes))
    df[,"gene"] <- factor(df[,"gene"], levels=c("Both", "Either", "None"))
    df[,"any_gene"] <- factor(df[,"any_gene"], levels=c("Either", "None"))
    pdf <- df %>% group_by(gene) %>% mutate(cs_gene = cumsum(count))
    pdf <- pdf %>% group_by(any_gene) %>% mutate(cs_any = cumsum(count))
    print(head(pdf))
    g <- ggplot(pdf, aes(x=rank, y=cs_gene, color=gene))+geom_point()+geom_line()+scale_color_viridis(discrete=TRUE)+theme_bw()
    png(paste0("ratio_cumsum_gene_annot_", header, ".png"))
    plot(g)
    dev.off()
    g <- ggplot(pdf, aes(x=rank, y=cs_any, color=any_gene))+geom_point()+geom_line()+scale_color_viridis(discrete=TRUE)+theme_bw()
    png(paste0("ratio_cumsum_gene_annot_any_", header, ".png"))
    plot(g)
    dev.off()
    whole_genes <- apply(a, c(1), function(x){
        if(nchar(x[["genes"]]) > 0) return("Gene"); return("None");
    })
    sink(paste0("statistics_result_", header, ".txt")) 
    print(test.wilcoxon.vector(top_genes_any, "Either"))
    print(test.wilcoxon.vector(whole_genes, "Gene"))
    uniq_list <- rbind(cbind(s[,"first"], s[,"lgenes"]), cbind(s[,"second"], s[, "rgenes"]))
    uniq_list <- uniq_list[!duplicated(uniq_list),]
    print(dim(uniq_list))
    #print(head(uniq_list))
    uniq_list <- unlist(sapply(uniq_list[,2], function(x){if(nchar(x) > 0)return("Gene"); return("None")}))
    N <- length(whole_genes); M <- length(which(whole_genes == "Gene")); n <- length(top_genes_any); m <- length(which(top_genes_any == "Either"))
    print(c(N, M, n, m))
    print(1.0-phyper(m-1, M, N-M, n))
    N <- length(whole_genes); M <- length(which(whole_genes == "Gene")); n <- length(uniq_list); m <- length(which(uniq_list == "Gene"))
    print(c(N, M, n, m))
    print(1.0-phyper(m-1, M, N-M, n))
    sink()
}
# assign.gene.names <- function(query=c(), dist=50000) {
#     if (length(query) == 0) {
#         df <- data.frame(chr=c("1", "2", "3"), start=c(1000, 10000, 100000), end=c(10000000, 11000, 101000))
#         query <- makeGRangesFromDataFrame(df)
#     }
# txdb <- makeTxDbFromEnsembl(organism="Mus musculus",
#                     release=NA,
#                     server="ensembldb.ensembl.org",
#                     username="anonymous", password=NULL, port=0L,
#                     tx_attrib=NULL)
#     subject <- promoters(txdb, upstream=dist, downstream=dist)
#     ov <- findOverlaps(query, subject, minoverlap=1L, type="any", select="all", ignore.strand=TRUE)
#     #print(ov)
#     ov <- data.frame(ov)
#     print(paste(names(subject)[ov[which(ov[,1] == 1),2]], collapse=","))
#     return(unlist(sapply(1:length(query), function(x){return(paste(names(subject)[ov[which(ov[,1] == x),2]], collapse=","))})))
# }



construct.go.mat <- function(fname) {
    require(org.Mm.eg.db)
    x <- org.Mm.egGO
    # Get the entrez gene identifiers that are mapped to a GO ID
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    gos <- unlist(lapply(names(xx), function(x){return(names(xx[[x]]))}))
    gos <- gos[order(gos)]
    gos <- gos[!duplicated(gos)]
    genes <- names(xx)
    genes <- genes[order(genes)]
    genes <- genes[!duplicated(genes)]
    matrix <- matrix(0, nrow=length(genes), ncol=length(gos))
    rownames(matrix) <- genes
    colnames(matrix) <- gos
    for (name in names(xx)) {
        for (go_candidate in names(xx[[name]])) {
            if (xx[[name]][[go_candidate]]$Evidence != "IEA") {
                matrix[name, go_candidate] <- 1
            }
        }
    }
    write.table(matrix, file=fname, sep=",", quote=F)
}

construct.go.desc <- function(fname) {
    require(GO.db)
    xx <- as.list(GOBPCHILDREN)
    vec <- names(xx)
    df <- data.frame(GO=vec, Term=unlist(mclapply(list(vec), function(x) {
       return(Term(x)) 
    }, mc.cores=8)))
    rownames(df) <- vec
    write.table(df, file=fname, sep=",", quote=T)
}

convert.en2eg <- function() {
    require(org.Mm.eg.db)
    require(biomaRt)
    xx <- mappedkeys(org.Mm.egENSEMBLTRANS2EG)
    print(length(xx))
    print(xx[1:5])
    mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast")
    #mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl", mirror = "useast"))
    #print(listAttributes(mart))
    genes <- getBM(
      filters="ensembl_transcript_id",
      attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "entrezgene_id"),
      values=xx, 
      mart=mart)
    colnames(genes) <- c("Ensembl", "En_gene", "Name", "Entrez")
    write.table(genes, file="mm.en.et.csv", sep=",", quote=F)
}

# construct.go.mat("mm.go.mat.csv")
# construct.go.desc("mm.go.description.csv")
# convert.en2eg("mm.en.et.csv")

genes.to.uniq.entrez <- function(x, en_eg_mat) {
    ens <- data.frame(Ensembl=unique(unlist(strsplit(x, ",", fixed=TRUE))))
    m <- merge(ens, en_eg_mat, by="Ensembl")
    if (dim(m)[1] == 0) return('')
    l <- list()
    entrez <- as.character(unique(m[,"Entrez"]))
    if (length(entrez) > 1) {
        return(names(sort(table(entrez), decreasing=TRUE))[1])
    }
    return(entrez[1])
}
select.cann.index <- function(bin_ann, index) {
    require(dplyr)
    require(multidplyr)
    annots_order = c(
    'mm10_genes_promoters',
    'mm10_genes_5UTRs',
    'mm10_genes_exons',
    'mm10_genes_intronexonboundaries',
    'mm10_genes_introns',
    'mm10_genes_3UTRs',
    'mm10_genes_1to5kb',
    'mm10_genes_intergenic')
    CORES = 10
    cluster <- new_cluster(CORES)
    if (any("annot.type" == colnames(bin_ann))) {
        bin_ann[,"annot.type"] <- factor(bin_ann[,"annot.type"], levels=annots_order)
        bin_ann[,"annot.type"] <- as.numeric(bin_ann[,"annot.type"])
    }
    df <- bin_ann %>% group_by_(index) %>% partition(cluster)
    df <- df %>% summarise(chr=first(chr), start=min(start, na.rm=TRUE), end=min(start, na.rm=TRUE)+bin, gene_id=na.omit(gene_id)[1], gene_name=na.omit(gene_name)[1], dist=min(abs(dist), na.rm=TRUE), enhancer=enhancer[1], id_proximal=na.omit(id_proximal)[1], id_order_distal=na.omit(id_order_distal)[1], annot.type=annots_order[min(annot.type, na.rm=TRUE)]) %>% collect()
    return(df)
}

compute.dist.pairs <- function(mat, bin_ann, index) {
    CORES <- 10
    print('compute dist pairs')
    left <- bin_ann[unlist(sapply(mat[,1], function(x){return(which(x == bin_ann[,1])[1])})),]
    print(head(left))
    right <- bin_ann[unlist(sapply(mat[,2], function(x){return(which(x == bin_ann[,1])[1])})),]
    print(head(right))
    dist <- unlist(mclapply(1:dim(mat)[1], function(x) {
        if (left[x,"chr"] == right[x,"chr"]) {
            if (left[x, "start"] < right[x, "start"])
                return(right[x, "start"]-left[x, "end"])
            else
                return(left[x, "start"]-right[x, "end"])                
        } else {
            return(0.001)
        }
    }, mc.cores=CORES))
    print('Distant signals')
    print(head(cbind(left[dist >= 1000 | dist < 1 ,c("chr", "start", "end")], right[dist >= 1000 | dist < 1,c("chr", "start", "end")])))
    colnames(left) = paste0("l_", colnames(left))
    colnames(right) = paste0("r_", colnames(right))
    return(cbind(dist=dist, cbind(left, right)))
}

dist.to.discrete <- function(dist) {
    return(unlist(sapply(dist, function(x){
        if(x < 0)return("inter");
        if(x == 0)return("neighbor");
        if(x <= 5)return("<5");
        return("intra");
        # if (x <= 10)return("<=10");
        # return(">10")
    })))
}


read.itemset.data <- function(fname, bin_ann, index, sup, max_limit=50000) {
    if (sup > 0) {
        df <- read.table(fname, header=F)
        df[,1] <- unlist(sapply(df[,1], function(x){return(gsub("\\(|\\)", "", x))}))
        colnames(df) <- c("count", "first", "second")
        df <- cbind(df, compute.dist.pairs(df[,2:3], bin_ann, index))
        df <- transform(df, count = as.numeric(count), first=as.numeric(first), second=as.numeric(second), dist=as.numeric(dist))
        print(sapply(df, class))
        print(head(df))
    } else {
        df <- read.table(text =  gsub(",|\\[|\\]"," ", readLines(fname)), header=F)
        print(head(df))
        colnames(df) <- c("observed", "expected", "first", "second")
        df <- cbind(df, compute.dist.pairs(df[,c("first", "second")], bin_ann, index))      
        df <- transform(df, observed = as.numeric(observed), expected=as.numeric(expected), first=as.numeric(first), second=as.numeric(second), dist=as.numeric(dist))
    }
    df <- df[order(df$count, decreasing=TRUE),]
    return(df[1:min(dim(df)[1], max_limit),])
}

plot.scatter.distance <- function(fname, header, bin_ann, index, sup) {
    print('plot_scatter')
    df <- read.itemset.data(fname, bin_ann, index, sup)
    print('read end')
    if (sup > 0) {
        png(paste0("abs_dist_", header, ".png"))
        g <- ggplot(df, aes(x=dist, y=count))+geom_point(alpha=0.3)+theme_bw()+scale_x_continuous(trans='log10')
        plot(g)
        dev.off()
        png(paste0("abs_dist_", header, "_inter.png"))
        g <- ggplot(subset(df, df$dist >= 1), aes(x=dist))+geom_histogram()+theme_bw()+scale_x_continuous(trans='log10')
        dev.off()
    } else {
        df <- cbind(df, log10_dist=log10(df[,"dist"]))
        png(paste0("lift_dist_", header, ".png"))
        g <- ggplot(df, aes(x=log10_dist, y=observed))+geom_point(alpha=0.3)+theme_bw()
        plot(g)
        dev.off()
        df <- cbind(df, disc_dist=dist.to.discrete(df$log10_dist))
        df[,"disc_dist"] <- factor(df[,"disc_dist"], levels=c("inter", "intra", "<5", "neighbor"))
        png(paste0("obs_lift_distcol_", header, ".png"))
        g <- ggplot(df, aes(x=expected, y=observed, color=disc_dist))+geom_point(alpha=0.3)+theme_bw()+scale_color_viridis(discrete=TRUE)
        plot(g)
        dev.off()
    }
    df <- df[order(df$count, decreasing=TRUE),]
    plot.ratio.dist(df, header)
    return(df)
}

visualize.lcm.result.pair <- function(file_list, plot=TRUE, sup=TRUE) {
    require(tools)
    a <- read.table(file_list, header=F, stringsAsFactor=FALSE, sep=" ")
    colnames(a) <- c("file", "bin", "threshold", "batch", "annfile")
    fall <- paste0(file_path_sans_ext(basename(file_list)), "_signal.rds")
    if (plot) {
        all_data <- NULL
        for (i in 1:dim(a)[1]) {
            file = a[i,1]
            header = paste(a[i,4], a[i,2], a[i,3], sep="_")
            index = paste0('global_index_', a[i,2])
            if (a[i,2] == 1000) index = 'global_index'
            bin_ann <- read.table(a[i,5], header=T, sep=",")
            bin_ann[,'chr'] <- gsub('chr', '', bin_ann[,'chr'], fixed=TRUE)
            fann <- paste0(file_path_sans_ext(basename(file_list)), "_", i, ".rds")
            if (file.exists(fann)) uniq_df <- readRDS(fann)
            else {
                uniq_df <- as.data.frame(select.cann.index(bin_ann, index))
                print(head(uniq_df))
                uniq_df[,"start"] = as.integer(uniq_df[,"start"]/a[i,2])*a[i,2]+1
                uniq_df[,"end"] = uniq_df[,"start"]+a[i,2]-1
                saveRDS(uniq_df, fann)
            }
            df <- plot.scatter.distance(file, header, uniq_df, index, sup)
            df <- df[,colnames(df) != index & colnames(df) != paste0("r_", index) & colnames(df) != paste0("l_", index)]
            all_data <- rbind(all_data, cbind(df, bin=rep(a[i,2], dim(df)[1])))
        }
        saveRDS(all_data, fall)
    } else if (file.exists(fall)) {
        all_data <- readRDS(all_data)
        header = paste0(file_path_sans_ext(basename(file_list)), "_signal")
        plot.ratio.dist.all(all_data, header)

    } else {
        df <- read.itemset.data(fname, sup)
        map.to.genes(df, header)
    }
}


args = commandArgs(trailingOnly=TRUE)
visualize.lcm.result.pair(args[1])