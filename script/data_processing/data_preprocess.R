
marker_dir = file.path(getwd(), "../../data/cluster_annotation/")
argvs = commandArgs(trailingOnly=TRUE)
GSE = argv[1]
data_dir = file.path(argv[2], GSE)
meta_dir = file.path(getwd(), "../../data/metadata/", GSE)
if (!file.exists(data_dir)){
    dir.create(data_dir)
}

if (GSE == 'GSE111586') {
    a <- read.table(paste0(meta_dir, "GSE111586/cell_metadata.txt"), header=T, sep="\t", stringsAsFactor=F)
    fname_list <- c("Wholebrain1", "Wholebrain2", "Prefrontal")
    start <- 0
    options(scipen=99)
    for (head in fname_list) {
        print(start)
        dir <- data_dir
        b <- read.table(paste0(data_dir, "GSE111586_cell_ng_", head, ".csv"), header=T, sep=",", stringsAsFactor=F)
        m <- subset(b, !is.na(b[,"local_index"]))
        colnames(m)[which(colnames(m) == "tsne_1")] = "tsne1"
        colnames(m)[which(colnames(m) == "tsne_2")] = "tsne2"
        colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
        m <- m[,colnames(m) != "X"]
        m$local_index <- as.integer(as.character(m$local_index))
        m <- m[order(m$local_index),]
        m$global_index <- as.numeric(m$local_index)+start
        rownames(m) <- m$local_index
        print(head(m))
        stopifnot(any(colnames(m) == "global_index"))
        dict_data = read.table(paste0(marker_dir, "GSE111586_icluster_celltype_annotation.csv"), header=T, sep=",")
        dict = dict_data[,2]
        names(dict) = dict_data[,1]
        m <- cbind(m, celltype=unlist(sapply(m[,"cell_label"], function(x){
            if(is.na(x)) return("NA");
            return(as.character(dict[x]))
        })))
        m[m == ""] <- NA
        write.table(m, paste0(dir, "GSE111586_cell_ng_", head, "_meta.csv"), sep=",", quote=F)
        start <- start + max(m$local_index, na.rm=TRUE) + 1
    }
}
if (GSE == 'GSE123576') {
    dir <- data_dir
    a <- read.table(paste0(dir, "GSE123576_cell_ng_mousebrain.csv"), sep=",", header=T)
    dict = c(rep("OT", 5), rep("IN", 5), rep("EX", 17))
    names(dict) =  c(c(1, 26, 24, 6, 27), c(9, 12, 14, 5, 11), c(2, 3, 4, 7, 8, 10, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25))
    a[,"celltype"] = unlist(sapply(a[,"cluster"], function(x) {
        if (is.na(x)) return("NA")
        return(getElement(dict, as.character(x)))
    }))
    write.table(a, paste0(dir, "GSE123576_cell_ng_mousebrain_meta.csv"), sep=",", quote=F)
}

if (GSE == 'GSE127257') {
    a <- read.table(file.path(meta_dir, "GSE127257_metadata.csv"), sep=",", header=T)
    fname_list <- c("Ts11", "Ts12", "Ts21", "Ts22", "N11", "N12", "N21", "N22")
    batch_list <- c("Ts1_R1", "Ts1_R2", "Ts2_R1", "Ts2_R2", "2n1_R1", "2n1_R2", "2n2_R1", "2n2_R2")
    for (i in 1:length(fname_list)) {
        head = fname_list[i]
        dir <- data_dir
        temp <- subset(a, a$sample == batch_list[i])
        b <- read.table(paste0(dir, "GSE127257_cell_ng_", head, ".csv"), header=T, sep=",")
        m <- merge(b, temp, by="barcode", all.x=TRUE)
        m <- m[,colnames(m) != "X"]
        m <- m[order(m$global_index),]
        write.table(m, paste0(dir, "GSE127257_cell_ng_", head, "_meta.csv"), sep=",", quote=F)
    }
}

add.cell.type <- function(md) {
    md[,"celltype"] = "NA"
    dict=c(rep("EX", 10), rep("IN", 4), rep("MG"), rep("OG", 3), "AC", "OT", "OT")
    names(dict)=c(c("E2Rasgrf2", "E3Rmst", "E3Rorb", "E4Il1rapl2", "E4Thsd7a","E5Galnt14", "E5Parm1", "E5Sulf1", "E5Tshz2", "E6Tle4"), c("InN", "InP", "InS", "InV", "Mic", "OPC", "OliI", "OliM", "Ast", "Clau", "Peri")) 
    md[,"celltype"] = unlist(sapply(md[,"Ident"], function(x) {
        x=as.character(x)
        if (is.na(x)) return("NA")
        if (all(names(dict) != x)) return("NA")
        return(getElement(dict, x))
    }))
    example <- md[,c("celltype", "Ident")]
    print(example[!duplicated(example),])
    return(md)
}

if (GSE == 'GSE126074') {
    a <- read.table(file.path(meta_dir, "AdCortex_ann_clust.csv", sep=",", header=T))
    dir <- data_dir
    b <- read.table(paste0(dir, "GSE126074_cell_ng_AdCortex.csv"), header=T, sep=",")
    stopifnot(all(a[,1] == b[,2]))
    m <- cbind(b, a[,-c(1)])
    m <- m[,colnames(m) != "X"]
    colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    m <- m[order(m$global_index),]
    d <- readRDS(file.path(meta_dir, "Ad_BrainCortex_metadata.rds"))
    d <- cbind(d, barcode=rownames(d))
    d <- d[,c("barcode", "Ident")]
    md <- merge(m, d, by="barcode", all.x=TRUE)
    md <- md[order(md$global_index),]
    md <- add.cell.type(md)
    write.table(md, paste0(dir, "GSE126074_cell_ng_", "AdCortex", "_meta.csv"), sep=",", quote=F)


    b <- read.table(paste0(dir, "GSE126074_cell_ng_AdCortexr.csv"), header=T, sep=",")
    b <- b[order(b$global_index_chrom),]
    stopifnot(all(a[,1] == b[,2]))
    m <- cbind(b, a[,-c(1)])
    m <- m[,colnames(m) != "X"]
    colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    m <- m[order(m$global_index),]
    d <- readRDS(file.path(meta_dir, "Ad_BrainCortex_metadata.rds"))
    d <- cbind(d, barcode=rownames(d))
    d <- d[,c("barcode", "Ident")]
    md <- merge(m, d, by="barcode", all.x=TRUE)
    md <- md[order(md$global_index),]
    md <- add.cell.type(md)
    write.table(md, paste0(dir, "GSE126074_cell_ng_", "AdCortexr", "_meta.csv"), sep=",", quote=F)
}

if (GSE == 'GSE130399') {
    gse_id <- c("GSE1303990", "GSE1303991")
    meta <- c("Adult_CTX", "Fetal_FB")
    i = 1
    dir <- file.path(paste0(data_dir, '0'))
    a <- read.table(file.path(data_dir, meta[i], paste0(meta[i], "_embed.csv")), sep=",", header=T, stringsAsFactor=F)
    b <- read.table(file.path(data_dir, paste0(gse_id[i], "_cell_ng_", c("Actxr", "Fb")[i], ".csv")), header=T, sep=",", stringsAsFactor=F)
    stopifnot(all(a[,1] == b[,2]))
    m <- cbind(b, a[,-c(1)])
    m <- m[,colnames(m) != "X"]
    colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    colnames(m)[which(colnames(m) == "Cluster")] = "cluster"
    m <- m[order(m$global_index),]
    png(paste0("test.png"))
    plot(m[,8], m[,9], col=m[,7])
    dev.off()
    d <- read.table(file.path(data_dir, meta[i], paste0(meta[i], "_embed_meta.csv")), sep=",", header=T)
    md <- merge(m, d, by="cluster", all.x=TRUE)
    md <- md[order(md$global_index),]
    md <- md[,colnames(md) != "batch"]

    colnames(md)[which(colnames(m) == "Rep")] = "batch"
    write.table(md, file.path(dir, paste0("GSE1303990_cell_ng_", "Actxr", "_meta.csv")), sep=",", quote=F)

    print(paste0(dir, '/', gse_id[i], "_cell_ng_", c("Actx", "Fb")[i], ".csv"))
    b <- read.table(paste0(dir, '/', gse_id[i], "_cell_ng_", c("Actx", "Fb")[i], ".csv"), header=T, sep=",")
    stopifnot(all(a[,1] == b[,2]))
    m <- cbind(b, a[,-c(1)])
    m <- m[,colnames(m) != "X"]
    colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    colnames(m)[which(colnames(m) == "Cluster")] = "cluster"
    m <- m[order(m$global_index),]
    md <- merge(m, d, by.x="cluster", all.x=TRUE)
    md <- md[order(md$global_index),]
    md <- md[,colnames(md) != "batch"]
    colnames(md)[which(colnames(m) == "Rep")] = "batch"
    write.table(md, file.path(dir, paste0("GSE1303990_cell_ng_", "Actx", "_meta.csv")), sep=",", quote=F)

    # i = 2

    # dir <-paste0("/data/rkawaguc/data/190813_meta_scATAC/processed/", gse_id[i])
    # a <- read.table(paste0("/data/rkawaguc/data/190813_meta_scATAC/GSE130399/Fetal_FB/", meta[i], "_embed.csv"), sep=",", header=T, stringsAsFactor=F)
    # b <- read.table(paste0(dir, '/', gse_id[i], "_cell_ng_", c("Actxr", "Fbr")[i], ".csv"), header=T, sep=",", stringsAsFactor=F)
    # stopifnot(all(a[,1] == b[,2]))
    # m <- cbind(b, a[,-c(1)])
    # m <- m[,colnames(m) != "X"]
    # colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    # colnames(m)[which(colnames(m) == "Ident")] = "cluster"
    # m <- m[order(m$global_index),]
    # d <- read.table("/data/rkawaguc/data/190813_meta_scATAC/GSE130399/Fetal_FB/Fetal_FB_embed_meta.csv", sep=",", header=T)
    # md <- merge(m, d, by="cluster", all.x=TRUE)
    # md <- md[order(md$global_index),]
    # md <- md[,colnames(md) != "batch"]
    # md <- cbind(md, batch=unlist(sapply(1:dim(md)[1], function(x){return(paste0(md[x,"Stage"], '_', md[x,"Rep"]))})))
    # write.table(md, paste0(dir, "/GSE1303991_cell_ng_", "Fbr", "_meta.csv"), sep=",", quote=F)

    # print(paste0(dir, '/', gse_id[i], "_cell_ng_", c("Actx", "Fb")[i], ".csv"))
    # b <- read.table(paste0(dir, '/', gse_id[i], "_cell_ng_", c("Actx", "Fb")[i], ".csv"), header=T, sep=",")
    # stopifnot(all(a[,1] == b[,2]))
    # m <- cbind(b, a[,-c(1)])
    # m <- m[,colnames(m) != "X"]
    # colnames(m)[which(colnames(m) == "barcodes")] = "barcode"
    # colnames(m)[which(colnames(m) == "Ident")] = "cluster"
    # m <- m[order(m$global_index),]
    # md <- merge(m, d, by="cluster", all.x=TRUE)

    # md <- md[order(md$global_index),]
    # md <- md[,colnames(md) != "batch"]
    # md <- cbind(md, batch=unlist(sapply(1:dim(md)[1], function(x){return(paste0(md[x,"Stage"], '_', md[x,"Rep"]))})))
    # write.table(md, paste0(dir, "/GSE1303991_cell_ng_", "Fb", "_meta.csv"), sep=",", quote=F)

    # stopifnot(all(a[,1] == b[,2]))
}

if (GSE == 'BICCN') {
    dir <- data_dir
    b <- read.table(paste0(dir, "ATAC/", "ATAC.cell_tidy_data.csv"), sep=",", header=T, stringsAsFactor=F)
    fname = c("CEMBA171206_3C", "CEMBA171207_3C", "CEMBA171212_4B", "CEMBA171213_4B", "CEMBA180104_4B", "CEMBA180409_2C", "CEMBA180410_2C", "CEMBA180612_5D", "CEMBA180618_5D")
    batch = c("3C1", "3C2", "4B3", "4B4", "4B5", "2C6", "2C7", "5D8", "5D9")
    celltype <- c("L6.CT", "L23.b", "ASC", "L5.IT.b", "Other", "L4", "Pv", "NP", "MGC", "L5.IT.a", "L23.a", "OPC", "OGC", "L5.PT", "Sst", "L6.IT", "CGE", "L23.c", "Endo", "Smc")
    celltype_neuron <- c("EX", "EX", "AC", "EX", "OT", "EX", "IN", "EX", "OT", "OT", "EX", "OG", "OG", "EX", "IN", "EX", "IN", "EX", "OT", "OT")
    b[,"batch"] <- unlist(sapply(b[,"sample"], function(x){return(batch[which(fname == x)][1])}))
    b[,"full_barcode"] <- unlist(sapply(1:dim(b)[1], function(x){return(paste0(b[x,"batch"], ".", b[x, "barcode"]))}))
    colnames(b)[colnames(b) == "tsne_1"] = "tsne_2"
    colnames(b)[colnames(b) == "tsne_0"] = "tsne_1"
    b <- b[,colnames(b) != "barcode"]
    b[,"celltype"] <- unlist(sapply(b[,"MajorCluster"], function(x){return(celltype_neuron[which(celltype == x)[1]])}))
    for (f in list.files(paste0(dir, "from_snap/global_ng_list"), "global_cell*")) {
        print(f)
        a <- read.table(paste0(dir, "from_snap/global_ng_list/", f), sep=",", stringsAsFactor=F, header=T)
        cols = colnames(a)[unlist(sapply(colnames(a), function(x){return(any(x == colnames(b)))}))]
        m <- merge(a, b, all.x=TRUE, by.x=c("barcode", "batch", cols), by.y=c("full_barcode", "batch", cols))
        m <- m[,colnames(m) != "X"]
        m <- m[order(m$global_index),]
        m <- cbind(m, global_index_1000=m)
        write.table(m, paste0(dir, "from_snap/global_ng_list/", "global_cell_ng_1_", m[1,"batch"], "_cluster.csv"), sep=",", quote=F)
    }
}