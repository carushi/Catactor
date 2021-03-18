# Get the expression matrix

# Sys.setenv(RETICULATE_PYTHON = "/home/rkawaguc/anaconda3/envs/BICCN/bin/python")

library("pROC")
library("reticulate")
py_config()

suppressPackageStartupMessages({
    library("ggplot2")
    library("SingleCellExperiment")
    library("Seurat")
})
source_python("read_pickle_data.py")

library(future)
# check the current active plan
plan("multiprocess", workers = 5)
plan()
MAX_CELL <<- 10000
options(future.globals.maxSize = 2500 * 1024^2 * 10)


split.test.dataset <- function(celltypes, a) {
    print('split')
    require(caret)
    set.seed(42)
    index <- which(unlist(sapply(celltypes, function(x){return(!'NA' %in% x)})))
    celltypes <- factor(celltypes[index])
    max_fold <- max(1, floor(length(celltypes)/MAX_CELL))
    print(celltypes)
    print(MAX_CELL)
    print(max_fold)
    fname = paste0('test_split_', a, '_', max_fold, '.rds')
    if (file.exists(fname)) {
        folds <- readRDS(fname)
    } else {
        folds <- createFolds(celltypes, k=max_fold)
        saveRDS(folds, fname)
    }
    return(list(max_fold, folds, index))
}

filter.train.dataset <- function(celltypes, b) {
    print('filter')
    require(caret)
    set.seed(42)
    index <- which(unlist(sapply(celltypes, function(x){return(!'NA' %in% x)})))
    celltypes <- factor(celltypes[index])
    fname = paste0('train_split_', b, '_', MAX_CELL, '.rds')
    p = 10000/length(celltypes)
    if (file.exists(fname)) {
        part <- readRDS(fname)
    } else {
        if (length(index) <= MAX_CELL) {
            part <- 1:length(index)
        } else {
            print(c(length(celltypes), p))
            part <- createDataPartition(celltypes, times = 1, p = p, list = TRUE)
        }
        saveRDS(part, fname)
    }
    return(list(part, index))
}
make.filtered.data <- function(data, index, index_of_index) {
    temp <- list()
    temp$X <- data$X[index[index_of_index],]
    temp$var <- data$var
    temp$obs <- data$obs[index[index_of_index],]
    return(temp)
}

run.for.all.combinations.wo.merge <- function(index=-1) 
{
    for (m in 1:2) {
        if (m == 1)
            fn <- extract_atac_rna_combination_wo_merge()
        else
            fn <- extract_atac_atac_combination_wo_merge()
        i = 1
        while (TRUE) { 
            results <- iter_next(fn)
            if (is.null(results[[1]])) break
            if (index > 0 && index != i) {
                i = i+1
                next
            }
            header <- results[[1]]
            data_a <- results[[2]]
            data_b <- results[[3]]
            new_data_a <- NULL
            new_data_b <- NULL
            a <- results[[4]]
            b <- results[[5]]
            print(c(a, b))
            if (m == 2 && dim(data_b$X)[1] > 10000) { # atac downsamlpnig
                print(data_b$obs)
                print(dim(data_b$X))
                print(dim(data_b$obs))
                cell_index <- filter.train.dataset(data_b$obs$celltype, b) 
                part <- cell_index[[1]]
                index <- cell_index[[2]]
                new_data_b <- make.filtered.data(data_b, index, unlist(part))
            }
            if (dim(data_a$X)[1] >  10000) {
                cell_index <- split.test.dataset(data_a$obs$celltype, a)
                max_index <- cell_index[[1]]
                folds <- cell_index[[2]]
                index <- cell_index[[3]]
                for (j in 1:max_index) {
                    new_data_a <- make.filtered.data(data_a, index, folds[[j]])
                    if (is.null(new_data_b))
                        sctransform.integrate.split(paste0(header, '_', j), new_data_a, data_b)
                    else
                        sctransform.integrate.split(paste0(header, '_', j), new_data_a, new_data_b)
                }
            } else {
                if (is.null(new_data_b))
                    sctransform.integrate.split(paste0(header, '_', 1), data_a, data_b)
                else
                    sctransform.integrate.split(paste0(header, '_', 1), data_a, new_data_b)
            }
            i = i+1
        }
    }
}

down.sample <- function(obj, size) {
    set.seed(57)
    print(c(colnames(obj), size))
    downsampled.obj <- obj[, sample(colnames(obj), size = size, replace=F)]    
    return(downsampled.obj)
}
read.marker.genes <- function() {
    a <- read.table('../data/marker_name_list.csv', header=T, sep=",", row.names=1)
    print(head(a))
    marker_list <- c()
    for (name in colnames(a)) {
        marker_list <- c(marker_list, unlist(strsplit(name, '_'))[1])
    }
    marker_list <- marker_list[!marker_list %in% c('SM')]
    marker_list <- c(marker_list, colnames(a)[grep('SM_', colnames(a))])
    marker_gene_list <- list()
    for (m in unique(marker_list)) {
        vec <- unlist(a[,grepl(paste0(m), colnames(a)),drop=F][1:100,])
        vec <- as.vector(vec[!is.na(vec)])
        vec <- vec[vec != '']
        marker_gene_list[[m]] <- vec
    }
    for (m in names(marker_gene_list)) {
        print(c(m, length(marker_gene_list[[m]])))
    }
    return(marker_gene_list)
}
make.seurat.count.matrix <- function(adata) {
    print('make')
    exprs <- t(as.matrix(adata$X))
    if (!is.null(adata$obs_names)) {
        colnames(exprs) <- adata$obs_names$to_list()
    } else {
        colnames(exprs) <- rownames(adata$obs)
    }
    if (!is.null(adata$var_names)) {
        rownames(exprs) <- adata$var_names$to_list()
    } else {
        rownames(exprs) <- rownames(adata$var)
    }
    seurat <- CreateSeuratObject(exprs)
    seurat <- AddMetaData(seurat, adata$obs)
    rm(exprs)
    return(seurat)
}
sctransform.integrate.split <- function(output, data_a, data_b, method='cca') {
    marker <- read.marker.genes()
    marker <- c(marker, all=c(''))
    print(marker)
    print(c('set exprs', output))
    seurat.list <- list()
    seurat.list[[1]] <- make.seurat.count.matrix(data_a)
    remove(data_a)
    seurat.list[[2]] <- make.seurat.count.matrix(data_b)
    remove(data_b)
    gc()
    for (i in 1:2) {
        seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]])
        seurat.list[[i]] <- NormalizeData(seurat.list[[i]])
        seurat.list[[i]] <- ScaleData(seurat.list[[i]])
        # seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
    }
    require(Signac)
    print('cca')
    gc()
    query <- seurat.list[[1]]
    reference <- seurat.list[[2]]
    remove(seurat.list)
    gc()
    for (m in names(marker)) {
        print(c(m, length(marker[[m]])))
        if (m == 'all') {
            reference <- RunPCA(reference, verbose = FALSE)
            query <- RunPCA(query, verbose = FALSE)
        } else {
            reference <- RunPCA(reference, pc.genes = marker[[m]], verbose = FALSE)
            query <- RunPCA(query, pc.genes = marker[[m]], verbose = FALSE)
        }
        print('find transfer anchors')
        anchors <- FindTransferAnchors(reference = reference, query = query, features = VariableFeatures(object = reference), reference.assay = "RNA", query.assay = "RNA", reduction='cca')
        print('transfer')
        print(anchors)
        if (m == 'all') {
            dims = 2:30
        } else {
            dims = 2:min(30, length(marker[[m]]))
        }
        predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = dims, weight.reduction='cca')
        seurat.query <- AddMetaData(object = query, metadata = predictions)
        mat <- cbind(predictions, answer=seurat.query$celltype)
        # compute.auroc(seurat.query$celltype, predictions)
        mat <- cbind(mat, cbind(predicted=seurat.query$predicted.id, celltype=seurat.query$celltype))
        # compute.auroc(mat, paste0('SCT_', output, '.csv'))
        mat <- cbind(mat, flag=(mat[,1] == mat[,2]))
        write.table(mat, paste0('SCT_', output, '_', m, '.csv'))
    }
    gc()
}


args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
    run.for.all.combinations.wo.merge(args[1])
} else {
    run.for.all.combinations.wo.merge()
}
