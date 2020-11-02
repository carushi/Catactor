# Get the expression matrix

Sys.setenv(RETICULATE_PYTHON = "/home/rkawaguc/anaconda3/envs/BICCN/bin/python")
# install.packages('reticulate')
library("pROC")
library("reticulate")
py_config()
# memory.limit(size=100000)

suppressPackageStartupMessages({
    library("ggplot2")
    library("SingleCellExperiment")
    # library("scater")
    library("Seurat")
})
source_python("read_pickle_data.py")

library(future)
# check the current active plan
plan("multiprocess", workers = 5)
plan()
MAX_CELL <<- 10000
options(future.globals.maxSize = 2500 * 1024^2 * 10)

# celltype.conversion <- function(problem, vec) {
#     if (problem == 'celltype') {

#     } else if (problem == 'neuron') {

#     } else if (problem == 'inex') {

#     }
#     return(new_vec)
#     if '_' in X[0]:
#         celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
#     else:
#         celltype_without_num = list(X)
#     print(set(X))
#     celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'OT' for x in celltype_without_num]
#     print(set(celltype_without_var))
#     return celltype_without_var
# }

run.for.all.combinations <- function(index=-1) 
{
    for (m in 2:2) {
        if (m == 1)
            fn <- extract_atac_rna_combination()
        else
            fn <- extract_atac_atac_combination()
        i = 1
        while (TRUE) { 
            print(i)
            results <- iter_next(fn)
            if (is.null(results[[1]])) break
            if (index > 0 && index != i) {
                i = i+1
                next
            }
            header <- results[[1]]
            data <- results[[2]]
            sctransform.integrate(header, data, filter=FALSE)
            i = i+1
        }
    }
}

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
    for (m in 2:2) {
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
            # i = i+1
            # next
            if (m == 2 && dim(data_b$X)[1] > 10000) {
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
    # return(list(all=c('')))
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
    }
    require(Signac)
    print('cca')
    # seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
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

sctransform.integrate <- function(output, adata, method='cca', filter=TRUE) {
    marker <- read.marker.genes()
    marker <- c(marker, all=c())
    exprs <- t(adata$X)
    colnames(exprs) <- adata$obs_names$to_list()
    rownames(exprs) <- adata$var_names$to_list()
    print('set exprs')
    metadata <- adata$obs
    remove(adata)
    gc()
    print('create seurat')
    # Create the Seurat object
    seurat <- CreateSeuratObject(exprs)
    print('made')
    # seurat <- SetAssayData(seurat, "data", exprs)
    remove(exprs)
    seurat <- AddMetaData(seurat, metadata)

    # seurat[["RNA"]][["n_cells"]] <- adata$var["n_cells"]
    print('split')
    seurat.list <- SplitObject(seurat, split.by = "batch")
    seurat.list[[2]] 
    remove(seurat)
    for (i in names(seurat.list)) {
        print(i)
        size = 1000
        if (ncol(seurat.list[[i]]) >= size) {
            print('downsample')
            print(length(colnames(seurat.list[[i]])))
            seurat.list[[i]] <- down.sample(seurat.list[[i]], size)
        }
        if (method == 'cca') {
            seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]])
            seurat.list[[i]] <- NormalizeData(seurat.list[[i]])
            seurat.list[[i]] <- ScaleData(seurat.list[[i]])
        } else {
            seurat.list[[i]] <- SCTransform(seurat.list[[i]], verbose = FALSE)
        }
    }
    gc()
    if (method == "sct") {
        print('sct')
        seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
        seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = seurat.features, 
            verbose = FALSE)
        seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
            anchor.features = seurat.features, verbose = FALSE)
        # seurat.anchors <- FindTransferAnchors(object.list = seurat.list, normalization.method = "SCT", 
        #     anchor.features = seurat.features, verbose = FALSE)
        print('anchors')
        gc()

        seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", 
            verbose = FALSE)
        print('pca')
        seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
        print('umap')
        seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30)
        print('integrated')
        pdf(paste0('SCT_', output, '.pdf'), width=10, height=7)
        plots <- DimPlot(seurat.integrated, group.by = c("batch", "celltype"))
        # plots <- plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
        #     override.aes = list(size = 3)))
        plot(plots)
        dev.off()
        seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
        seurat.integrated <- FindClusters(seurat.integrated, resolution=0.8)
        remove(seurat.integrated)
        gc()
        query <- seurat.list[[1]]
        reference <- seurat.list[[2]]
        rm(seurat.list)
        gc()
        reference <- RunPCA(reference, verbose = FALSE)
        query <- RunPCA(query, verbose = FALSE)
        anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30)
        print(anchors)
        predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = 1:30, k.weight=5)
        print('query')
        seurat.query <- AddMetaData(object = query, metadata = predictions)
        seurat.query$prediction.match <- seurat.query$predicted.id == seurat.query$celltype
        print(table(seurat.query$prediction.match))
        print(table(seurat.query$predicted.id))
    # }
    # commongenes_tokeep <- Reduce(intersect, lapply(SCT_Integrated.anchors@object.list, rownames))
    } else if (method == "cca") {
        require(Signac)
        print('cca')
        # seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
        gc()
        query <- seurat.list[[1]]
        reference <- seurat.list[[2]]
        rm(seurat.list)
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
            if (m == 'all')
                dims = 2:30
            else
                dims = 2:min(30, length(marker[[m]]))
            predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = dims, weight.reduction='cca')
            seurat.query <- AddMetaData(object = query, metadata = predictions)
            mat <- cbind(predictions, answer=seurat.query$celltype)
            # compute.auroc(seurat.query$celltype, predictions)
            mat <- cbind(seurat.query$predicted.id, celltype=seurat.query$celltype)
            # compute.auroc(mat, paste0('SCT_', output, '.csv'))
            mat <- cbind(mat, flag=(mat[,1] == mat[,2]))
            write.table(mat, paste0('SCT_', output, '_', m, '.csv'))
            # seurat.query$prediction.match <- seurat.query$predicted.id == seurat.query$celltype
            # print(table(seurat.query$prediction.match))
            # print(table(seurat.query$predicted.id))
        }
    }
    gc()
}
# compute.auroc <- function(mat) {
#     nmat <- mat[mat[,'answer'] != 'NA']
#     auc_vec <- NULL
#     for (col in colnames(mat)) {
#         if (grepl('prediction.score', col)) {
#             celltype = unlist(gsub("prediction.score.", '', col))[1]
#             if (celltype == 'NA') next
#             answer <- unlist(sapply(mat[,'answer'], function(x){if(x == celltype)return(1); return(0)}))
#             score <- mat[,col]
#             rocobj <- roc(answer, score)
#             auc_val <- rbind(auc_vec, c('celltype', celltype, auc(rocobj)))
#             if (any(celltype = c('IN', 'EX'))) {

#             }
#         }
#     }
# }
test.load <- function(fname='', method='SCT', output='output') {
    if (fname == '') 
        fname = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/bbknn_object/bbknn_GSE123576_BICCN2_atac_rna_10_TN.pyn"
    adata <- read_pickle_file(fname)
    sce <- SingleCellExperiment(
        assays      = list(logcounts = t(adata$X)),
        colData     = adata$obs,
        rowData     = adata$var,
        # reducedDims = list(umap = adata$obsm["X_umap"])
    )
    exprs <- t(as.matrix(adata$X))
    colnames(exprs) <- adata$obs_names$to_list()
    rownames(exprs) <- adata$var_names$to_list()
    # Create the Seurat object
    seurat <- CreateSeuratObject(exprs)
    seurat <- SetAssayData(seurat, "data", exprs)
    seurat <- AddMetaData(seurat, adata$obs)
    # seurat[["RNA"]][["n_cells"]] <- adata$var["n_cells"]
    # Add embedding
    embedding <- adata$obsm["X_umap"]
    rownames(embedding) <- adata$obs_names$to_list()
    colnames(embedding) <- c("umap_1", "umap_2")
    seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")  
    seurat.list <- SplitObject(seurat, split.by = "batch")
    for (i in names(seurat.list)) {
        seurat.list[[i]] <- SCTransform(seurat.list[[i]], verbose = FALSE)
    }
    # if (method == "SCT") {
    #     seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
    #     seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = seurat.features, 
    #         verbose = FALSE)
    #     seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
    #         anchor.features = seurat.features, verbose = FALSE)
    #     seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", 
    #         verbose = FALSE)
    #     seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
    #     seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30)
    #     plots <- DimPlot(seurat.integrated, group.by = c("tech", "celltype"))
    #     png(paste0('SCT_', output, '.png'))
    #     plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    #         override.aes = list(size = 3)))
    #     dev.off()
    #     seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
    #     seurat.integrated <- FindClusters(seurat.integrated, resolution=0.8)

    #     query <- seurat.list[[1]]
    #     reference <- seurat.list[[2]]
    #     anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30)
    #     predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = 1:30)
    #     seurat.query <- AddMetaData(object = query, metadata = predictions)
    #     seurat.query$prediction.match <- seurat.query$predicted.id == seurat.query$celltype
    #     print(table(seurat.query$prediction.match))
    #     print(table(seurat.query$predicted.id))
    # }
    
    
    # commongenes_tokeep <- Reduce(intersect, lapply(SCT_Integrated.anchors@object.list, rownames))
    # SCT_integrated <- IntegrateData(anchorset = SCT_Integrated.anchors, normalization.method = "SCT", features.to.integrate = commongenes_tokeep)
    # SCT_integrated <- RunPCA(SCT_integrated)
    # SCT_integrated <- RunUMAP(SCT_integrated, dims = 1:30)
    # SCT_integrated <- FindNeighbors(SCT_integrated, dims = 1:30)
    # SCT_integrated <- FindClusters(SCT_integrated, resolution=0.8)

    # DefaultAssay(Sample_subset) <- "RNA"
    # Sample_subset <- subset(SCT_integrated, idents = c(0, 2, 3, 7))
    # Sample_subset <- ScaleData(Sample_subset, vars.to.regress = c("percent.mt")
    # Sample_subset <- RunPCA(Sample_subset)
    # Sample_subset <- RunUMAP(Sample_subset, dims = 1:30)
    # Sample_subset <- FindNeighbors(Sample_subset, dims = 1:30)
    # Sample_subset <- FindClusters(Sample_subset, resolution = 0.8)



    # } else if (method == "CCA") {

    # transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
    #     reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
    # celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype, 
    #     weight.reduction = pbmc.atac[["lsi"]])
    # pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
    # hist(pbmc.atac$prediction.score.max)
    # abline(v = 0.5, col = "red")
    # table(pbmc.atac$prediction.score.max > 0.5)
    # }
}

# read.rna.data <- function() {
    
# }


# pancreas <- CreateSeuratObject(pancreas.data, meta.data = metadata)
# pancreas.list <- SplitObject(pancreas, split.by = "tech")

# pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
# pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
#     verbose = FALSE)

#     pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
#     anchor.features = pancreas.features, verbose = FALSE)
# pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
#     verbose = FALSE)
# pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
# pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
# plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"), combine = FALSE)
# plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
#     byrow = TRUE, override.aes = list(size = 3))))
# CombinePlots(plots)    

# integrate.same.feature.dim <- function(sce) {
# # pbmc.data <- readRDS("../data/pbmc_ssc_mat.rds")
# # pbmc.metadata <- readRDS("../data/pbmc_ssc_metadata.rds")
# # pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data = pbmc.metadata)
# # pbmc <- subset(pbmc, subset = nFeature_RNA > 200)
# pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
# pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
# pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
#     anchor.features = pbmc.features)
# pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")

# pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
# pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)
# plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"), combine = FALSE)
# plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
#     byrow = TRUE, override.aes = list(size = 3))))
# CombinePlots(plots)
# DimPlot(pbmc.integrated, group.by = "CellType", split.by = "Method", ncol = 3)
# DefaultAssay(pbmc.integrated) <- "RNA"
# # Normalize RNA data for visualization purposes
# pbmc.integrated <- NormalizeData(pbmc.integrated, verbose = FALSE)
# FeaturePlot(pbmc.integrated, c("CCR7", "S100A4", "GZMB", "GZMK", "GZMH", "TCL1A"))
# # exprs <- t(py$adata$X)
# # colnames(exprs) <- py$adata$obs_names$to_list()
# # rownames(exprs) <- py$adata$var_names$to_list()
# # # Create the Seurat object
# # seurat <- CreateSeuratObject(exprs)
# # # Set the expression assay
# # seurat <- SetAssayData(seurat, "data", exprs)
# # # Add observation metadata
# # seurat <- AddMetaData(seurat, py$adata$obs)
# # # Add fetaure metadata
# # seurat[["RNA"]][["n_cells"]] <- py$adata$var["n_cells"]
# # # Add embedding
# # embedding <- py$adata$obsm["X_umap"]
# # rownames(embedding) <- py$adata$obs_names$to_list()
# # colnames(embedding) <- c("umap_1", "umap_2")
# # seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")


# # pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
# # pbmc <- RunSVD(pbmc)
# # pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
# # pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
# # pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
# # DimPlot(object = pbmc, label = TRUE) + NoLegend()
# # gene.activities <- GeneActivity(pbmc)
# # add the gene activity matrix to the Seurat object as a new assay and normalize it

# counts <- Read10X_h5(filename = "/home/stuartt/github/chrom/vignette_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
# metadata <- read.csv(
#   file = "/home/stuartt/github/chrom/vignette_data/atac_v1_pbmc_10k_singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )

# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":", "-"),
#   genome = 'hg19',
#   fragments = '/home/stuartt/github/chrom/vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
#   min.cells = 10,
#   min.features = 200
# )

# pbmc <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks",
#   meta.data = metadata
# )
# peakdata 
# pbmc <- readPeak()
# pbmc <- RunTFIDF(pbmc)

# pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
# pbmc <- NormalizeData(
#   object = pbmc,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(pbmc$nCount_RNA)
# )
# DefaultAssay(pbmc) <- 'RNA'

# FeaturePlot(
#   object = pbmc,
#   features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )
# # Load the pre-processed scRNA-seq data for PBMCs
# pbmc_rna <- readRDS("/home/stuartt/github/chrom/vignette_data/pbmc_10k_v3.rds")

# transfer.anchors <- FindTransferAnchors(
#   reference = pbmc_rna,
#   query = pbmc,
#   reduction = 'cca'
# )

# predicted.labels <- TransferData(
#   anchorset = transfer.anchors,
#   refdata = pbmc_rna$celltype,
#   weight.reduction = pbmc[['lsi']],
#   dims = 2:30
# )

# pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
# plot1 <- DimPlot(
#   object = pbmc_rna,
#   group.by = 'celltype',
#   label = TRUE,
#   repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

# plot2 <- DimPlot(
#   object = pbmc,
#   group.by = 'predicted.id',
#   label = TRUE,
#   repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

# plot1 + plot2
# pbmc <- subset(pbmc, idents = 14, invert = TRUE)
# pbmc <- RenameIdents(
#   object = pbmc,
#   '0' = 'CD14 Mono',
#   '1' = 'CD4 Memory',
#   '2' = 'CD8 Effector',
#   '3' = 'CD4 Naive',
#   '4' = 'CD14 Mono',
#   '5' = 'DN T',
#   '6' = 'CD8 Naive',
#   '7' = 'NK CD56Dim',
#   '8' = 'pre-B',
#   '9' = 'CD16 Mono',
#   '10' = 'pro-B',
#   '11' = 'DC',
#   '12' = 'NK CD56bright',
#   '13' = 'pDC'
# )
# DefaultAssay(pbmc) <- 'peaks'

# da_peaks <- FindMarkers(
#   object = pbmc,
#   ident.1 = "CD4 Naive",
#   ident.2 = "CD14 Mono",
#   min.pct = 0.2,
#   test.use = 'LR',
#   latent.vars = 'peak_region_fragments'
# )

# head(da_peaks)

# plot1 <- VlnPlot(
#   object = pbmc,
#   features = rownames(da_peaks)[1],
#   pt.size = 0.1,
#   idents = c("CD4 Memory","CD14 Mono")
# )
# plot2 <- FeaturePlot(
#   object = pbmc,
#   features = rownames(da_peaks)[1],
#   pt.size = 0.1
# )

# plot1 | plot2

# fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14 Mono")
# head(fc)
# open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
# open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

# closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
# closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)


args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
#     fname = args[1]
#     data.load(fname)
# } else if (length(args) > 1) {
#     fname = args[1]
#     peak_fname = args[2]
#     data.peak.load(fname, peak_fname)
# } else  {
    run.for.all.combinations.wo.merge(args[1])
    # test.load()
} else {
    run.for.all.combinations.wo.merge()
}
