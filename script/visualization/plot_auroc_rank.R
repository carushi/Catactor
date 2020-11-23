library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

convert.row.names.to.annot <- function(name_list) {
    find.marker.position <- function(contents) {
        for (i in 1:length(contents)) {
            if (any(contents[i] == c('all', 'SF', 'CU', 'TA', 'TN', 'SC', 'SM')))
                return(i)   
        }
    }
    annot <- NULL
    for (n in name_list) {
        contents <- unlist(strsplit(n, '_'))
        index <- find.marker.position(contents)
        temp <- c(paste(contents[1:(index-1)], collapse='_'), contents[index], paste(contents[(index+1):length(contents)], collapse='_'))
        annot <- rbind(annot, temp)
    }
    return(annot)
}

col = list(Marker=c("all"="darkgrey", "SF"='#E64B35FF', "CU"='#4DBBD5FF', 'TA'='#00A087FF', 'TN'='#91D1C2FF', 'SC'='#3C5488FF'))
# col[['Training']] = c(grey(0.2), grey(0.4), grey(0.6), grey(0.8))
col[['Training']] = rev(brewer.pal(4, 'Reds'))
col[['Supervised']] = rev(brewer.pal(4, 'Blues'))
names(col[['Training']]) = c("ranking", "concensus", "rna_atlas", "rna")
names(col[['Supervised']]) = c('ML', 'Expression', 'Seurat', 'BBKNN')
print(col)

for (fname in list.files("./", pattern="^mat_*")) {
    if (!grepl('csv', fname)) next
    rank=grepl('rank', fname)
    print(fname)
    header = unlist(strsplit(fname, '.', fixed=TRUE))[1]
    a <- read.table(fname, sep=",", header=T)
    method_label <- a[,1]
    a <- a[,-c(1)]
    print(method_label)
    print(head(a))
    for (label in c('all', 'top', 'rna')) {
        if (label == 'top') {
            a <- a[1:15,]
            method_label <- method_label[1:15]
        }
        annot_df <- convert.row.names.to.annot(method_label)
        if (label == 'rna') {
            a <- a[,c(4, 6)]
            a <- cbind(a, "mean AUROC"=rowMeans(a))
            ord <- order(a[,3], decreasing=TRUE)[1:15]
            annot_df <- annot_df[ord,]
            a <- a[ord,]
        }
        rownames(annot_df) <- rownames(a)
        print(head(annot_df))
        colnames(annot_df) <- c('Method', 'Marker', 'Training')
        annot_df <- cbind(annot_df, 'Supervised'=unlist(sapply(annot_df[,'Method'], function(x){
            if (any(x == c('svm', 'logistic', 'la', 'rf')))
                return('ML')
            if (x == 'exp')
                return('Expression')
            if (x == 'seurat_rna')
                return('Seurat')
            return('BBKNN')
        })))
        annot_df <- annot_df[,c("Marker", "Supervised", "Training")]
        annot_df <- data.frame(annot_df)
        annot_df[,'Marker'] <- factor(annot_df[,'Marker'], levels=c('SF', 'CU', 'TA', 'TN', 'SC', 'all'))
        annot_df[,'Supervised'] <- factor(annot_df[,'Supervised'], levels=c('Expression', 'ML', 'Seurat', 'BBKNN'))
        annot_df[,'Training'] <- factor(annot_df[,'Training'], levels=c('ranking', 'concensus', 'rna_atlas', 'rna'))
        print(head(annot_df))
        ha <- HeatmapAnnotation(df=annot_df, col=col, which="row")
        pdf(paste0(header, '_', label, ".pdf"))
        if (label == 'all')
            draw(Heatmap(a, col=viridis(100), column_title="Test data", cluster_rows=FALSE, cluster_columns=FALSE, row_title="Classifier", right_annotation=ha, show_row_names=FALSE))
        else {
            if (rank)
            draw(Heatmap(a, col=viridis(100), column_title="Test data", cluster_rows=FALSE, cluster_columns=FALSE, row_title="Classifier", right_annotation=ha, show_row_names=FALSE, 
                cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%2f", a[i, j]), x, y, gp = gpar(fontsize = 10))}))
            else
            draw(Heatmap(a, col=viridis(100), column_title="Test data", cluster_rows=FALSE, cluster_columns=FALSE, row_title="Classifier", right_annotation=ha, show_row_names=FALSE, 
                cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", a[i, j]), x, y, gp = gpar(fontsize = 10))}))
        }
        dev.off()
    }
}
