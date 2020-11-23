require(ggplot2)
require(viridis)
require(UpSetR)
require(ComplexHeatmap)
require(plyr)
df <- read.table("marker_name_list.csv", header=T, sep=",")
df <- df[,-c(1)]
df <- df[1:100,]
print(head(df))
gdf <- read.table("rank_gene_list_celltype.csv", header=T, sep=",")
gdf <- gdf[,-c(1)]
gdf <- gdf[1:1000,]
print(head(gdf))
for (marker_type in c('major', 'minor', 'data')) {
    if (any(marker_type == c('major', 'data'))) {
        tdf <- df[,!grepl('SM', colnames(df), fixed=TRUE)]
        if (marker_type == 'data') {
            # mdf <- t(rbind.fill(as.data.frame(t(tdf)), as.data.frame(t(gdf))))
            # colnames(mdf) <- c(colnames(tdf), colnames(gdf))
            tdf <- gdf
            print(head(tdf))
            # exit()
        }
    } else {
        tdf <- df[,grepl('SM', colnames(df), fixed=TRUE)]
    }
    if (any(marker_type == c('major', 'data', 'minor'))) {
        for (celltype in c('NN', 'IN', 'EX')) {
            temp <- as.list(as.data.frame(tdf))
            print(grepl(celltype, names(temp)))
            print(which(grepl(celltype, names(temp))))
            temp <- temp[which(grepl(celltype, names(temp)))]
            temp <- lapply(temp, function(x) {
                x <- sapply(x, function(x){if(!is.na(x) && x != '') return(as.character(x))})
                return(unlist(x))
            })
            all_factors <- unique(unlist(temp))
            temp <- lapply(temp, function(x) {
                return(as.numeric(factor(x, levels=all_factors)))
            })
            print(temp)
            # tdf <- list_to_matrix(tdf)
            # png(paste0('upset_', marker_type, '.png'), width=1000, height=1000)
            # pdf(paste0('upset_', marker_type, '.pdf'), width=10, height=7, onefile=FALSE)
            upset(fromList(temp), nsets=length(temp), order.by="freq")
            # ggsave(paste0('upset_', marker_type, '_', celltype, '.pdf'))
            # dev.off()
            # break
            # exit()

            # labels = venn.get_labels([range(10), range(5, 15), range(3, 8), range(8, 17), range(10, 20), range(13, 25)], fill=['number', 'logic'])
            # fig, ax = venn.venn6(labels, names=['list 1', 'list 2', 'list 3', 'list 4', 'list 5', 'list 6'])
            # fig.show()
        }
    }
}