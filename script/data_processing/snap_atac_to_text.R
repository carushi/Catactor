library(SnapATAC)
library(Matrix)
dir <- ""
output_dir <- "/data/rkawaguc/data/191210_new_BICCN/from_snap/"

for (sp in c("2C", "3C", "4B", "5D")[2:3]) {
    x1.sp <- createSnap(file=paste0(sp, ".snap"), sample=sp)
    x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    write.table(data.frame(barcode=paste0(sp, ".", x1.sp@barcode)), paste0(output_dir, "barcodes_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
    write.table(data.frame(x1.sp@feature), paste0(output_dir, "bin_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
    write.table(x1.sp@metaData, paste0(output_dir, "qc_", sp, "_1000.tsv"), quote=F)
    writeMM(x1.sp@bmat, file=paste0(output_dir, "sparse_mat_", sp, "_1000.mtx"))
}
count <- 1
global_samples <- 0
for (f in list.files("./", "C*.snap")) {
    print(f)
    region = unlist(strsplit(unlist(strsplit(f, "_"))[2], ".", fixed=TRUE))[1]
    sp = paste0(region, count)
    x1.sp <- createSnap(file=f, sample=sp)
    x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    write.table(data.frame(barcode=paste0(sp, ".", x1.sp@barcode)), paste0(output_dir, "barcodes_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
    write.table(data.frame(x1.sp@feature), paste0(output_dir, "bin_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
    data = x1.sp@metaData
    data[,"local_index"] = 1:dim(data)[1]
    data[,"global_index"] = 1:dim(data)[1]+global_samples
    data[,"batch"] = sp
    write.table(data, paste0(output_dir, "qc_", sp, "_1000.tsv"), quote=F)
    writeMM(x1.sp@bmat, file=paste0(output_dir, "sparse_mat_", sp, "_1000.mtx"))
    count <- count+1
    global_samples <- global_samples+dim(data)[1]
}