library(SnapATAC)
library(Matrix)
args = commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_dir <- args[2]

# Older version
# for (sp in c("2C", "3C", "4B", "5D")) {
#     x1.sp <- createSnap(file=paste0(sp, ".snap"), sample=sp)
#     x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
#     write.table(data.frame(barcode=paste0(sp, ".", x1.sp@barcode)), paste0(output_dir, "barcodes_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
#     write.table(data.frame(x1.sp@feature), paste0(output_dir, "bin_", sp, "_1000.tsv"), quote=F, col.names=F, row.names=F)
#     write.table(x1.sp@metaData, paste0(output_dir, "qc_", sp, "_1000.tsv"), quote=F)
#     writeMM(x1.sp@bmat, file=paste0(output_dir, "sparse_mat_", sp, "_1000.mtx"))
# }

count <- 1
global_samples <- 0
for (f in list.files(input_dir, "*\.snap")) {
    print(file.path(input_dir, f))
    region = unlist(strsplit(unlist(strsplit(f, "_"))[2], ".", fixed=TRUE))[1]
    sp = paste0(region, count)
    x1.sp <- createSnap(file=file.path(input_dir, f), sample=sp)
    x1.sp <- addBmatToSnap(obj=x1.sp, bin.size=1000)
    write.table(data.frame(barcode=paste0(sp, ".", x1.sp@barcode)), file.path(output_dir, paste0("barcodes_", sp, "_1000.tsv")), quote=F, col.names=F, row.names=F)
    write.table(data.frame(x1.sp@feature), file.path(output_dir, paste0("bin_", sp, "_1000.tsv")), quote=F, col.names=F, row.names=F)
    data = x1.sp@metaData
    data[,"local_index"] = 1:dim(data)[1]
    data[,"global_index"] = 1:dim(data)[1]+global_samples
    data[,"batch"] = sp
    write.table(data, file.path(output_dir, paste0("qc_", sp, "_1000.tsv")), quote=F)
    writeMM(x1.sp@bmat, file=file.path(output_dir, paste0("sparse_mat_", sp, "_1000.mtx")))
    count <- count+1
    global_samples <- global_samples+dim(data)[1]
}