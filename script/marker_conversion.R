for (fhead in c("gabaergic", "excitatory")) {
    fname = paste0(fhead, "_markers_fc.txt")
    a <- read.table(fname, header=T, sep=" ")
    for (n in unique(a[,1])) {
        temp <- subset(a, a[,1] == n)
        b <- temp[,2,drop=F]
        write.table(b, paste0(fhead, "_", n, ".txt"), col.names=FALSE, row.names=FALSE, quote=F)
    }
}