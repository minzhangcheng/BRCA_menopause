rm(list = ls())
library("DEDS")

memory.limit(16000)

expr = read.table("/home/minzhang/workspace/BRCA_meno/expr/all.txt")
cn = sort(colnames(expr))
expr = as.matrix(expr[, cn])
genes = row.names(expr)

clinic = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/all.csv")
clinic = clinic[order(clinic$GSM), ]
clinic = clinic$MenoTag

deds = deds.stat.linkC(expr, clinic, B=1300)
r = topgenes(deds, number=10000, sort.by="fc", genelist=genes)
write.table(r, "/home/minzhang/workspace/BRCA_meno/DEDS/deds.txt")
