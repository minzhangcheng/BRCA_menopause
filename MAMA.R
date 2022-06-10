rm(list = ls())
library("MAMA")
library("metaMA")
library("affyPLM")
library("affy")

memory.limit(16000)

GSE76124 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE76124.txt")
cn = sort(colnames(GSE76124))
GSE76124 = as.matrix(GSE76124[, cn])
GSE50948 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE50948.txt")
cn = sort(colnames(GSE50948))
GSE50948 = as.matrix(GSE50948[, cn])
GSE71258 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE71258.txt")
cn = sort(colnames(GSE71258))
GSE71258 = as.matrix(GSE71258[, cn])
GSE43365 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE43365.txt")
cn = sort(colnames(GSE43365))
GSE43365 = as.matrix(GSE43365[, cn])
GSE140494 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE140494.txt")
cn = sort(colnames(GSE140494))
GSE140494 = as.matrix(GSE140494[, cn])
GSE76274 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE76274.txt")
cn = sort(colnames(GSE76274))
GSE76274 = as.matrix(GSE76274[, cn])
GSE58792 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE58792.txt")
cn = sort(colnames(GSE58792))
GSE58792 = as.matrix(GSE58792[, cn])
GSE43502 = read.table("/home/minzhang/workspace/BRCA_meno/expr/GSE43502.txt")
cn = sort(colnames(GSE43502))
GSE43502 = as.matrix(GSE43502[, cn])
GEDM=list(GSE76124, GSE50948, GSE71258, GSE43365, GSE140494, GSE76274, GSE58792, GSE43502)
rm(GSE76124, GSE50948, GSE71258, GSE43365, GSE140494, GSE76274, GSE58792, GSE43502)
gc()

setwd("/home/minzhang/workspace/BRCA_meno/MAMA")
GSE76124 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE76124.csv")
GSE76124 = GSE76124[order(GSE76124$GSM), ]
row.names(GSE76124) = GSE76124$GSM
GSE76124 = GSE76124["MenoTag"]
GSE76124["MenoTag"] = as.factor(GSE76124$MenoTag)
GSE50948 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE50948.csv")
GSE50948 = GSE50948[order(GSE50948$GSM), ]
row.names(GSE50948) = GSE50948$GSM
GSE50948 = GSE50948["MenoTag"]
GSE50948["MenoTag"] = as.factor(GSE50948$MenoTag)
GSE71258 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE71258.csv")
GSE71258 = GSE71258[order(GSE71258$GSM), ]
row.names(GSE71258) = GSE71258$GSM
GSE71258 = GSE71258["MenoTag"]
GSE71258["MenoTag"] = as.factor(GSE71258$MenoTag)
GSE43365 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE43365.csv")
GSE43365 = GSE43365[order(GSE43365$GSM), ]
row.names(GSE43365) = GSE43365$GSM
GSE43365 = GSE43365["MenoTag"]
GSE43365["MenoTag"] = as.factor(GSE43365$MenoTag)
GSE140494 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE140494.csv")
GSE140494 = GSE140494[order(GSE140494$GSM), ]
row.names(GSE140494) = GSE140494$GSM
GSE140494 = GSE140494["MenoTag"]
GSE140494["MenoTag"] = as.factor(GSE140494$MenoTag)
GSE76274 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE76274.csv")
GSE76274 = GSE76274[order(GSE76274$GSM), ]
row.names(GSE76274) = GSE76274$GSM
GSE76274 = GSE76274["MenoTag"]
GSE76274["MenoTag"] = as.factor(GSE76274$MenoTag)
GSE58792 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE58792.csv")
GSE58792 = GSE58792[order(GSE58792$GSM), ]
row.names(GSE58792) = GSE58792$GSM
GSE58792 = GSE58792["MenoTag"]
GSE58792["MenoTag"] = as.factor(GSE58792$MenoTag)
GSE43502 = read.csv("/home/minzhang/workspace/BRCA_meno/clinic/GSE43502.csv")
GSE43502 = GSE43502[order(GSE43502$GSM), ]
row.names(GSE43502) = GSE43502$GSM
GSE43502 = GSE43502["MenoTag"]
GSE43502["MenoTag"] = as.factor(GSE43502$MenoTag)
Clinical = list(GSE76124, GSE50948, GSE71258, GSE43365, GSE140494, GSE76274, GSE58792, GSE43502)
rm(GSE76124, GSE50948, GSE71258, GSE43365, GSE140494, GSE76274, GSE58792, GSE43502)
gc()

datanames = c("GSE76124", "GSE50948", "GSE71258", "GSE43365", "GSE140494", "GSE76274", "GSE58792", "GSE43502")

setClass("esets",slots=list(GEDM="list",clinical="list",datanames="character"),package = "MAMA")
meta.data = new("esets", GEDM=GEDM, clinical=Clinical, datanames=datanames)
rm(Clinical, datanames)
gc()

pval = metaMA(meta.data,"MenoTag",which="pval")
gc()
es2 = ES.GeneMeta(meta.data,"MenoTag",nperm=1000)

results1 = join.results(pval, type=1, genenames=rownames(GEDM(meta.data)[[1]]))
p_value = as.data.frame(results1)
rawpval = 2 * (1 - pnorm(abs(pval$TestStatistic)))
FDR_pval = p.adjust(rawpval, method="BY", n=length(rawpval))
p_value$c_pval = rawpval
p_value$FDR = FDR_pval
rm(rawpval, FDR_pval)
gc()
es2_theScores = es2$theScores
es2_ScoresFDR = es2$ScoresFDR
es2_ScoresFDR = es2_ScoresFDR$two.sided
write.table(p_value, file="/home/minzhang/workspace/BRCA_meno/MAMA/p_value", sep="\t", col.names = T)
write.table(es2_ScoresFDR, file="/home/minzhang/workspace/BRCA_meno/MAMA/es", sep="\t", col.names = T)

