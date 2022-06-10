library("affy")
library("affyPLM")
memory.limit(16000)

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE76124")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE76124.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE50948")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE50948.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE71258")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE71258.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE43365")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE43365.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE140494")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE140494.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE76274")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE76274.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE58792")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE58792.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/GSE43502")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/GSE43502.txt")

setwd("/home/minzhang/workspace/BRCA_meno/cel/all")
raw.set = ReadAffy()
rma.data = rma(raw.set)
expr = exprs(rma.data)
write.table(expr, "/home/minzhang/workspace/BRCA_meno/expr/all.txt")


