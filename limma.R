library(limma)

expr <- read.table("/home/minzhang/workspace/BRCA_meno/expr/all.txt")
groupList <- read.csv("/home/minzhang/workspace/BRCA_meno/clinic/all.csv")
groupList <- groupList$MenoTag

design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c("Control", "Treatment")
fit.lm = lmFit(expr, design)
contrast = makeContrasts("Treatment-Control", levels=design)
fit.cts = contrasts.fit(fit.lm, contrast)
fit.ebs = eBayes(fit.cts)
tempOutput = topTable(fit.ebs, coef=1, n=Inf)
write.csv(tempOutput, "/home/minzhang/workspace/BRCA_meno/limma/limma.csv")

