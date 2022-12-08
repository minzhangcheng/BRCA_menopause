import tarfile
import os
import shutil
import pandas as pd
from GEO.utility import mkdir
from GEO.utility import extractGzipFile


rootDir = '/home/minzhang/workspace/BRCA_meno'
clincDir = '{}/clinic'.format(rootDir)
rawDir = '{}/raw'.format(rootDir)
celDir = '{}/cel'.format(rootDir)
clinicDir = '{}/clinic'.format(rootDir)
exprDir = '{}/expr'.format(rootDir)
mamaDir = '{}/MAMA'.format(rootDir)
dedsDir = '{}/DEDS'.format(rootDir)
limmaDir = '{}/limma'.format(rootDir)
tag = 'MenoTag'

"""
gseList = ['GSE76124', 'GSE50948', 'GSE71258', 'GSE23177', 'GSE43365',
           'GSE140494', 'GSE76274', 'GSE58792', 'GSE22035', 'GSE43502']
"""
gseList = ['GSE76124', 'GSE50948', 'GSE71258', 'GSE43365', 'GSE140494',
           'GSE76274', 'GSE58792', 'GSE43502']
gseListMama = ['GSE76124', 'GSE50948', 'GSE71258', 'GSE43365', 'GSE140494',
           'GSE76274', 'GSE58792', 'GSE43502']

"""
mkdir(rawDir)
getGseRawList(gseList, rawDir, log=True, saveFilenameList=True, 
              maxTrial=4, thread=2, downloadMethod='wget')
"""


# prepare cel files and generate R script to caculate expression matrix
mkdir(celDir)
mkdir(exprDir)

exprR = 'library("affy")\n' \
          'library("affyPLM")\n' \
          'memory.limit(16000)\n\n'


for gse in gseList:
    clinic = pd.read_table('{}/{}.txt'.format(clinicDir, gse), sep='\t',
                           index_col=0)
    clinic = clinic[clinic['MenoTag'] < 2]['MenoTag']
    clinic = clinic.astype(int)
    clinic.to_csv('{}/{}.csv'.format(clinicDir, gse))
    gsmList = list(clinic.index)

    d = '{}/{}'.format(celDir, gse)
    mkdir(d)
    mkdir('{}/all'.format(celDir))
    with tarfile.open('{}/{}_RAW.tar'.format(rawDir, gse), 'r') as tf:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tf, d)
    for file in os.listdir(d):
        gsm = file.split('.')[0].split('_')[0].split('-')[0].upper()
        if gsmList:
            if gsm not in gsmList:
                os.remove('{}/{}'.format(d, file))
                continue
        f = extractGzipFile('{}/{}'.format(d, file), rm=True)
        os.rename(f, '{}/{}.CEL'.format(d, gsm))
        shutil.copy('{}/{}.CEL'.format(d, gsm), '{}/all'.format(celDir))

    exprR += 'setwd("{}/{}")\n' \
             'raw.set = ReadAffy()\n' \
             'rma.data = rma(raw.set)\n' \
             'expr = exprs(rma.data)\n' \
             'write.table(expr, "{}/{}.txt")\n\n'.format(celDir, gse, exprDir, gse)

clinic = pd.concat([pd.read_csv('{}/{}.csv'.format(clinicDir, gse), index_col=0)
                    for gse in gseList])
clinic.to_csv('{}/all.csv'.format(clinicDir))

exprR += 'setwd("{}/all")\n' \
         'raw.set = ReadAffy()\n' \
         'rma.data = rma(raw.set)\n' \
         'expr = exprs(rma.data)\n' \
         'write.table(expr, "{}/all.txt")\n\n'.format(celDir, exprDir)
with open('{}/expr.R'.format(exprDir), 'w') as wf:
    print(exprR, file=wf)


# generate R script for meta-analysis with MAMA package
mkdir(mamaDir)
mamaR = 'rm(list = ls())\nlibrary("MAMA")\nlibrary("metaMA")\n' \
        'library("affyPLM")\nlibrary("affy")\n\n' \
        'memory.limit(16000)\n\n'
for gse in gseListMama:
    mamaR += '{} = read.table("{}/{}.txt")\n'.format(gse, exprDir, gse)
    mamaR += 'cn = sort(colnames({}))\n'.format(gse)
    mamaR += '{} = as.matrix({}[, cn])\n'.format(gse, gse)
mamaR += 'GEDM=list({})\n'.format(', '.join(gseListMama))
mamaR += 'rm({})\ngc()\n\n'.format(', '.join(gseListMama))
mamaR += 'setwd("{}")\n'.format(mamaDir)
for gse in gseListMama:
    mamaR += '{} = read.csv("{}/{}.csv")\n'.format(gse, clinicDir, gse)
    mamaR += '{} = {}[order({}$GSM), ]\n'.format(gse, gse, gse)
    mamaR += 'row.names({}) = {}$GSM\n'.format(gse, gse)
    mamaR += '{} = {}["{}"]\n'.format(gse, gse, tag)
    mamaR += '{}["{}"] = as.factor({}${})\n'.format(gse, tag, gse, tag)
mamaR += 'Clinical = list({})\n'.format(', '.join(gseListMama))
mamaR += 'rm({})\n' \
         'gc()\n\n'.format(', '.join(gseListMama))
mamaR += 'datanames = c({})\n\n'.format(', '.join(['"{}"'.format(gse) for gse in gseListMama]))
mamaR += 'setClass("esets",slots=list(GEDM="list",clinical="list",datanames="character"),package = "MAMA")\n' \
         'meta.data = new("esets", GEDM=GEDM, clinical=Clinical, datanames=datanames)\n' \
         'rm(Clinical, datanames)\ngc()\n\n'
mamaR += 'pval = metaMA(meta.data,"{}",which="pval")\n'.format(tag)
mamaR += 'gc()\n'
mamaR += 'es2 = ES.GeneMeta(meta.data,"{}",nperm=1000)\n\n'.format(tag)
mamaR += 'results1 = join.results(pval, type=1, genenames=rownames(GEDM(meta.data)[[1]]))\n' \
         'p_value = as.data.frame(results1)\n' \
         'rawpval = 2 * (1 - pnorm(abs(pval$TestStatistic)))\n' \
         'FDR_pval = p.adjust(rawpval, method="BY", n=length(rawpval))\n' \
         'p_value$c_pval = rawpval\n' \
         'p_value$FDR = FDR_pval\n' \
         'rm(rawpval, FDR_pval)\n' \
         'gc()\n' \
         'es2_theScores = es2$theScores\n' \
         'es2_ScoresFDR = es2$ScoresFDR\n' \
         'es2_ScoresFDR = es2_ScoresFDR$two.sided\n'
mamaR += 'write.table(p_value, file="{}/p_value", sep="\\t", col.names = T)\n'.format(mamaDir)
mamaR += 'write.table(es2_ScoresFDR, file="{}/es", sep="\\t", col.names = T)\n'.format(mamaDir)

with open('{}/MAMA.R'.format(mamaDir), 'w') as wf:
    print(mamaR, file=wf)


# generate R script for meta-analysis with DEDS package
mkdir(dedsDir)
dedsR = 'rm(list = ls())\n' \
        'library("DEDS")\n\n'
dedsR += 'memory.limit(16000)\n\n'
dedsR += 'expr = read.table("{}/all.txt")\n'.format(exprDir)
dedsR += 'cn = sort(colnames(expr))\n' \
         'expr = as.matrix(expr[, cn])\n' \
         'genes = row.names(expr)\n\n'
dedsR += 'clinic = read.csv("{}/all.csv")\n'.format(clincDir)
dedsR += 'clinic = clinic[order(clinic$GSM), ]\n' \
         'clinic = clinic${}\n\n'.format(tag)
dedsR += 'deds = deds.stat.linkC(expr, clinic, B=1300)\n' \
         'r = topgenes(deds, number=10000, sort.by="fc", genelist=genes)\n' \
         'write.table(r, "{}/deds.txt")'.format(dedsDir)

with open('{}/DEDS.R'.format(dedsDir), 'w') as wf:
    print(dedsR, file=wf)


# generate R script for meta-analysis with limma package
mkdir(limmaDir)
limmaR = 'library(limma)\n\n'
limmaR += 'expr <- read.table("{}/all.txt")\n'.format(exprDir)
limmaR += 'groupList <- read.csv("{}/all.csv")\n'.format(clincDir)
limmaR += 'groupList <- groupList${}\n\n'.format(tag)
limmaR += 'design = model.matrix(~ 0 + factor(groupList))\n' \
          'colnames(design) = c("Control", "Treatment")\n' \
          'fit.lm = lmFit(expr, design)\n' \
          'contrast = makeContrasts("Treatment-Control", levels=design)\n' \
          'fit.cts = contrasts.fit(fit.lm, contrast)\n' \
          'fit.ebs = eBayes(fit.cts)\n' \
          'tempOutput = topTable(fit.ebs, coef=1, n=Inf)\n'
limmaR += 'write.csv(tempOutput, "{}/limma.csv")\n'.format(limmaDir)

with open('{}/limma.R'.format(limmaDir), 'w') as wf:
    print(limmaR, file=wf)


