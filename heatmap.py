import random
import pandas as pd
import os

rootDir = 'C:/Users/minzhang/Downloads'
exprFile = '{}/expr.txt'.format(rootDir)
clinicFile = '{}/clinic.csv'.format(rootDir)
heatmapDir = '{}/heatmap'.format(rootDir)


def randomSamples(clinicFile, count, key='', outputFile='', seed=0,
                  stratification='equal'):
    results = list()
    clinic = pd.read_csv(clinicFile, index_col=0)
    if not key:
        key = clinic.columns[-1]
    clinic = pd.DataFrame(clinic, columns=[key])
    clinic.dropna(axis=0, how='any')

    if stratification == 'random':
        random.seed(seed)
        sampleSelected = random.sample(list(clinic.index), count)
    else:
        groups = set(clinic[key])
        groupCount = {
            i: len(clinic[clinic[key] == i])
            for i in groups
        }
        allCount = len(clinic)
        countSelected = dict()
        if stratification == 'ratio':
            countSelected = {
                i: round(groupCount[i] / allCount * count)
                for i in groupCount
            }
        elif stratification == 'equal':
            countLeft = count
            groupLeft = len(groups)
            countSorted = sorted(groupCount.items(), key=lambda x: x[1])
            for i in countSorted:
                avg = countLeft / groupLeft
                k, value = i
                if value < avg:
                    countSelected[k] = value
                    countLeft = countLeft - value
                else:
                    countSelected[k] = round(avg)
                    countLeft = countLeft - round(avg)
                groupLeft = groupLeft - 1
        sampleSelected = list()
        for i in groups:
            random.seed(seed)
            population = list(clinic[clinic[key] == i].index)
            selected = random.sample(population, countSelected[i])
            sampleSelected.extend(selected)

    clinicSelected = pd.DataFrame(clinic, index=sampleSelected)
    if outputFile:
        clinicSelected.to_csv(outputFile)

    return clinicSelected


def exprSelected(exprFile, probeSelected, sampleSelected, outputFile=''):
    expr = pd.read_table(exprFile, index_col=0, sep=' ')
    exprSelected = pd.DataFrame(expr,
                                index=probeSelected,
                                columns=sampleSelected)
    if outputFile:
        exprSelected.to_csv(outputFile)
    return exprSelected


def randomHeatmap(clinicFile, exprFile, probeSelected, runningDir, count,
                  symbolSelected=[], key='', seed=0, stratification='equal',
                  runningR=False):
    clinicSelectedFile = '{}/clinic_random_{}.csv'.format(runningDir, seed)
    exprSelectedFile = '{}/expr_random_{}.csv'.format(runningDir, seed)

    clinicSelected = randomSamples(clinicFile, count, key,
                                   clinicSelectedFile, seed, stratification)
    sampleSelected = list(clinicSelected.index)
    expr = exprSelected(exprFile, probeSelected,
                        sampleSelected, exprSelectedFile)

    rScript = 'library(gplots)\nlibrary(pheatmap)\nlibrary(RColorBrewer)\n\n'
    rScript += '\nexpr <- read.csv("{}", row.names=1)\n'.format(exprSelectedFile)
    if symbolSelected:
        rScript += 'symbols = c("{}")\n'.format('", "'.join(symbolSelected))
        rScript += 'row.names(expr) <- symbols\n'
    rScript += 'expr <- log(expr + 1, base=2)\n'
    rScript += 'x <- as.matrix(expr)\n'
    rScript += 'x <- t(scale(t(x)))\n'

    rScript += 'png(file="{}/pheatmap_random_{}.png", width=1024, height=1024, bg="transparent")\n' \
        .format(runningDir, seed)
    rScript += 'pheatmap(x, cutree_rows=2, cutree_cols=2, color=greenred(75), border_color=NA)\n'
    rScript += 'dev.off()\n'
    rScript += 'png(file="{}/heatmap2_random_{}.png", width=1024, height=1024, bg="transparent")\n' \
        .format(runningDir, seed)
    rScript += 'heatmap.2(x, col=greenred, scale="row", trace="none")\n'
    rScript += 'dev.off()\n'

    rScript += '\nclinic <- read.csv("{}", head=T, row.names=1)\n'.format(clinicSelectedFile)
    rScript += 'c <- clinic\n'
    rScript += 'annotation_c <- data.frame(c)\n'
    rScript += 'rownames(annotation_c) <- colnames(x)\n'
    for i in ['correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski']:
        rScript += 'png(file="{}/pheatmap2_random_{}_{}.png", width=1024, height=1024, bg="transparent")\n' \
            .format(runningDir, seed, i)
        rScript += 'pheatmap(as.matrix(x), annotation_col=annotation_c, color=bluered(200), border_color=NA, cutree_rows=2, cutree_cols=2, clustering_distance_cols="{}", scale="column")\n' \
            .format(i)
        rScript += 'dev.off()\n'

    rFile = '{}/heatmap_{}.R'.format(runningDir, seed)
    with open(rFile, 'w') as wf:
        print(rScript, file=wf)

    if runningR:
        if runningR:
            cmd = 'C:\\Open\\R\\bin\\x64\\Rscript.exe {} > {}.log'.format(rFile, rFile)
            r = os.system(cmd)


def randomHeatmapMultiple(clinicFile, exprFile, probeSelected, runningDir,
                          count, symbolSelected=[], repeat=10, key='', seed=[],
                          stratification='equal', runningR=False):
    if not seed:
        seed = list(range(0, repeat))
    for s in seed:
        randomHeatmap(clinicFile, exprFile, probeSelected, runningDir,
                      count, symbolSelected, key, s, stratification, runningR)


probes = ['238578_at', '42361_g_at', '244383_at', '243736_at', '224686_x_at', '227477_at', '37425_g_at', '243209_at', '243149_at', '238706_at', '239597_at', '242407_at', '238576_at', '225657_at', '242770_at', '240052_at', '31807_at', '215942_s_at', '218586_at', '238587_at', '239673_at', '227371_at', '239802_at', '242143_at', '32042_at', '238462_at', '234788_x_at', '243561_at', '221906_at', '239886_at', '242787_at', '225612_s_at', '219490_s_at', '229551_x_at', '218868_at', '223492_s_at', '209825_s_at', '224712_x_at', '220018_at', '227484_at', '238712_at', '232113_at', '239232_at', '243584_at', '241472_at', '238595_at', '242769_at', '221747_at', '230419_at', '238075_at']
randomHeatmap(clinicFile, exprFile, probes, heatmapDir, 50, seed=7, stratification='equal', runningR=True)
