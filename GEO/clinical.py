##############################################################################
#
# Copyright (C) 2018 Minzhang Cheng
# Contact: minzhangcheng@gmail.com
#
# GNU Lesser General Public License Usage
# This file may be used under the terms of the GNU Lesser General Public
# License version 3 as published by the Free Software Foundation and
# appearing in the file LICENSE included in the packaging of this file.
# Please review the following information to ensure the GNU Lesser
# General Public License version 3 requirements will be met:
# http://www.gnu.org/licenses/gpl-3.0.html
#
##############################################################################

##############################################################################
# Fetch and summary clinical data of each datasets.
#
# Only works on Python 3.3+
# Only works on linux (Debian 10 was used) as some characters not supported
#   in windows console.
##############################################################################


from Bio import Entrez as e
import multiprocessing
import multiprocessing.dummy
import json
import requests

from GEO.utility import mkdir

def search(db, term, email, save='', max=0, thread=8, once=50, start=0):
    e.email = email
    if max != 0 and max < once:
        once = max
    try:
        handle = e.esearch(db=db, term=term, retmax=once, retstart=start)
        record = e.read(handle)
    except:
        return search(db, term, email, save, max, thread, once, start)
    if max ==0 or max > int(record['Count']):
        max = int(record['Count'])

    if max != 0 and max <= once:
        if save:
            with open(save, 'w') as wf:
                json.dump(record, wf, indent=2)
        return record['IdList']

    r = record['IdList']

    if thread > 1:
        p = [(db, term, email, '', once, 1, once, i) for i in range(once, max, once)]
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(search, p)
        for i in results:
            r.extend(i)
    else:
        for i in range(once, max, once):
            r.extend(search(db, term, email, '', once, thread, once, i))
    r = r[0:max]
    if save:
        with open(save, 'w') as wf:
            json.dump(r, wf, indent=2)
    return r[0:max]


def summary(db, ids, email, save='', thread=8, once=20):
    e.email = email
    if len(ids) <= once:
        try:
            handle = e.esummary(db=db, id=','.join(ids), retmax=once)
            record = e.read(handle)
        except:
            return summary(db, ids, email, save, thread, once)
        if save:
            with open(save, 'w') as wf:
                json.dump(record, wf, indent=2)
        return record

    groupedIds = [ids[i:i+once] for i in range(0, len(ids), once)]
    r = list()
    if thread > 1:
        p = [(db, i, email, '', 1, once) for i in groupedIds]
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(summary, p)
        for i in results:
            r.extend(i)
    else:
        for i in groupedIds:
            r.extend(summary(db, i, email, '', thread, once))
    if save:
        with open(save, 'w') as wf:
            json.dump(r, wf, indent=2)
    return r


def _gseClinical(gdsInfo, dataTableDir='', saveJsonDir='',
                 valueKeep=True, download=True):
    gdsId = gdsInfo['Accession']
    id = gdsInfo['GSE']
    id = 'GSE%s' % id
    gpl = gdsInfo['GPL']
    if ';' in gpl:
        gplList = ['GPL%s' % i for i in gpl.split(';')]
        gpl = ';'.join(gplList)
    else:
        gpl = 'GPL%s' % gpl
        gplList = [gpl]
    returnValue = list()
    for gpl in gplList:
        if download:
            url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/'
            r = requests.get(url, {'acc': id})
            gse = json.loads(r.text)
            gpl2 = gse['GeoMetaData'][0]['entity']['series']['platforms']
            if gpl in gpl2:
                r = requests.get(url, {'type': 'samples', 'series': id, 'platform': gpl})
                samples = json.loads(r.text)['GeoMetaData']
                if saveJsonDir:
                    with open('%s/%s_%s_clinical.json' % (saveJsonDir, id, gpl), 'w') as wf:
                        json.dump(samples, wf, indent=2)
        else:
            try:
                with open('%s/%s_%s_clinical.json' % (saveJsonDir, id, gpl), 'r') as rf:
                    samples = json.load(rf)
            except:
                continue
        sample_count = len(samples)
        tags = dict()
        data = dict()
        for sample in samples:
            gsm = sample['acc']
            data[gsm] = dict()
            valid_sample = False
            if 'entity' in sample:
                if 'sample' in sample['entity']:
                    if 'channels' in sample['entity']['sample']:
                        if sample['entity']['sample']['channels']:
                            if 'characteristics'in sample['entity']['sample']['channels'][0]:
                                valid_sample = True
            if valid_sample:
                c = sample['entity']['sample']['channels'][0]['characteristics']
                for i in c:
                    if 'tag' in i:
                        tag = i['tag'].strip()
                        value = i['value'].strip()
                        data[gsm][tag] = value
                        if tag not in tags:
                            tags[tag] = dict()
                        if value not in tags[tag]:
                            tags[tag][value] = 0
                        tags[tag][value] += 1
                    else:
                        if 'value' in i:
                            if ':' in i['value']:
                                tag, value = i['value'].split(':', 1)
                                tag = tag.strip()
                                value = value.strip()
                                data[gsm][tag] = value
                                if tag not in tags:
                                    tags[tag] = dict()
                                if value not in tags[tag]:
                                    tags[tag][value] = 0
                                tags[tag][value] += 1
        tagList = tags.keys()
        if dataTableDir:
            s = 'GSM\t'
            s += '\t'.join(tagList)
            for sample in data:
                s += '\n%s' % sample
                for tag in tagList:
                    if tag in data[sample]:
                        s += '\t%s' % data[sample][tag]
                    else:
                        s += '\t'
            with open('%s/%s_%s_clinical.tsv' % (dataTableDir, id, gpl), 'w') as wf:
                print(s, file=wf, end='')
        if valueKeep:
            returnValue.append((gdsId, id, gpl, sample_count, ['%s:%s' % (tag, tags[tag]) for tag in tagList]))
        else:
            returnValue.append((gdsId, id, gpl, sample_count, tagList))
    return returnValue


def getClinical(gdsInfos, dataTableDir='', saveClinicalTag='',
                saveJsonDir='', valueKeep=True, thread=8, download=True):
    p = [(gdsInfo, dataTableDir, saveJsonDir, valueKeep, download)
         for gdsInfo in gdsInfos]
    with multiprocessing.dummy.Pool(processes=thread) as pool:
        r = pool.starmap(_gseClinical, p)
    r = [j for i in r for j in i]
    s = 'GDS id\tGSE id\tGPL id\tSample Count\tCharacteristics Tags\n'
    for tag in r:
        s += '%s\t%s\t%s\t%d\t' % (tag[0], tag[1], tag[2], tag[3])
        s += '\t'.join(tag[4])
        s += '\n'
    if saveClinicalTag:
        with open(saveClinicalTag, 'w') as wf:
            print(s, file=wf, end='')
    return r


def fetchClinicalData(searchTerm, resultDir, email, min_n_samples=20,
                      gdsType='Expression profiling by array',
                      thread=8, download=True, saveJson=True,
                      saveTable=True, valueKeep=True):

    """
    Search and fetch GEO datasets according to the key words, and then obtain
      all the clinical data contained in these datasets.
    :param searchTerm: the key words, "breast+cancer" as an example
    :param resultDir: the directory to store the results
    :param email: your e-mail
    :param min_n_samples: the datasets containing less samples would be ignored
    :param gdsType: the type of dataset, 'Expression profiling by array' as default
    :param thread: thread to be used for downloading data
    :param download: keep it as true
    :param saveJson: keep it as true
    :param saveTable: keep it as true
    :param valueKeep: keep it as true
    :return: None
    """

    mkdir(resultDir)
    mkdir('%s/gse_clinical_tsv' % resultDir)
    mkdir('%s/gse_clinical_json' % resultDir)

    # search and summary GDS
    gdsIdList = search('gds', searchTerm, email,
                       save='%s/gds_id.json' % resultDir)
    gdsInfos = summary('gds', gdsIdList, email,
                       save='%s/gds_summary.json' % resultDir, once=20)

    # filter
    filteredGds = list()
    for gds in gdsInfos:
        if gds['gdsType'] == gdsType:
            if gds['n_samples'] > min_n_samples:
                if gds['taxon'] == 'Homo sapiens':
                    filteredGds.append(gds)
    with open('%s/filtered_gds_summary.json' % resultDir, 'w') as wf:
        json.dump(filteredGds, wf, indent=4)
    print(len(filteredGds))

    # fetch detailed clinical data tags
    with open('%s/filtered_gds_summary.json' % resultDir, 'r') as rf:
        gdsInfos = json.load(rf)

    getClinical(gdsInfos, dataTableDir='%s/gse_clinical_tsv' % resultDir,
                saveClinicalTag='%s/gds_clinical_tags.txt' % resultDir,
                saveJsonDir='%s/gse_clinical_json' % resultDir,
                valueKeep=True, thread=8)
