#!/usr/bin/python

import os

#das_client.py --query "run dataset=/DoubleMuon/Run2016C-PromptReco-v2/AOD | grep run.run_number"

def getRunsInDataset(dataset):
    dbs = '''das_client.py --query "run dataset='''+dataset+''' | grep run.run_number" --limit 999'''
    dbsOut = os.popen(dbs)
    files = []
    for line in dbsOut:
        line = line.rstrip()
        try:
            files.append(int(line))
        except ValueError:
            continue
    return files

def getFillFromRun(run):
    dbs = '''das_client.py --query "run run=%i | grep run.lhcFill"''' % (run)
    dbsOut = os.popen(dbs)
    fill = 0
    for line in dbsOut:
        try:
            fill = int(line)
        except ValueError:
            continue
    return fill
    

if __name__=='__main__':
    runs = getRunsInDataset('/DoubleMuon/Run2016C-PromptReco-v2/AOD')
    for run in runs:
        print '%i:%i' % (run, getFillFromRun(run))
        
