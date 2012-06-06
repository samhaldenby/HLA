import logging
import re

'''
Created on 18 May 2012

@author: sh695
'''

def pick_Nx(hitCounts, proportion=50):
    '''
    Takes raw results of Aligner and returns top hits based on a percentage of the TOTAL score.
    e.g. if proportion = 50, takes top N50 of hits based on total score of each hit
    
    args:
        hitCounts:    results generated from Aligner (dict{geneName:[score,score,score,score])
        proportion:   percent [50]. Function will grab top hits containing [proportion] percent of score. Similar to N50 in assemblies
    
    return:
        top results in dict{geneName:[score,score,score,score]}
    '''
    
    logging.info('pick_Nx: Selecting top results (Results containing top %d%% of total score)', proportion)
    #get total score for entire results set
    totalScore = 0.0
    totalScoreDict = {}
    for entry in hitCounts.items():
        #print entry[0],entry[1]
        scoreForGene = sum(entry[1])  #sum of scores across all reads for this gene
        totalScoreDict[entry[0]] = scoreForGene
        totalScore += scoreForGene
    
    #calculate cut off point
    cutOffScore = totalScore * (proportion/100.0)
    
    #go through results again from top -> bottom scorer and grab results until [proportion] of the total score has been reached
    runningTotal=0.0
    topResults = {}
    drb1RegEx=re.compile("DRB1")
    compensationForOtherDrbs=0
    for entry in sorted(totalScoreDict, key=totalScoreDict.get, reverse=True):
        #print entry,hitCounts[entry],totalScoreDict[entry]
        if drb1RegEx.match(entry) == None:
            compensationForOtherDrbs+=1
        else:
            topResults[entry] = hitCounts[entry]
        runningTotal+=totalScoreDict[entry]
        if runningTotal > cutOffScore:
            if compensationForOtherDrbs == 0:
                break
            elif compensationForOtherDrbs > 0:
                compensationForOtherDrbs-=1

    logging.info('Trimmed results from %d to %d', len(hitCounts), len(topResults))
    return topResults
