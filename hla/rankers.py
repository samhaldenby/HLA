import logging
import math

'''
Created on 18 May 2012

@author: sh695
'''

def rank_results(results):
    """Ranks genes based on number of each of the four reads/amplicons that hit the gene
    
               
        Step 1: The results are summed
        Gene        Counts             Sum
        DRB1*01:01  [120,10 ,500,0 ]   195
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5
        DRB1*03:04  [10 ,100,50 ,35]   95
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25
    
    
        Step 2: The scores are weighted if any reads reported no hits, i.e. weighting of 1 for no misses, 0.75 for 1, 0.5 for 2, 0.25 for 3 misses
        Gene        Counts             Sum    NumMisses    MissWeighting
        DRB1*01:01  [120,10 ,500,0 ]   195    1              0.75
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5   2              0.5
        DRB1*03:04  [10 ,100,50 ,35]   95     0              1
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25   2              0.5
        
        
        Step 3: The scores are then weighted by multiplying all the (read totals+1) together and log10-ing it
        Gene        Counts             Sum    MissWeighting    ReadWeighting
        DRB1*01:01  [120,10 ,500,0 ]   195    0.75             5.82
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5   0.5              3.05
        DRB1*03:04  [10 ,100,50 ,35]   95     1                6.31
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25   0.5              1.48
        
        
        Step 4: Multiply Sum x MissWeighting x ReadWeighting
        Gene        Counts             Sum    MissWeighting    ReadWeighting    Final Score
        DRB1*01:01  [120,10 ,500,0 ]   195    0.75             5.82             851.18
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5   0.5              3.05             41.94
        DRB1*03:04  [10 ,100,50 ,35]   95     1                6.31             599.45
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25   0.5              1.48             2.41
        
        So in this example, your best scorers are 01:01 and 03:04
        

        Args:
            results: A dict mapping keys as gene name and values as lists (size 4) containing the hit count for each read.
    
        Returns:
            dict of {geneTargets:scores}
    
        Raises:
            Nothing. It probably should though!
    """    
    
    logging.info('rank_results: Ranking %d results',len(results))
    
    

    newScoreMap = {}
    
    #for each result (e.g. {DRB1*01:01, [100,200,0,300]} calculate an overall score and rank it
    for entry in results.items():
        name = entry[0]
        originalCounts = entry[1]

        #Step 1: Calculate total
        summedScore = sum(originalCounts)
        
        #Step 2: Calculate miss modifier
        misses =0
        missModifier = 1.0        
        penaltyPerMiss = 1.0 / (len(originalCounts)*1.0)
        for i in range(0,len(originalCounts)):
            if originalCounts[i]==0:
                missModifier-= penaltyPerMiss
                misses+=1
                
        #Step 3: Calculate log-multiple modifier
        mod = 1.0
        for i in range(0,len(originalCounts)):
            mod*=(originalCounts[i]+1)  #add 1 so that you're not multiplying by 0 TODO: Does this need tweaking?
        modlog = math.log10(mod)
        
        #Step 4: Multiply together
        newScoreMap[name] = summedScore * missModifier * modlog
    
    
    #Sort results in order and return
    returnResults = {}
    for entry in sorted(newScoreMap, key=newScoreMap.get, reverse=True):
        name = entry
        score = newScoreMap[entry]
        print "%s\t%s\t%.2f"%(name,'\t'.join(map(str,results[name])),score)
        returnResults[name]= score

        
    return returnResults
            

