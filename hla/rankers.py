import logging
import math

'''
Created on 18 May 2012

@author: sh695
'''

def rank_by_cross_division(results):
    #TODO:make this accept more than 4 reads!
    """Ranks genes based on number of each of the four reads that hit the gene
    
        Step1: For each read (2_2_F, 2_4_F, 2_5_R, 2_4_R), each hit gene will be ranked according to how many hits it received. For example, if 2_2_F hit the following genes in this pattern, they would be ranked as shown in column 3
        Gene        Count   Rank
        DRB1*01:01  120     1
        DRB1*03:02  10      3
        DRB1*03:04  10      3
        DRB1*04:01  5       4
        
        This is repeated for all reads and genes.
        
        Gene        Counts             Ranks
        DRB1*01:01  [120,10 ,500,0 ]   [1,2,1,4]
        DRB1*03:02  [10 ,0  ,100,0 ]   [3,4,2,4]
        DRB1*03:04  [10 ,100,50 ,35]   [3,1,3,1]
        DRB1*04:01  [5  ,0  ,0   4 ]   [4,4,4,2]
        
        
        Step 2: For each gene, the best count is divided by the worst rank, the second best count is divided by the second worst rank, etc 
        Gene        Counts             Ranks      Scores
        DRB1*01:01  [120,10 ,500,0 ]   [1,2,1,4]  [60  ,10  ,125 ,0 ] (i.e. 500/4, 120/2, 10/1, 0/1)
        DRB1*03:02  [10 ,0  ,100,0 ]   [3,4,2,4]  [2.5 ,0   ,25  ,0 ]
        DRB1*03:04  [10 ,100,50 ,35]   [3,1,3,1]  [10  ,33.3,16.7,35]
        DRB1*04:01  [5  ,0  ,0   4 ]   [4,4,4,2]  [1.25,0   ,0   ,2 ]
        
        
        Step 3: The results are then summed
        Gene        Counts             Ranks      Scores               Sum
        DRB1*01:01  [120,10 ,500,0 ]   [1,2,1,4]  [60  ,10  ,125 ,0 ]  195
        DRB1*03:02  [10 ,0  ,100,0 ]   [3,4,2,4]  [2.5 ,0   ,25  ,0 ]  27.5
        DRB1*03:04  [10 ,100,50 ,35]   [3,1,3,1]  [10  ,33.3,16.7,35]  95
        DRB1*04:01  [5  ,0  ,0   4 ]   [4,4,4,2]  [1.25,0   ,0   ,2 ]  3.25
    
    
        Step 4: The scores are weighted if any reads reported no hits, i.e. weighting of 1 for no misses, 0.75 for 1, 0.5 for 2, 0.25 for 3 misses
        Gene        Counts             Sum    NumNoMisses    MissWeighting
        DRB1*01:01  [120,10 ,500,0 ]   195    1              0.75
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5   2              0.5
        DRB1*03:04  [10 ,100,50 ,35]   95     0              1
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25   2              0.5
        
        
        Step 5: The scores are then weighted by multiplying all the (read totals+1) together and log10-ing it
        Gene        Counts             Sum    MissWeighting    ReadWeighting
        DRB1*01:01  [120,10 ,500,0 ]   195    0.75             5.82
        DRB1*03:02  [10 ,0  ,100,0 ]   27.5   0.5              3.05
        DRB1*03:04  [10 ,100,50 ,35]   95     1                6.31
        DRB1*04:01  [5  ,0  ,0   4 ]   3.25   0.5              1.48
        
        
        Step 6: Multiply Sum x MissWeighting x ReadWeighting
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
    
    logging.info('logging_by_cross_division: Ranking %d results',len(results))
    
    
    #Manipulate data into 4 lists of dicts{geneName: countForThatRead}
    countsPerRead = [{}] *4
    for readNum in range(0,4):
        #print readNum
        countsPerRead[readNum]={}
        for entry in results.items():
            name = entry[0]
            counts = entry[1]
            countsPerRead[readNum][name] = counts[readNum]
            #print readNum,countsPerRead[readNum][name]
        
        
    #now rank each gene based on number of reads mapping
    finalRanks = {}
    for readNum in range(0,4):

        #read through all readCounts in sorted order and assign a rank
        currRank =1
        prevScore = -1
        allSameScoreList = []
        rankDifferential = 0 
        for entry in sorted(countsPerRead[readNum], key=countsPerRead[readNum].get, reverse=True): 
            score = countsPerRead[readNum][entry]
            name = entry 
            if score != prevScore:
                currRank = currRank+rankDifferential     
                rankDifferential=1
                allSameScoreList = []
                allSameScoreList.append(name)
                #add rank
                if name not in finalRanks:
                    finalRanks[name] = [0]*4
                finalRanks[name][readNum]= currRank                
            else:
 
                if name not in finalRanks:
                    finalRanks[name] = [0]*4
                finalRanks[name][readNum]= rankDifferential + currRank  
                allSameScoreList.append(name)
                for name in allSameScoreList:
                    finalRanks[name][readNum]=rankDifferential + currRank
                #add rank
                rankDifferential+=1
            prevScore = score

            
    #calculate raw score
    overallRanks = {}
    for entry in finalRanks.items():
        #print entry
        name = entry[0]
        rankings = entry[1]
        overallRanks[name] = sum(rankings)
    
    newScoreMap = {}
    for entry in sorted(overallRanks, key=overallRanks.get, reverse=True):
        name = entry
        rank = overallRanks[entry]
        individualRanks = finalRanks[name]
        originalCounts = results[name]

        
        #new ranking idea multiply highest count by worst rank.....lowest count by best rank
        orderedRanks=[]
        orderedCounts=[]
        for i in sorted(individualRanks):
            orderedRanks.append(i)
        for i in sorted(originalCounts):
            orderedCounts.append(i)
        misses =0
        newIdeaScore = 0.0
        for i in range(0,4):
            if orderedCounts[i]>0:
                newIdeaScore+=orderedCounts[i]/orderedRanks[i]
            else:
                misses+=1
            
        newIdeaScore = newIdeaScore/(1+misses)
        
       
        calibratedScore = 0
        missModifier = 1.0        
        for i in range(0,4):
            calibratedScore+=originalCounts[i]/individualRanks[i]
            if originalCounts[i]==0:
                missModifier-=0.25
                misses+=1
                
        newIdeaScore = newIdeaScore/(1+misses)
        calibratedScore = calibratedScore / missModifier
        #TODO this was original = 141
#        newScoreMap[name]=newIdeaScore
        #END 
        
        
        #TODO (multiplied by sum(original counts) for results 2
#        newScoreMap[name]=newIdeaScore* sum(originalCounts) = 146
        #END
        
        
        #TODO this is what was done for results 3 = 143
#        mod = 1.0
#        for i in range(0,4):
#            mod*=(originalCounts[i]+1)
#        modlog = math.log10(mod)
#        newScoreMap[name]=newIdeaScore*(modlog)
        #END
        
        
        #TODO this is for results 4 = 147
#        mod = 1.0
#        for i in range(0,4):
#            mod*=(originalCounts[i]+1)
#        modlog = math.log10(mod)
#        newScoreMap[name] = modlog * sum(originalCounts)       
        #END
        
        
        #TODO this is for results 5 = 149
        mod = 1.0
        for i in range(0,4):
            mod*=(originalCounts[i]+1)
        modlog = math.log10(mod)
        newScoreMap[name] = modlog * sum(originalCounts) * missModifier  
        #END
        
        
        #    print "size",len(newScoreMap), name, newIdeaScore, newScoreMap[name]
        #add to map based on newIdeaScore
        
        #print name,"\t",individualRanks,"\t",originalCounts,"\t", rank,"\t", sum(originalCounts)/rank,"\t", calibratedScore,"\t", newIdeaScore
    
    top10 =0
   # for entry in newScoreMap.items():
      #  print "@",entry
    #for entry in sorted(finalScores, key=finalScores.get, reverse=True):
    returnResults = {}
    for entry in sorted(newScoreMap, key=newScoreMap.get, reverse=True):
      #  print "ENTRY:",entry, len(newScoreMap)
        name = entry
        score = newScoreMap[entry]
        
        if top10 < 1000:
            print name,finalRanks[name],results[name],score
            returnResults[name]= score
        top10+=1
        
    return returnResults
            
    
    #translate back to score map i.e. (name,[1][10][13][1])
#    for entry in finalRanks.items():
#        print entry[0],entry[1]
