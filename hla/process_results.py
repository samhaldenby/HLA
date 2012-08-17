'''
Created on 16 Aug 2012

@author: sh695
'''
import math
import hla.final_results as final_results
from hla.allele_status import Status
import hla.all_ranked_results as AllRankedResults


def process_all_results2(rawResults, results, openList): #results is [{name:score}], openList is [name]
    rejectionList = set() #{geneName}
    for gene in openList:
        resultsForThisGene = results[gene]
        #find genes in openList
        for target in openList:
            if target!=gene: #don't add to rejection list if querying the same gene that was removed from the set for reanalyse - would be daft!
                if target not in resultsForThisGene:
                    print target,"not in",gene 
                    rejectionList.add(target)
                elif resultsForThisGene[target] ==0.0:
                    print target,"has 0 score with",gene
                    rejectionList.add(target)
                else:
                    print target,"has score",resultsForThisGene[target],"with",gene
                    
    for entry in rejectionList:
        print "HERE",entry
        
    #now create accepted list
    newOpenList = []
    for gene in openList:
        if gene not in rejectionList:
            newOpenList.append(gene)
            
    openList = newOpenList
    
    #print "OPENLIST: ", openList
    
    #now remove entries from results so that only those on the open list remain
    newResults = {}
    for entry in results.items():
        geneName = entry[0]
        if geneName in openList:
            newResults[geneName] = entry[1]
    
    #add initial results
    newResults['First']=results['First']
    results = newResults
        
    #create map
    scores = {}
    originalScores = {}
    for gene in openList:
        scores[gene]=[]
        for entry in results.items():
            omittedGene = entry[0]
            result = entry[1]
            #print "Result for %s = %s"%(entry, result)
            if gene in result:  #i.e. if it's not the one being omited, OR!!! if it got no hits.... deal with this later TODO:
                scores[gene].append(result[gene])
            if omittedGene =="First":
                originalScores[gene]= result[gene]

   
    
    #now, make sure all modifier lists are same length (reason they might not be is due to getting 0 hits)
    longestList = 0
    for s in scores.items():
        if len(s[1]) > longestList:
            longestList = len(s[1])
        
    for s in scores.items():
        while len(s[1]) < longestList:
            s[1].append(0.0)

    modifiers = {}
    for entry in scores.items():

        gene = entry[0]
        modifiers[gene]=[]

        for result in entry[1]:
            modifiers[gene].append(result/originalScores[gene])

    #multiply all modifiers by all scores
    modifiedScores = {}
    unmodifiedScores = {}
    #for each gene
    for entry in scores.items():
        thisScore = 1.0
        unmodifiedThisScore = 1.0
        geneName = entry[0]
        geneScores = entry[1]
        #multiply all scores
        for s in geneScores:
            thisScore *= s
            if s>1:
                unmodifiedThisScore *=s
        #mow multiply by modifiers
        

        for m in modifiers[geneName]:
            thisScore *= m
        
        #log it
        if unmodifiedThisScore > 0:
            unmodifiedThisScore = math.log10(unmodifiedThisScore)
        if thisScore >0:
            thisScore = math.log10(thisScore)
            
        #divide by number of top hits
        thisScore /= len(scores)
        unmodifiedThisScore /=len(scores)
        modifiedScores[geneName] = thisScore
        unmodifiedScores[geneName] = unmodifiedThisScore
    
    print "***************"
    for gene in scores.items():
        modifiersForThisOne = modifiers[gene[0]]
        print "%s\t%s\t%s"%(gene[0],'\t'.join(map(str,gene[1])), '\t'.join(map(str,modifiersForThisOne)))
    print "***************"   
    
    prevScore = 0
    scoreNum=0
    finalResults = {} # list of FinalResult class instances
    for entry in sorted(modifiedScores, key=modifiedScores.get, reverse=True):
        scoreNum +=1
        thisScore = modifiedScores[entry]
        status =""
        if thisScore >3:
            #is there too big a difference between this score and the next highest score?
            if prevScore - thisScore > 1:
                status = Status.Fail_score_diff() #"FAIL - Score differential failure"
            #Is this one of the top 2 scorers?
            elif scoreNum <= 2:
                status = Status.Pass() #"PASS"
            else:
                status = Status.Warn_not_top_two() #"WARNING - Pass threshold but not in top 2 alleles"

        else:
            status = Status.Fail_low_score() #"FAIL - Score too low"
        
        
        #print "RAWRESULTS:",rawResults
        #raw_input()
        #print "RESULTS:",results
        #print ">>%s\t%s\t%s"%(entry, modifiedScores[entry], status) #unmodifiedScores[entry],status
        f = final_results.FinalResult()
        f.score = modifiedScores[entry]
        f.targetName = entry
        f.status = status
        f.rawCounts = rawResults[entry]
        finalResults[f.score]=(f)
        
        print ">>%s\t%s\t%s\t%s"%(entry,modifiedScores[entry],status,rawResults[entry])
        prevScore = thisScore
    resultsBundle = final_results.FinalResultBundle()
    resultsBundle.results = finalResults
            
    return resultsBundle        





def get_overlap_score(openList): #results is [{name:score}], openList is [name]
    results = AllRankedResults.results
    rejectionList = set() #{geneName}
    for gene in openList:
        resultsForThisGene = results[gene]
        #find genes in openList
        for target in openList:
            if target!=gene: #don't add to rejection list if querying the same gene that was removed from the set for reanalyse - would be daft!
                if target not in resultsForThisGene:
                    print target,"not in",gene 
                    rejectionList.add(target)
                elif resultsForThisGene[target] ==0.0:
                    print target,"has 0 score with",gene
                    rejectionList.add(target)
                else:
                    print target,"has score",resultsForThisGene[target],"with",gene
                    
    for entry in rejectionList:
        print "HERE",entry
        
    #now create accepted list
    newOpenList = []
    for gene in openList:
        if gene not in rejectionList:
            newOpenList.append(gene)
            
    openList = newOpenList
    
    #print "OPENLIST: ", openList
    
    #now remove entries from results so that only those on the open list remain
    newResults = {}
    for entry in results.items():
        geneName = entry[0]
        if geneName in openList:
            newResults[geneName] = entry[1]
    
    #add initial results
    newResults['First']=results['First']
    results = newResults
        
    #create map
    scores = {}
    originalScores = {}
    for gene in openList:
        scores[gene]=[]
        for entry in results.items():
            omittedGene = entry[0]
            result = entry[1]
            #print "Result for %s = %s"%(entry, result)
            if gene in result:  #i.e. if it's not the one being omited, OR!!! if it got no hits.... deal with this later TODO:
                scores[gene].append(result[gene])
            if omittedGene =="First":
                originalScores[gene]= result[gene]

   
    
    #now, make sure all modifier lists are same length (reason they might not be is due to getting 0 hits)
    longestList = 0
    for s in scores.items():
        if len(s[1]) > longestList:
            longestList = len(s[1])
        
    for s in scores.items():
        while len(s[1]) < longestList:
            s[1].append(0.0)

    modifiers = {}
    for entry in scores.items():

        gene = entry[0]
        modifiers[gene]=[]

        for result in entry[1]:
            modifiers[gene].append(result/originalScores[gene])

    #multiply all modifiers by all scores
    modifiedScores = {}
    unmodifiedScores = {}
    #for each gene
    overlapScore = {}
    for entry in scores.items():
        totScoreForGene =0.0
        thisScore = 1.0
        unmodifiedThisScore = 1.0
        geneName = entry[0]
        geneScores = entry[1]
        #multiply all scores
        for s in geneScores:
            thisScore *= s
            if s>1:
                unmodifiedThisScore *=s
        #mow multiply by modifiers
        
        print "--- Modifier for %s"%(geneName)
        for m in modifiers[geneName]:
            print "\t--- %.2f"%(m) 
            
            thisScore *= m
            totScoreForGene +=m
        
        #deduct one from score (the 1 comes from self crossing being a 1 - don't care about this so get rid
        totScoreForGene -=1
        print "--- Tot Score for %s = %.2f"
        
        overlapScore[geneName] = totScoreForGene
        
    #convert overlapScores to list
   
    overlapScoreLists = []
    for o in overlapScore.items():
        overlapScoreLists.append({"name":o[0],"oScore":o[1]})
    return overlapScoreLists

#        #log it
#        if unmodifiedThisScore > 0:
#            unmodifiedThisScore = math.log10(unmodifiedThisScore)
#        if thisScore >0:
#            thisScore = math.log10(thisScore)
#            
#        #divide by number of top hits
#        thisScore /= len(scores)
#        unmodifiedThisScore /=len(scores)
#        modifiedScores[geneName] = thisScore
#        unmodifiedScores[geneName] = unmodifiedThisScore
#    
#        
#    print "***************"
#    for gene in scores.items():
#        modifiersForThisOne = modifiers[gene[0]]
#        print "%s\t%s\t%s"%(gene[0],'\t'.join(map(str,gene[1])), '\t'.join(map(str,modifiersForThisOne)))
#    print "***************"   
#    
#    
#    returnBundle = list()
#    prevScore = 0
#    scoreNum=0
#    finalResults = {} # list of FinalResult class instances
#    for entry in sorted(modifiedScores, key=modifiedScores.get, reverse=True):
#        scoreNum +=1
#        thisScore = modifiedScores[entry]
#        status =""
#        if thisScore >3:
#            #is there too big a difference between this score and the next highest score?
#            if prevScore - thisScore > 1:
#                status = Status.Fail_score_diff() #"FAIL - Score differential failure"
#            #Is this one of the top 2 scorers?
#            elif scoreNum <= 2:
#                status = Status.Pass() #"PASS"
#            else:
#                status = Status.Warn_not_top_two() #"WARNING - Pass threshold but not in top 2 alleles"
#
#        else:
#            status = Status.Fail_low_score() #"FAIL - Score too low"
#        
#        
#      
#        #raw_input()
#        #print "RESULTS:",results
#        #print ">>%s\t%s\t%s"%(entry, modifiedScores[entry], status) #unmodifiedScores[entry],status
#        f = final_results.FinalResult()
#        f.score = modifiedScores[entry]
#        f.targetName = entry
#        f.status = status
#      
#        finalResults[f.score]=(f)
#        
#        print ">>&& %s\t%s\t%s"%(entry,modifiedScores[entry],status)
#        returnBundle.append({"name":entry, "newScore":modifiedScores[entry]})
#        prevScore = thisScore
#    resultsBundle = final_results.FinalResultBundle()
#    resultsBundle.results = finalResults
#            
#    return returnBundle       
