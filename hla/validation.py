import logging
import hla.align as align
from hla.process_results import *

'''
Created on 21 May 2012

@author: sh695
'''


PROP_DIFF_IN_OVERLAP = 4
FOLD_DIFF_IN_COUNTS = 4
def compare_with_cam2(resultsBundle, wellId, cam2FileName):
    logging.info('Comparing to Cam2 results')
    
    print ">R> %s ***"%wellId
    
    try:
        cam2File = open(cam2FileName)
    except IOError:
        logging.error("Unable to open %s",cam2FileName)
    
    
    #grab correct line from dynal results
    dynalResultsLine = ""
    for line in cam2File:
        if wellId in line:
            dynalResultsLine = line.strip()
            
    dynalResultsLineSplit = dynalResultsLine.split("\t")
    firstDynalHit = dynalResultsLineSplit[1]
    secondDynalHit = dynalResultsLineSplit[2]
    
    #compile dynal results
    #dynalResults = set()
    #dynalResults.add(firstDynalHit)
    #dynalResults.add(secondDynalHit)
    dynalResults = list()
    dynalResults.append(firstDynalHit)
    dynalResults.append(secondDynalHit)
    
    print ">R>",dynalResults
    
    
    score = 0
    for allele in resultsBundle.top_alleles:
        found = False
        #trim down to major type
        majorType = "CT"
        print "allele: ",allele
    
        majorType = allele[:2]
        print "Major type: %s -> %s"%(allele,majorType) 
        for i in range(0,len(dynalResults)):
            print "DYNAL_RESULTS: ",dynalResults
            if majorType in dynalResults[i]:
                dynalResults[i] = "DONE"
                
                found = True
                score+=1
                print ">R> %s\tYes"%allele
                break
        if not found:
            print ">R> %s\tNo"%allele
    print ">R>>\t%s\t%s"%(wellId,score)
        
#def compare_with_dynal(resultsBundle, wellId, dynalFileName):
#    
#    logging.info('Comparing to dynal results')
#                 
#    print "*** %s ***"%wellId
#                 
#    #open dynal results file
#    try:
#        dynalFile = open(dynalFileName)
#    except IOError:
#        logging.error("Unable to open %s",dynalFileName)
#    
#    
#    #grab correct line from dynal results
#    dynalResultsLine = ""
#    for line in dynalFile:
#        if wellId in line:
#            dynalResultsLine = line.strip()
#            
#    dynalResultsLineSplit = dynalResultsLine.split("\t")
#    firstDynalHit = dynalResultsLineSplit[1]
#    secondDynalHit = dynalResultsLineSplit[2]
#    
#    #compile dynal results
#    dynalResults = set()
#    dynalResults.add(firstDynalHit)
#    dynalResults.add(secondDynalHit)
#    print ">>",dynalResults
#    
#    
#    score = 0
#    for allele in resultsBundle.accepted_alleles:
#        if allele in dynalResults:
#            score+=1
#            print ">> %s\tYes"%allele
#        else:
#            print ">> %s\tNo"%allele
#            
#    cam2File.close()
    
        
    
#    #compile these results
#    myResults = set()
#    score = 0
#    top2 = 0
#    for entry in sorted(results, key=results.get, reverse=True):
#      #  print "ENTRY:",entry, len(newScoreMap)
#        name = entry
#        #score = results[entry]
#        
#        if top2 < 2:
#            modifiedName = "%s%s"%(name[5:7],name[8:])
#            
#            if modifiedName in dynalResults:
#                score+=1
#                print ">> %s\tYes"%modifiedName
#            else:
#                print ">> %s\tNo"%modifiedName
#        top2+=1
#        
#    print ">> Score:\t\t\t\t\t\t\t\t\t%d"% score
        
def polish_results(finalBundle):
    #first, process bundle into corresponding genes, major types, counts and scores
    scores= []
    names = []
    majorTypes = []
    countsPerAmp = []
    for entry in sorted(finalBundle.results.iteritems(), key=lambda (x,y): float(x), reverse=True):
        score = entry[0]
        result = entry[1]
        scores.append(score)
        names.append(result.targetName)
        majorTypes.append(result.targetName.split(':')[0])
        countsPerAmp.append(result.rawCounts)
        print score,result.targetName,result.targetName.split(':')[0]
        
        
    
    bundle = []
    for na,sc,mt,cpa in zip(names,scores,majorTypes,countsPerAmp):
        bundle.append({"name" : na, "score" : sc, "majorType" : mt, "countsPerAmp" : cpa})
        
    #Now, go through all results in bundle and determine whether any merging (due to ambiguity) should occur
    print "-> Will now compare [0] bundles with [1] bundles iteratively"
    while len(bundle)>1:
        action = determine_action(bundle[0],bundle[1])
        print "Returned Action: %s"%action
        if action == "sep":
            break
        elif action == "merge":
            bundle[0] = merge_bundles(bundle[0],bundle[1])
            del bundle[1]
        elif "del" in action:
            nameOfBundleToDelete = action.split('\t')[1]
            bundleNum = 0
            print "Going to delete %s"%nameOfBundleToDelete
            print "Searching..."
            for b in bundle:
                print bundleNum,b
                if type(b["name"]) is list:
                    print "b is a list"
                    for n in b["name"]:
                        print "\t",n
                        if nameOfBundleToDelete == n:
                            print "Deleting!"
                            del bundle[bundleNum]
                            break
                else:
                    print "not a list"
                    print "\t",b["name"]
                    if nameOfBundleToDelete == b["name"]:
                        print "Deleting!"
                        del bundle[bundleNum]
                        break
                bundleNum+=1   
        elif action == "diffMt":
            break
        
    print "-> Done comparing [0] and [1] bundles"
    #now repeat, but with 2nd allele
    print "-> Will now compare [1] bundles with [2] bundles iteratively"
    while len(bundle)>2:
        action = determine_action(bundle[1],bundle[2])
        print "Action to carry out : %s"%(action)
        if action == "sep":
            break
        elif action == "merge":
            bundle[1] = merge_bundles(bundle[1],bundle[2])
            del bundle[2]
        elif "del" in action:
            nameOfBundleToDelete = action.split('\t')[1]
            bundleNum = 0
            print "Going to delete %s"%nameOfBundleToDelete
            print "Searching..."
            for b in bundle:
                print bundleNum,b
                if type(b["name"]) is list:
                    print "b is a list"
                    for n in b["name"]:
                        print "\t",n
                        if nameOfBundleToDelete == n:
                            print "Deleting!"
                            del bundle[bundleNum]
                            break
                else:
                    print "not a list"
                    print "\t",b["name"]
                    if nameOfBundleToDelete == b["name"]:
                        print "Deleting!"
                        del bundle[bundleNum]
                        break
                bundleNum+=1   
#            print "Deleting bundle[1]: %s"%(get_name(bundle[1]))
#            del bundle[1]
#        elif action == "2del":
#            print "Deleting bundle[2]: %s"%(get_name(bundle[2]))
#            del bundle[2]
        elif action == "diffMt":
            break
    print "-> Done comparing [1] and [2] bundles"            
    print "FINAL BUNDLES:"
    for b in bundle:
        print b
        
    display_final_results(bundle, finalBundle.wellId)
        

def display_final_results(bundle, wellId):
    topNames = list()
    for bun in bundle:
        name = ""
        if type(bun["name"]) is list:
            for n in bun["name"]:
                name = "%s\%s"%(name,n)
            name = name[1:]
        else:
            name = bun["name"]
        topNames.append(name)
        
    top2 = 0
    resultStr = "<><>\t%s"%(wellId)
    for tn in topNames:
        resultStr = "%s\t%s"%(resultStr,tn)
        top2+=1
        if top2 ==2:
            break
    print resultStr
    
    
def determine_action(bunA, bunB):
    print "## Determining action for %s vs %s"%(get_name(bunA),get_name(bunB))
    #get major types, names
    majorTypeA = get_major_type(bunA)
    majorTypeB = get_major_type(bunB)
    nameA = get_name(bunA)
    nameB = get_name(bunB)
    countsA = get_counts(bunA)
    countsB = get_counts(bunB)
    
    ##########################################################
#    print "0) DEBUG: Do 1st and 2nd overlap with each other?"
#    #Reanalyse both bundles to determine if one dominates the other, i.e. get overlap score
#    openList = [nameA,nameB]
#    overlapScores = get_overlap_score(openList)
#    #if scores dependent of each other, i.e. different scores
#    if scores_similar(overlapScores)== False:
#        if overlapScores[0]["oScore"] > overlapScores[1]["oScore"]:
#            print "Del2 : Recommend removing %s"% overlapScores[1]["name"]
#            print "Return phrase is del\t%s"%overlapScores[1]["name"]
#            return "del\t%s"%overlapScores[1]["name"]
#        else:
#            print "Del1: Recommend removing %s"% overlapScores[0]["name"]
#            print "Return phrase is del\t%s"%overlapScores[0]["name"]
#            return "del\t%s"%overlapScores[0]["name"]
          
    ######################
    
    
    
    print "1) Same major type? %s vs %s"%(majorTypeA, majorTypeB)
    #if same major type
    if majorTypeA == majorTypeB:
        print "\tYes"
        #Reanalyse both bundles to determine if one dominates the other, i.e. get overlap score
        openList = [nameA,nameB]
        overlapScores = get_overlap_score(openList)
         #if scores independent of each other, i.e. similar scores:
        if scores_similar(overlapScores)== True:
            #if scores both very low:
            if is_high_overlap(overlapScores) == True:
                print "Recommend merging"
                return "merge"
            #if scores both very high, i.e. not a lot of cross hybing
            else:
                #if read count of second is a lot lower than that of first (e.g. a quarter?) discard it
                print "Assessing counts %s vs %s"%(countsA, countsB)
                if (sum(countsA)*1.0) / (sum(countsB)*1.0) > FOLD_DIFF_IN_COUNTS:
                    print "Too big a disparity in counts - chuck away %s"%nameB
                    print "Return statement: del\t%s"%nameB
                    return "del\t%s"%nameB
                else:
                    #keep both seperate
                    print "Recommend keeping seperate"
                    return "sep"
                    
                
                #merge both
        #if one dominates the other:
        else:
            #remove dominated one!
            if overlapScores[0]["oScore"] > overlapScores[1]["oScore"]:
                print "Del2 : Recommend removing %s"% overlapScores[1]["name"]
                print "Return phrase is del\t%s"%overlapScores[1]["name"]
                return "del\t%s"%overlapScores[1]["name"]
            else:
                print "Del1: Recommend removing %s"% overlapScores[0]["name"]
                print "Return phrase is del\t%s"%overlapScores[0]["name"]
                return "del\t%s"%overlapScores[0]["name"]
            
        
        
    
    #if different major type
    else:
        print "\tNo"
        return "diffMt"
      
def scores_similar(overlapScores):
    print "Are scores similar? (%.5f,%.5f)"%(overlapScores[0]["oScore"],overlapScores[1]["oScore"])
    if overlapScores[0]["oScore"] / overlapScores[1]["oScore"] > PROP_DIFF_IN_OVERLAP:
        print "No"
        return False
    elif overlapScores[1]["oScore"] / overlapScores[0]["oScore"] > PROP_DIFF_IN_OVERLAP:
        print "No"
        return False
    else:
        print "Yes"
        return True

            
def is_high_overlap(overlapScores): 
    print "Is there a lot of cross-hyb/overlap?(%.5f,%.5f)"%(overlapScores[0]["oScore"],overlapScores[1]["oScore"])
    if overlapScores[0]["oScore"] < 0.25 and overlapScores[1]["oScore"] < 0.25:
        print "Yes"
        return True
    else:
        print "No"
        return False
    
    
def get_major_type(bundle):
    if type(bundle["majorType"]) is list:
        return bundle["majorType"][0]
    else:
        return bundle["majorType"]
    
def get_counts(bundle):
    if type(bundle["countsPerAmp"][0]) is list:
        return bundle["countsPerAmp"][0]
    else:
        return bundle["countsPerAmp"]
    
def get_name(bundle):
    if type(bundle["name"]) is list:
        return bundle["name"][0]
    else:
        return bundle["name"]


    #check major type
def calc_score_ratios(finalBundle):
    scores= []
    names = []
    majorTypes = []
    countsPerAmp = []
    for entry in sorted(finalBundle.results.iteritems(), key=lambda (x,y): float(x), reverse=True):
        score = entry[0]
        result = entry[1]
        scores.append(score)
        names.append(result.targetName)
        majorTypes.append(result.targetName.split(':')[0])
        countsPerAmp.append(result.rawCounts)
        print score,result.targetName,result.targetName.split(':')[0]
        
        
    
    bundle = []
    for na,sc,mt,cpa in zip(names,scores,majorTypes,countsPerAmp):
        bundle.append({"name" : na, "score" : sc, "majorType" : mt, "countsPerAmp" : cpa})
    
    print bundle
    totScore = scores[0] + scores[1]
    
    #compare bundles - this step merges ambiguous calls together
    # e.g. if top scores are 0701 (100), 0717 (98), 0404 ( 67), 0456 ( 64)
    # it will merge to generate 0701/0717, 0404/0456 as top scores
    # it does NOT take into account difference between 0701/0717 scores and 0404/0456 scores to determine zygosity - it does that further down 
    
    while len(bundle)>1 and should_be_merged(bundle[0],bundle[1]):
        bundle[0]= merge_bundles(bundle[0],bundle[1])
        del bundle[1]
        
        
    while len(bundle)>2 and should_be_merged(bundle[1],bundle[2]):
        bundle[1]= merge_bundles(bundle[1],bundle[2])
        del bundle[2]
        
    print bundle

    #check if homozygous
    if is_homozygous(bundle):
        bundle = make_homozygous(bundle)
    
    
    print "FINAL BUNDLE: ",bundle
    
    #print "&& %s"%(finalBundle.wellId)
    display_results(bundle, finalBundle.wellId)
    
    #check if homozygous    
    if majorTypes[0]==majorTypes[1]:
        print ">R> WARNING: HOMOZYGOUS"
    
    #check if amps missing
    top2 = 0
    for c,n in zip(countsPerAmp,names):
        
        numMissingAmps = 0
        missingAmps = set()
        counter=0
        for a in c:
            if a==0:
                numMissingAmps +=1
                missingAmps.add(align.Aligner.Order[counter]) # grab amplicon id
                
            counter+=1
        #determine type of warning
        warningMsg = ""
        if numMissingAmps ==1:
            warningMsg = "Warning lvl 1: Missing 1 PCR"
        elif numMissingAmps ==2:
            print ">R> Missing",missingAmps
            if "2_4_F" in missingAmps and "2_5_R" in missingAmps:
                warningMsg = "Warning lvl 3: Missing 2 PCRs from same Amp. i.e. Half info missing!"
            elif "2_1_F" in missingAmps and "2_4_R" in missingAmps:
                warningMsg = "Warning lvl 3: Missing 2 PCRs from same Amp. i.e. Half info missing!"
            else:
                warningMsg = "Warning lvl 2: Missing 2 PCRs but from different amps"
        elif numMissingAmps == 3:
            warningMsg = "Warning lvl 4: Missing 3 PCRs. Total failure!"
            
            
        
            
            

        print ">R>",n,"\t",numMissingAmps,"\t",warningMsg

        


        
def should_be_merged(bunA,bunB):
    mtA = ""
    mtB = ""
    nameA = ""
    nameB = ""
    nameAGroup = None
    nameBGroup = None
    countsA = 0
    countsB = 0
    scoreA = 0.0
    scoreB = 0.0
    if type(bunA["majorType"]) is list:
        mtA = bunA["majorType"][0]
        countsA = bunA["countsPerAmp"][0] 
        scoreA = bunA["score"][0]
        nameA = bunA["name"][0]
    else:
        mtA = bunA["majorType"]
        countsA = bunA["countsPerAmp"]
        scoreA = bunA["score"]
        nameA = bunA["name"]
    nameAGroup = bunA["name"]
    
    if type(bunB["majorType"]) is list:
        mtB = bunB["majorType"][0]
        countsB = bunB["countsPerAmp"][0]
        scoreB = bunB["score"][0]
        nameB = bunB["name"][0]
    else:
        mtB = bunB["majorType"]
        countsB = bunB["countsPerAmp"]
        scoreB = bunB["score"]
        nameB = bunB["name"]
    nameBGroup = bunB["name"]
        
    print "&&& Should these be merged? %s vs %s"%(nameAGroup, nameBGroup)
    openList = list()
    openList.append(nameA)
    openList.append(nameB)
    
    #really, only merge if major types are the same
    if mtA == mtB:
        #if, when cross-aligned, one dominates over the other it's best not to merge but rather ditch the lower one
        overlapScores = get_overlap_score(openList)
        print "Overlap Scores:"
        for o in overlapScores.items():
            print "%s\t%.5f"%(o[0],o[1])
        
        #verify only 2 results obtained - debugging
        if len(overlapScores) !=2:
            print "^^^ ERROR: Too many overlap Scores"
            
        #convert to lists
        overlapScoreLists = []
        for o in overlapScores.items():
            overlapScoreLists.append({"name":o[0],"oScore":o[1]})
            
        if overlapScoreLists[0]["oScore"] / overlapScoreLists[1]["oScore"] > 10:
            print "Recommend ditching",overlapScoreLists[1]["name"]
        elif overlapScoreLists[1]["oScore"] / overlapScoreLists[0]["oScore"] > 10:
            print "Recommend ditching",overlapScoreLists[0]["name"]
        else:
            print "Keep both"
        #Check score disparity here ^^^^ i.e. not what score is returned but the modifiers, e.g. 0404 score goes from 1000 to 200 when 0408 reads are removed. This is a 0.2 modifier. Whereas 0408 score goes from 2000 to 1000 which is only a 0.5 modifier

#        totNewScore = 0
#        for r in result:
#            totNewScore+=r["newScore"]
#        
#        newRatio =  result[0]["newScore"]/ (totNewScore *1.0)
#        #if within tolerance
#        print "New ratio",newRatio
#        if newRatio > 0.3 and newRatio < 0.7:
#            print "Ratios are within tolerance (%.2f) - returning true - merger!"%newRatio
#            return True
#    
#    print "Returning false - no merge"
    
    return False
            
    
    
#    print "RESULT: ",result
#    print "Calculating ratios"
#    print sum(countsA)," / (",sum(countsA)," + ",sum(countsB),")"
#    ratioAB = (sum(countsA) *1.0) / (sum(countsB)*1.0)
#    totScore = scoreA + scoreB
#    scoreRatio = scoreA / (totScore*1.0)
#    
#    print ratioAB
#    if mtA == mtB:
#        if ratioAB >= 0.75 and ratioAB <= 1.33:
#            if scoreRatio > 0.4 and scoreRatio < 0.6:
#                print "Returning true"
#                return True
#    
#    print "Returning false"
#    return False
    
    
def merge_bundles(bunA, bunB):
    newBun = {}
    newBun["name"] = list()
    newBun["score"] = list()
    newBun["majorType"] = list()
    newBun["countsPerAmp"] = list()
    
    if type(bunA["name"]) is list:
        for n,s,m,c in zip(bunA["name"],bunA["score"],bunA["majorType"],bunA["countsPerAmp"]):
            newBun["name"].append(n)
            newBun["score"].append(s)
            newBun["majorType"].append(m)
            newBun["countsPerAmp"].append(c)
    else:
        newBun["name"].append(bunA["name"])
        newBun["score"].append(bunA["score"])
        newBun["majorType"].append(bunA["majorType"])
        newBun["countsPerAmp"].append(bunA["countsPerAmp"])
        
        
        if type(bunB["name"]) is list:
            for n,s,m,c in zip(bunB["name"],bunB["score"],bunB["majorType"],bunB["countsPerAmp"]):
                newBun["name"].append(n)
                newBun["score"].append(s)
                newBun["majorType"].append(m)
                newBun["countsPerAmp"].append(c)
        else:
            newBun["name"].append(bunB["name"])
            newBun["score"].append(bunB["score"])
            newBun["majorType"].append(bunB["majorType"])
            newBun["countsPerAmp"].append(bunB["countsPerAmp"])
    
    
    print "Merged into", newBun
    return newBun
    
def is_homozygous(bundle):
    print "Checking Zygosity"
    aveScoreA = 0.0
    aveScoreB = 0.0
    majorA = ""
    majorB = ""
    if len(bundle)==1:
        print "Only one entry -> homoyzygous"
        return True
    
    bunA = bundle[0]
    bunB = bundle[1]
    
    if type(bunA["score"]) is list:
        totScoreA=0.0
        majorA = bunA["majorType"][0]
        for s in bunA["score"]:
            totScoreA += s
        aveScoreA = totScoreA / (len(bunA["score"])*1.0)
    else:
        aveScoreA = bunA["score"]
        majorA = bunA["majorType"]
    if type(bunB["score"]) is list:
        totScoreB=0.0
        majorB = bunB["majorType"][0]
        for s in bunB["score"]:
            totScoreB += s
        aveScoreB = totScoreB / (len(bunB["score"])*1.0)
    else:
        aveScoreB = bunB["score"]
        majorB = bunB["majorType"]
        
    print "aveScoreA",aveScoreA
    print "aveScoreB",aveScoreB
    
    
    #check ratio
    totScoreAB = aveScoreA + aveScoreB
    ratio = aveScoreA / totScoreAB
    print "Ratio",ratio
    if ratio >= 0.3 and ratio <= 0.7:
        print " -> Potential Het"
        #check that second on list isn't same as first on list in terms of major group
        print " -> Comparing major types: %s vs %s"%(majorA, majorB)
        if majorA == majorB: 
            print "-> Same type: Likely Hom"

            return True
        else:
            print "-> Different type: Likely Het"
            return False
    else:
        print " -> Hom"
        return True
    
def make_homozygous(bundle):
    #when making homozygous, pick the one with the highest number of counts
    print "\nMAKING HOM\n"
    
    if len(bundle)==1:
        return bundle
    
    bunA = bundle[0]
    bunB = bundle[1]
    scoreA = 0.0
    scoreB = 0.0
    nameA = ""
    nameB =""
    if type(bunA["majorType"]) is list:
        iter = 0
        topCounts = 0
        bestIter = -1
        for c in bunA["countsPerAmp"]:
            if sum(c) > topCounts:
                bestIter = iter
                topCounts = sum(c)
            iter+=1
        countsA = bunA["countsPerAmp"][bestIter]
        nameA = bunA["name"][bestIter]
        
        iter = 0
        topScore = 0
        bestIter = -1
        for s in bunA["score"]:
            if s > topScore:
                bestIter = iter
                topScore = s
            iter+=1 
        scoreA = bunA["score"][bestIter]
    else:
        countsA = bunA["countsPerAmp"]
        scoreA = bunA["score"]
        nameA = bunA["name"]
        
        
    if type(bunB["majorType"]) is list:
        iter = 0
        topCounts = 0
        bestIter = -1
        for c in bunB["countsPerAmp"]:
            if sum(c) > topCounts:
                bestIter = iter
                topCounts = sum(c)
            iter+=1
        countsB = bunB["countsPerAmp"][bestIter]
        nameB = bunB["name"][bestIter]
        
        iter = 0
        topScore = 0
        bestIter = -1
        for s in bunB["score"]:
            if s > topScore:
                bestIter = iter
                topScore = s
            iter+=1 
        scoreB = bunB["score"][bestIter]
    else:
        countsB = bunB["countsPerAmp"]
        scoreB = bunB["score"]
        nameB = bunB["name"]
    newBundle = []    
    
    
    print 
    print
    print
    print "%%%%%%%%%%%%%%%%%%%%"
    print "doing mini validation"
    
    openList = list()

    openList.append(nameA)
    openList.append(nameB)
    get_overlap_score(openList)
    

    if sum(countsA) > sum(countsB):    
        newBundle.append(bundle[0])
    else:
        newBundle.append(bundle[1])
    return newBundle



def display_results(bundle, wellId):
    
    names=list()
    names.append("")
    names.append("")
    
    top2=0
    for bun in bundle:
        name = ""
        if type(bun["name"]) is list:
            for n in bun["name"]:
                name = "%s\%s"%(name,n)
            name = name[1:]
        else:
            name = bun["name"]
        names[top2]=name
        #print "&&\t%d\t%s"%(top2+1,name)
        top2+=1
        if top2 >1:
            break
        
    print "%s\t%s\t%s"%(wellId, names[0],names[1])
        
    
            
        