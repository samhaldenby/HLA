import logging
import hla.align as align
'''
Created on 21 May 2012

@author: sh695
'''

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
    
    countsA = 0
    countsB = 0
    scoreA = 0.0
    scoreB = 0.0
    if type(bunA["majorType"]) is list:
        mtA = bunA["majorType"][0]
        countsA = bunA["countsPerAmp"][0] 
        scoreA = bunA["score"][0]
    else:
        mtA = bunA["majorType"]
        countsA = bunA["countsPerAmp"]
        scoreA = bunA["score"]
    if type(bunB["majorType"]) is list:
        mtB = bunB["majorType"][0]
        countsB = bunB["countsPerAmp"][0]
        scoreB = bunB["score"][0] 
    else:
        mtB = bunB["majorType"]
        countsB = bunB["countsPerAmp"]
        scoreB = bunB["score"]
        
    print "&&& Should these be merged? %s vs %s"%(mtA, mtB)
    
    print "Calculating ratios"
    print sum(countsA)," / (",sum(countsA)," + ",sum(countsB),")"
    ratioAB = (sum(countsA) *1.0) / (sum(countsB)*1.0)
    totScore = scoreA + scoreB
    scoreRatio = scoreA / (totScore*1.0)
    
    print ratioAB
    if mtA == mtB:
        if ratioAB >= 0.75 and ratioAB <= 1.33:
            if scoreRatio > 0.4 and scoreRatio < 0.6:
                print "Returning true"
                return True
    
    print "Returning false"
    return False
    
    
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
    
    newBundle = []    
    
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
        
    
            
        