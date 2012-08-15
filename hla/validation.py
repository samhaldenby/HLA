import logging
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
    top2 =0
    scores= []
    for entry in sorted(finalBundle.results.iteritems(), key=lambda (x,y): float(x), reverse=True):
        score = entry[0]
        result = entry[1]
        scores.append(score)
        print score,result
        
        
        top2 +=1
        if top2==2:
            break
        
    totScore = sum(scores)
    for score in scores:
        print ">R>>\t",score/totScore
#        self.resultsSortedList.append(result)
#    for result in finalBundle.results:
#        print "ASDFASFSV: ",result
        
        
    
    
        
        