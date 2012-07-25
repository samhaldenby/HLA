import logging
'''
Created on 21 May 2012

@author: sh695
'''


def compare_with_dynal(resultsBundle, wellId, dynalFileName):
    
    logging.info('Comparing to dynal results')
                 
    print "*** %s ***"%wellId
                 
    #open dynal results file
    try:
        dynalFile = open(dynalFileName)
    except IOError:
        logging.error("Unable to open %s",dynalFileName)
    
    
    #grab correct line from dynal results
    dynalResultsLine = ""
    for line in dynalFile:
        if wellId in line:
            dynalResultsLine = line.strip()
            
    dynalResultsLineSplit = dynalResultsLine.split("\t")
    firstDynalHit = dynalResultsLineSplit[1]
    secondDynalHit = dynalResultsLineSplit[2]
    
    #compile dynal results
    dynalResults = set()
    dynalResults.add(firstDynalHit)
    dynalResults.add(secondDynalHit)
    print ">>",dynalResults
    
    
    score = 0
    for allele in resultsBundle.accepted_alleles:
        if allele in dynalResults:
            score+=1
            print ">> %s\tYes"%allele
        else:
            print ">> %s\tNo"%allele
    
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
        
    
    
    
    
        
        