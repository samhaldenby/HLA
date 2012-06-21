import logging
import hla.cluster as cluster
import hla.reference as reference
import hla.align as align
from hla.top_pickers import *
from hla.rankers import *
from hla.validation import *
import hla.final_results as final_results
from hla.allele_status import Status
import hla.contaminants as contam
from optparse import OptionParser
import sys
import copy

###NOTES
# 1 - The program has trouble in detecting 0408/0901....is this a primer issue? Compare 04XX with 0408 and see
#Current changes: 
#    1) 1 point per read hit - no dividing if ambiguous
#    2) Everything else normalish? Don't forget iterative approach is included though!
#    3) Appears that ranking by cross division is not being carried out, i.e. not multiplying by rankings. Check this by changing N30 to N50+
#        3a) Trying this now - Done! Makes no difference!
def run_analysis(referenceName,  inputPath, wellId, logTag, outputTag, dynalResults=None):
   
    '''
    Main function for running hla analysis
    
    args:
        referenceName:    tab-separated reference file name (name \t sequence)
        inputPath:        path to directory containing well subdirectories
        logTag:           name for log file (will be appended with .log)
        outputTag:        name for output file (will be appended with relevant suffix)
        dynalResults:     name of file containing dynal results. Only needed if comparing
    
    return:
        Nothing
    '''

    #temp var
    nX = 100
    #open logger
    
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p', filename="%s.log"%logTag, level=logging.DEBUG)
    logging.info('\n***** Commencing run *****')

    #prepare paths
    wellPath="%s/%s/"%(inputPath, wellId)
    

    #load references
    r = reference.DnaReference()     #TODO: Currently, Use this even if it is an AA reference
    r.load_reference_tabbed(referenceName)#("/home/sh695/Documents/Scripts/PythonScripts/HLA/modules/subject.txt")


    #load FQs into clusters
    c1 = cluster.DnaClusters()
    c1.set_read_region(0,160) #orig[0,160]
    c1.load_from_fq("%s2_2_F.fq"%wellPath)
    c2 = cluster.DnaClusters()
    c2.set_read_region(10,150) #orig [0,150]
    c2.load_from_fq("%s2_4_F.fq"%wellPath)
    c3 = cluster.DnaClusters()
    c3.set_read_region(33,200)  #orig [33,200]
    c3.load_from_fq("%s2_4_R.fq"%wellPath)
    c4 = cluster.DnaClusters()
    c4.set_read_region(70,232) #orig [70,232]
    c4.load_from_fq("%s2_5_R.fq"%wellPath)
   
   
    #create aligner and do initial alignment
    
    a = align.DnaAligner()
    a.align(c1, r)
    a.align(c2, r)
    a.align(c3, r)
    a.align(c4, r)
    a.compile_results()
    rawResults = copy.deepcopy(a.get_results())
    a.modify_scores_for_read_imbalances()
    
    
    #determine contaminant levels
    print ">>",wellId
    a.report_basic_stats()
    contamInfo = contam.ContaminantInfo(a.get_results())
    contamProp = contamInfo.calculate_contamination_levels()
    realResults = contamInfo.get_real_mapped_counts()
    print ">>______________"

    
    #grab top results
    print "* Running 1st alignment"
    topHits = pick_Nx(a.get_results(),nX)
    
   # topHits = normalise_hits(topHits,realResults)
    #prepare list of results
    results =[]
    #rank them
    results.append(rank_results(topHits))
    

    openList = []
    topHitNames = set()
    
    #grab top 2 hits and add to open list
    for entry in sorted(results[-1], key=results[-1].get, reverse=True):
        name = entry
        score = results[-1][entry]
        print name, score
        if len(topHitNames)<2:
            openList.append(name)
            topHitNames.add(name)
        else:
            break
        
    print "Open list after 1st run: %s"%openList
        
        
    
    #realign, omiting r1 top hit
    print "\n* Running 2nd alignment, omitting reads mapping to top scorer from 1st (%s)"%openList[0]
    a.align(c1,r, openList[0])
    a.align(c2,r, openList[0])
    a.align(c3,r, openList[0])
    a.align(c4,r, openList[0])
    a.compile_results()
 #   a.modify_scores_for_read_imbalances()
    topHits = pick_Nx(a.get_results(),nX)
    #topHits = normalise_hits(topHits,realResults)
    results.append(rank_results(topHits))
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(results[-1], key=results[-1].get, reverse=True):
        name = entry
        score = results[-1][entry]
        if firstLine:
            if name not in openList:
                print "%s is new top hit. Not already in open list. Adding"%name
                openList.append(name)
                topHitNames.add(name)
            firstLine = False         
        else:
            break

    print "Open list after 2nd run: %s"%openList
        
        
        
    #realign, omiting r1 second-top hit
    print "\n* Running 3rd alignment, omitting reads mapping to second-top scorer from 1st (%s)"%openList[1]
    a.align(c1,r, openList[1])
    a.align(c2,r, openList[1])
    a.align(c3,r, openList[1])
    a.align(c4,r, openList[1])
    a.compile_results()
#    a.modify_scores_for_read_imbalances()
    topHits = pick_Nx(a.get_results(),nX)
    #topHits = normalise_hits(topHits,realResults)
    results.append(rank_results(topHits))
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(results[-1], key=results[-1].get, reverse=True):
        name = entry
        score = results[-1][entry]
        if firstLine:
           # print "WE ARE NOW CHECKING IF %s IS ON OPEN LIST %s"%(name,openList)
            if name not in openList:
                print "%s is new top hit. Not already in open list. Adding"%name
                openList.append(name)
                topHitNames.add(name)
            #print "CHECK DONE"
            firstLine = False         
        else:
            break
        
    print "Open list after 3rd run: %s"%openList
    
    
    remainingToCheck = len(openList) - 2
    finalResults =[]
    #print "remainingToCheck = %d"%remainingToCheck
    while remainingToCheck != 0:
        print "\n* Running final alignments, omitting next target from open list (%s)"%(openList[len(openList)-remainingToCheck])
        a.align(c1,r, openList[len(openList)-remainingToCheck])
        a.align(c2,r, openList[len(openList)-remainingToCheck])
        a.align(c3,r, openList[len(openList)-remainingToCheck])
        a.align(c4,r, openList[len(openList)-remainingToCheck])
        a.compile_results()
#        a.modify_scores_for_read_imbalances()
        topHits = pick_Nx(a.get_results(),nX)
      #  topHits = normalise_hits(topHits,realResults)
        results.append(rank_results(topHits))
        
        #now see if there is a new top hit not already on open list
        
        print "Open list after run: %s"%openList
        
        remainingToCheck -= 1
        
        
    print ">> *** %s ***"%wellId
    resultsBundle = process_all_results(rawResults, results, openList)
    resultsBundle.wellId = wellId
    
    #prepare output file for result line
    outFile = open(outputTag,"a")
    outFile.write("%s\n"%resultsBundle.result_string())
    outFile.close()
        
    #comparing to dynal results? 
    #if dynalResults != None:
    #    compare_with_dynal(modifiedResults, wellId, dynalResults)

   
    print ">> _______________"
        

def process_all_results(rawResults, results, openList): #results is [{name:score}], openList is [name]
    #create map
    scores = {}
    originalScores = {}
    for gene in openList:
        scores[gene]=[]
        for result in results:
            if gene in result:  #i.e. if it's not the one being omited, OR!!! if it got no hits.... deal with this later TODO:
                scores[gene].append(result[gene])
                if gene not in originalScores:
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


    
    
        

    
   
#def normalise_hits(topHits,realResults):
#    retHits = {}
#    for entry in topHits.items():
#        name = entry[0]
#        origHits = entry[1]
#        newHits = []
#        for i in range(0,len(origHits)):
#           if realResults[i]==0:
#               newHits.append(0.0)
#            else:
#                mod =  (realResults[i]*1.0)/(sum(realResults)*1.0)
#                #print "Sum(realResults): ",sum(realResults)
#                #print "realResults[i]: ",realResults[i]
#                #print "DenomMod = ",denominatorMod
#                newHits.append(((origHits[i]*100000)/(realResults[i])) *mod)
#        retHits[name] = newHits
#        
#    return retHits

#Parse command line options

parser = OptionParser()
parser.add_option("-R", "--reference", dest="referenceName", action="store", help="tab-separated (name seq) reference file name", metavar="FILE")
parser.add_option("-d", "--inputPath", dest="inputPath", action="store", help="path to directory containing well directories, e.g. /data/A04/MiSeq/", metavar="PATH")
parser.add_option("-w", "--well", dest="wellId", action="store", help="well ID (e.g. D04)", metavar="STRING")
parser.add_option("-L", "--logTag", dest="logTag", help="log file tag", metavar="STRING")
parser.add_option("-o", "--output", dest="outputTag", help="output file tag", metavar="STRING")
parser.add_option("-r", "--dynal-results", dest="dynalResults", help="Previous dynal results to compare results to", metavar="FILE")
(options,args) = parser.parse_args()

run_analysis(referenceName=options.referenceName, inputPath=options.inputPath, wellId=options.wellId, logTag=options.logTag, outputTag=options.outputTag, dynalResults=options.dynalResults)



