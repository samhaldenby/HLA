import sys
#sys.path.append('/home/sh695/Documents/Scripts/PythonScripts/newHla')
import logging
import hla.cluster as cluster
import hla.reference as reference
import hla.align as align
from hla.top_pickers import *
from hla.rankers import *
import hla.primer_info as primer_info
from hla.validation import *
import hla.final_results as final_results
from hla.allele_status import Status
import hla.contaminants as contam
from optparse import OptionParser

import copy
import os




##Info on how reads should be trimmed
#class PrimerInfo:
#    name = ""
#    trimFrom = None
#    trimTo = None
#    
#    def __init__(self,name,trimFrom,trimTo):
#        self.name = name
#        self.trimFrom = trimFrom
#        self.trimTo = trimTo
        





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

    # 1 - grab fastq file names
    logging.info('Grabbing file names from %s',wellPath)
    fqFileNames = []
    for fqFileName in os.listdir(wellPath):
        if fqFileName.endswith(".fq"):# and "2_1" in fqFileName and "2_3" not in fqFileName:
            fqFileNames.append(fqFileName)
            logging.info(" -> %s",fqFileName)
    
    clusts = {}
    #for each fastq file, prepare clusters
    for fqFileName in fqFileNames:
        #determine what primer was used
        
        for primerInfo in primerMap:
            logging.info('Looking for %s in %s',primerInfo,fqFileName)
            if primerInfo in fqFileName:
                logging.info('found!')
                clusts[primerInfo] = cluster.DnaClusters()
                currCluster = clusts[primerInfo]
                info = primerMap[primerInfo]
                currCluster.set_read_region(info.trimFrom, info.trimTo)
                currCluster.load_from_fq("%s%s"%(wellPath,fqFileName));
            else:
                logging.info('not found!!!')
                
        
    
    
    #export clusters
    minClusterSizeForExport = 2
    for entry in clusts.items():
        id = entry[0]
        clust = entry[1]
        clust.export_clusters("%s_%s_clusters_%s.fa"%(wellId,id,minClusterSizeForExport),minClusterSizeForExport);



    #create aligner and do initial alignment
    a = align.DnaAligner()
    for entry in clusts.items():
        id = entry[0]
        clust = entry[1]
        a.align(clust, r, unmappedFileName = "NM_%s_%s.fa"%(wellId, id));
    
    a.compile_results()
    
   
    
    
    #determine contaminant levels
    print ">>",wellId
    a.report_basic_stats()
    contamInfo = contam.ContaminantInfo(a.get_results())
    contamProp = contamInfo.calculate_contamination_levels()
    realResults = contamInfo.get_real_mapped_counts()
    print ">>______________"
    
     #print raw results 
    print "* Raw Counts"
    rawTopHits = pick_Nx(a.get_results(),nX)
    rank_results(rawTopHits)
    
    #copy raw and modify 
    rawResults = copy.deepcopy(a.get_results())

    
    #grab top results
    print "* Running 1st alignment"
    topHits = pick_Nx(a.get_results(),nX)
    
    
    #prepare list of results
    #results =[]
    #rank them
    #results.append(rank_results(topHits))
    #allRankedResults = [] #will be a list of all results after all top pickers have been selected and reanalysis is complete
    allRankedResults = {} #{omitted target : Results}
    allRankedResults["First"] =rank_results(topHits)

    openList = []
    topHitNames = set()
    
    #grab top 2 hits and add to open list
    for entry in sorted(allRankedResults["First"], key=allRankedResults["First"].get, reverse=True):
        name = entry
        score = allRankedResults["First"][entry]
        print name, score
        if len(topHitNames)<2:
            openList.append(name)
            topHitNames.add(name)
        else:
            break
        
    print "Open list after 1st run: %s"%openList
    
    
 #   allSepResults = {} # {omitted target : [0,122,12314,124]}
 #   allSepResults["First"]= firstResults
        
        
    
    #realign, omiting r1 top hit
    print "\n* Running 2nd alignment, omitting reads mapping to top scorer from 1st (%s)"%openList[0]
    
    for entry in clusts.items():
        id = entry[0]
        clust = entry[1]
        a.align(clust, r, openList[0])

    a.compile_results()
    topHits = pick_Nx(a.get_results(),nX)
    #topHits = normalise_hits(topHits,realResults)
    #results.append(rank_results(topHits))
    allRankedResults[openList[0]] =rank_results(topHits)
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(allRankedResults[openList[0]], key=allRankedResults[openList[0]].get, reverse=True):
        name = entry
        score = allRankedResults[openList[0]][entry]
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
    for entry in clusts.items():
        id = entry[0]
        clust = entry[1]
        a.align(clust, r, openList[1])
    
#    a.align(c1,r, openList[1])
#    a.align(c2,r, openList[1])
#    a.align(c3,r, openList[1])
#    a.align(c4,r, openList[1])
    a.compile_results()
#    a.modify_scores_for_read_imbalances()
    topHits = pick_Nx(a.get_results(),nX)
    #topHits = normalise_hits(topHits,realResults)
    #results.append(rank_results(topHits))
    allRankedResults[openList[1]] =rank_results(topHits)
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(allRankedResults[openList[1]], key=allRankedResults[openList[1]].get, reverse=True):
        name = entry
        score = allRankedResults[openList[1]][entry]
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
        for entry in clusts.items():
            id = entry[0]
            clust = entry[1]
            a.align(clust, r, openList[len(openList)-remainingToCheck])
#        a.align(c1,r, openList[len(openList)-remainingToCheck])
#        a.align(c2,r, openList[len(openList)-remainingToCheck])
#        a.align(c3,r, openList[len(openList)-remainingToCheck])
#        a.align(c4,r, openList[len(openList)-remainingToCheck])
        a.compile_results()
#        a.modify_scores_for_read_imbalances()
        topHits = pick_Nx(a.get_results(),nX)
      #  topHits = normalise_hits(topHits,realResults)
        #results.append(rank_results(topHits))
        allRankedResults[openList[len(openList)-remainingToCheck]] =rank_results(topHits)
        
        #now see if there is a new top hit not already on open list
        
        print "Open list after run: %s"%openList
        
        remainingToCheck -= 1
        
   
        
    print ">> *** %s ***"%wellId
    resultsBundle = process_all_results2(rawResults, allRankedResults, openList)
    resultsBundle.wellId = wellId
    
    #prepare output file for result line
    outFile = open(outputTag,"a")
    outFile.write("%s\n"%resultsBundle.result_string())
    outFile.close()
        
    #comparing to dynal results? 
    if dynalResults != None:
        #compare_with_dynal(resultsBundle, wellId, dynalResults)
        compare_with_cam2(resultsBundle, wellId, dynalResults)

    #check ratios
    calc_score_ratios(resultsBundle)
   
    print ">> _______________"
        
        
        
        
        
        
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
        



        
        
        

#def process_all_results(rawResults, results, openList): #results is {omitedGene: {name:score}}, openList is [name]
#    print
#    print "////////////////Processing results////////////////"
#    print "OpenList: ",openList
#    print
#    print "results: ", results
#    print
#    #create map
#    #scores = {}
#    originalScores = {}
#    rejectionList = set()
#    for gene in openList: #string in []
#        #scores[gene]=[] 
#        resultsForThisGene = results[gene]
#        
#        for target in openList:  #iterate again to do comparisons
#            if target !=gene:
#                if target not in resultsForThisGene:
#                    print target,"not in",resultsForThisGene
#                    rejectionList.add(target)
#                elif resultsForThisGene[target]==0.0:
#                    print target,"has zero score with",gene
#                    rejectionList.add(target)
#                else:
#                    print target,"has score",resultsForThisGene[target],"with",gene
#    
#    #now create accepted list
#    newOpenList = []
#    for gene in openList:
#        if gene not in rejectionList:
#            newOpenList.append(gene)
#            
#    openList = newOpenList
#                
#                
#        for result in results:
#            if gene in result:  #i.e. if it's not the one being omited, OR!!! if it got no hits.... deal with this later TODO:
#                scores[gene].append(result[gene])
#                if gene not in originalScores:
#                    originalScores[gene]= result[gene]
#
#   
#    
#    #now, make sure all modifier lists are same length (reason they might not be is due to getting 0 hits)
#    longestList = 0
#    for s in scores.items():
#        if len(s[1]) > longestList:
#            longestList = len(s[1])
#        
#    for s in scores.items():
#        while len(s[1]) < longestList:
#            s[1].append(0.0)
#
#    modifiers = {}
#    for entry in scores.items():
#
#        gene = entry[0]
#        modifiers[gene]=[]
#
#        for result in entry[1]:
#            modifiers[gene].append(result/originalScores[gene])
#
#    #multiply all modifiers by all scores
#    modifiedScores = {}
#    unmodifiedScores = {}
#    #for each gene
#    for entry in scores.items():
#        thisScore = 1.0
#        unmodifiedThisScore = 1.0
#        geneName = entry[0]
#        geneScores = entry[1]
#        #multiply all scores
#        for s in geneScores:
#            thisScore *= s
#            if s>1:
#                unmodifiedThisScore *=s
#        #mow multiply by modifiers
#        
#
#        for m in modifiers[geneName]:
#            thisScore *= m
#        
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
#    print "***************"
#    for gene in scores.items():
#        modifiersForThisOne = modifiers[gene[0]]
#        print "%s\t%s\t%s"%(gene[0],'\t'.join(map(str,gene[1])), '\t'.join(map(str,modifiersForThisOne)))
#    print "***************"   
#    
#    prevScore = 0
#    scoreNum=0
#    finalResults = {} # list of FinalResult class instances
#    for entry in sorted(modifiedScores, key=modifiedScores.get, reverse=True):
#        scoreNum +=1
#        thisScore = modifiedScores[entry]
#
#        
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
#        #print "RAWRESULTS:",rawResults
#        #raw_input()
#        #print "RESULTS:",results
#        #print ">>%s\t%s\t%s"%(entry, modifiedScores[entry], status) #unmodifiedScores[entry],status
#        f = final_results.FinalResult()
#        f.score = modifiedScores[entry]
#        f.targetName = entry
#        f.status = status
#        f.rawCounts = rawResults[entry]
#        finalResults[f.score]=(f)
#        
#                #if score is 0, re-run but omit from open list
#        if thisScore == 0:
#            print "ZERO SCORE FOUND - re-analysing"
#            openList.remove(entry)
#            #results.remove(entry)
#            #print "Need to remove %s from raw results which are %s"%(entry,rawResults)
#            #del rawResults[entry]
#            resultsBundle = process_all_results(rawResults, results, openList)
#            return resultsBundle
#        
#        print ">>%s\t%s\t%s\t%s"%(entry,modifiedScores[entry],status,rawResults[entry])
#        prevScore = thisScore
#        
#    
#    resultsBundle = final_results.FinalResultBundle()
#    resultsBundle.results = finalResults
#    print "FINAL RESULTS: %s"%(finalResults)
#            
#    return resultsBundle


    
    
        

    
   
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
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-R", "--reference", dest="referenceName", action="store", help="tab-separated (name seq) reference file name", metavar="FILE")
    parser.add_option("-d", "--inputPath", dest="inputPath", action="store", help="path to directory containing well directories, e.g. /data/A04/MiSeq/", metavar="PATH")
    parser.add_option("-w", "--well", dest="wellId", action="store", help="well ID (e.g. D04)", metavar="STRING")
    parser.add_option("-L", "--logTag", dest="logTag", help="log file tag", metavar="STRING")
    parser.add_option("-o", "--output", dest="outputTag", help="output file tag", metavar="STRING")
    parser.add_option("-r", "--dynal-results", dest="dynalResults", help="Previous dynal results to compare results to", metavar="FILE")
    (options,args) = parser.parse_args()
    
    
    
    

    
    primerMap = {}
    #primerMap["2_2_F"] = primer_info.PrimerInfo("2_2_F",0,160)
    primerMap["2_4_F"] = primer_info.PrimerInfo("2_4_F",1,150) 
    primerMap["2_4_R"] = primer_info.PrimerInfo("2_4_R",50,190) #50,190
    primerMap["2_5_R"] = primer_info.PrimerInfo("2_5_R",70,200)
    primerMap["2_1_F"] = primer_info.PrimerInfo("2_1_F",70,180) #80,180
    
    
    run_analysis(referenceName=options.referenceName, inputPath=options.inputPath, wellId=options.wellId, logTag=options.logTag, outputTag=options.outputTag, dynalResults=options.dynalResults)

    


