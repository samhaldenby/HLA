import sys
#sys.path.append('/home/sh695/Documents/Scripts/PythonScripts/newHla')
import logging
import hla.cluster as cluster
import hla.reference as reference
import hla.align as align
from hla.top_pickers import *
from hla.rankers import *
import hla.primer_info as primer_info
import hla.all_ranked_results as AllRankedResults
from hla.validation import *
from hla.process_results import *


import hla.contaminants as contam
from optparse import OptionParser

import copy
import os



    
primerRunningOrder = []


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
                primerRunningOrder.append(primerInfo)
                logging.info('found!')
                clusts[primerInfo] = cluster.DnaClusters()
                currCluster = clusts[primerInfo]
                info = primerMap[primerInfo]
                print ">>SETTING READ REGION"
                currCluster.set_read_region(info.trimFrom, info.trimTo)
                print "LOADING FROM FQ"
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
    
    AllRankedResults.results = allRankedResults
    
    #prepare output file for result line
    outFile = open(outputTag,"a")
    outFile.write("%s\n"%resultsBundle.result_string())
    outFile.close()
        
    #comparing to dynal results? 
    if dynalResults != None:
        #compare_with_dynal(resultsBundle, wellId, dynalResults)
        compare_with_cam2(resultsBundle, wellId, dynalResults)

    #check ratios
    polish_results(resultsBundle)
    #calc_score_ratios(resultsBundle)
   
    print ">> _______________"
        
        
        
        
        
        

        



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
    primerMap["2_2_F"] = primer_info.PrimerInfo("2_2_F",0,160)#0,160
    primerMap["2_4_F"] = primer_info.PrimerInfo("2_4_F",1,150)#1,150
    primerMap["2_4_R"] = primer_info.PrimerInfo("2_4_R",50,190) #50,190
    primerMap["2_5_R"] = primer_info.PrimerInfo("2_5_R",70,210) #70,200
    #primerMap["2_1_F"] = primer_info.PrimerInfo("2_1_F",-70,-187) #70,180
    #primerMap["2_1_F"] = primer_info.PrimerInfo("2_1_F",70,180) #70,180
    
    run_analysis(referenceName=options.referenceName, inputPath=options.inputPath, wellId=options.wellId, logTag=options.logTag, outputTag=options.outputTag, dynalResults=options.dynalResults)
    #print "DRUM ROLL.....",align.Aligner.Order
    

    #print "PRIMER RUNNING ORDER",primerRunningOrder
