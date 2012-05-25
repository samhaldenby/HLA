import logging
import hla.cluster as cluster
import hla.reference as reference
import hla.align as align
from hla.top_pickers import *
from hla.rankers import *
from hla.validation import *
from optparse import OptionParser

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
    r = reference.DnaReference()
    r.load_reference_tabbed(referenceName)#("/home/sh695/Documents/Scripts/PythonScripts/HLA/modules/subject.txt")


    #load FQs into clusters
    c1 = cluster.DnaClusters()
    c1.set_read_region(0,160)
    c1.load_from_fq("%s2_2_F.fq"%wellPath)
    c2 = cluster.DnaClusters()
    c2.set_read_region(0,150)
    c2.load_from_fq("%s2_4_F.fq"%wellPath)
    c3 = cluster.DnaClusters()
    c3.set_read_region(33,200)
    c3.load_from_fq("%s2_4_R.fq"%wellPath)
    c4 = cluster.DnaClusters()
    c4.set_read_region(70,232)
    c4.load_from_fq("%s2_5_R.fq"%wellPath)
   
   
    #create aligner and do initial alignment
    a = align.DnaAligner()
    a.align(c1, r)
    a.align(c2, r)
    a.align(c3, r)
    a.align(c4, r)
    a.compile_results()
    
    #grab top results
    topHits = pick_Nx(a.get_results(),nX)   #TODO: Change this back to 50

    #rank them
    results1 = rank_by_cross_division(topHits)
    

    openList = []
    topHitNames = set()
    
    #grab top 2 hits and add to open list
    for entry in sorted(results1, key=results1.get, reverse=True):
        name = entry
        score = results1[entry]
        print name, score
        if len(topHitNames)<2:
            openList.append(name)
            topHitNames.add(name)
        else:
            break
        
    print "Open list from first run: %s"%openList
        
        
    
    #realign, omiting r1 top hit
    print "\n* Running second alignment, omiting top scorer from 1st (%s)"%openList[0]
    a.align(c1,r, openList[0])
    a.align(c2,r, openList[0])
    a.align(c3,r, openList[0])
    a.align(c4,r, openList[0])
    a.compile_results()
    topHits = pick_Nx(a.get_results(),nX)
    results2_1 = rank_by_cross_division(topHits)
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(results2_1, key=results2_1.get, reverse=True):
        name = entry
        score = results2_1[entry]
        if firstLine:
            if name not in openList:
                print "%s is new top hit. Not already in open list. Adding"%name
                openList.append(name)
                topHitNames.add(name)
            firstLine = False         
        else:
            break

    print "Open list after run2, omiting top hit: %s"%openList
        
        
        
    #realign, omiting r1 second-top hit
    print "\n* Running second alignment, omiting second-top scorer from 1st (%s)"%openList[1]
    a.align(c1,r, openList[1])
    a.align(c2,r, openList[1])
    a.align(c3,r, openList[1])
    a.align(c4,r, openList[1])
    a.compile_results()
    topHits = pick_Nx(a.get_results(),nX)
    results2_2 = rank_by_cross_division(topHits)
    #now see if there is a new top hit not already on open list
    firstLine = True
    for entry in sorted(results2_2, key=results2_2.get, reverse=True):
        name = entry
        score = results2_2[entry]
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
        
    print "Open list after run2, omiting second-top hit: %s"%openList
    
    
    remainingToCheck = len(openList) - 2
    finalResults =[]
    print "remainingToCheck = %d"%remainingToCheck
    while remainingToCheck != 0:
        print "\n* Running final alignments, omiting #%d from open list (%s)"%(len(openList)-remainingToCheck,openList[len(openList)-remainingToCheck])
        a.align(c1,r, openList[len(openList)-remainingToCheck])
        a.align(c2,r, openList[len(openList)-remainingToCheck])
        a.align(c3,r, openList[len(openList)-remainingToCheck])
        a.align(c4,r, openList[len(openList)-remainingToCheck])
        a.compile_results()
        topHits = pick_Nx(a.get_results(),nX)
        finalResults.append(rank_by_cross_division(topHits))
        #now see if there is a new top hit not already on open list
        
        print "Open list after run2, omiting #%d from open list (%s)"%(len(openList)-remainingToCheck, openList[len(openList)-remainingToCheck])
        
        remainingToCheck -= 1
        
  
        
    
        
    
    

#    #grab top hit again
#    topHitNames = set()
#    for entry in sorted(results, key=results.get, reverse=True):
#        name = entry
#        score = results[entry]
#        print name, score
#        if len(topHitNames)==0:
#            topHitNames.add(name)
#            no2Hit = name
#            break
#        
#    #realign
#    topHitNames.add(no1Hit)
#    topHitNames.add(no2Hit)
#    realigner = align.DnaAligner()
#    realigner.align(c1,r, topHitNames)
#    realigner.align(c2,r, topHitNames)
#    realigner.align(c3,r, topHitNames)
#    realigner.align(c4,r, topHitNames)
#    realigner.compile_results()
#    topHits = pick_Nx(realigner.get_results(),nX)
#    results = rank_by_cross_division(topHits)
#    
#     #grab top hit again
#    topHitNames = set()
#    for entry in sorted(results, key=results.get, reverse=True):
#        name = entry
#        score = results[entry]
#        print name, score
#        if len(topHitNames)==0:
#            topHitNames.add(name)
#          #  no1Hit = name
#            break
#        
#    print
#    print "_____________"
#    print "1: %s\n2: %s"%(no1Hit,no2Hit)

    
    
        

    
    #comparing to dynal results? 
    if dynalResults != None:
        compare_with_dynal(results, wellId, dynalResults)


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



