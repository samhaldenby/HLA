import logging
import hla.cluster as cluster
import hla.reference as reference
import hla.align as align
from hla.top_pickers import *
from hla.rankers import *
from hla.validation import *
from optparse import OptionParser

###NOTES
# 1 - The program has trouble in detecting 0408....is this a primer issue? Compare 04XX with 0408 and see
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
    
    #open logger
    
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p', filename="%s.log"%logTag, level=logging.DEBUG)
    logging.info('\n***** Commencing run *****')

    #prepare paths
    wellPath="%s/%s/"%(inputPath, wellId)
    
    #load references
    r = reference.Reference()
    r.load_reference_tabbed(referenceName)#("/home/sh695/Documents/Scripts/PythonScripts/HLA/modules/subject.txt")


    #create aligner
    a = align.Aligner()
    
    #load reads into clusters 
    c1 = cluster.Clusters()
    c1.set_read_region(0,160)
    c1.load_from_fq("%s2_2_F.fq"%wellPath)
    a.align(c1, r)
    #a.align_unique_only(c1, r)

    
    
    c2 = cluster.Clusters()
    c2.set_read_region(0,150)
    c2.load_from_fq("%s2_4_F.fq"%wellPath)
    a.align(c2, r)
    #a.align_unique_only(c2, r)
    
    
    c3 = cluster.Clusters()
    c3.set_read_region(33,200)
    c3.load_from_fq("%s2_4_R.fq"%wellPath)
    #a.calculate_alignment_cross_talk(c3, r)
    #raw_input("Press key")
    a.align(c3, r)
    #a.align_unique_only(c3, r)
    
    
    c4 = cluster.Clusters()
    c4.set_read_region(70,232)
    c4.load_from_fq("%s2_5_R.fq"%wellPath)
    a.align(c4, r)
    #a.align_unique_only(c4, r)

   
    a.compile_results()


    #grab top results
    topHits = pick_Nx(a.get_results(),50)   #TODO: Change this back to 50

    
    #rank them
    results = rank_by_cross_division(topHits)
    
        

    
    #comparing to dynal results? 
    if dynalResults != None:
        compare_with_dynal(results, wellId, dynalResults)


#Parse command line options
parser = OptionParser()
parser.add_option("-R", "--reference", dest="referenceName", action="store", help="tab-separated (name seq) reference file name", metavar="STRING")
parser.add_option("-d", "--inputPath", dest="inputPath", action="store", help="path to directory containing well directories", metavar="STRING")
parser.add_option("-w", "--well", dest="wellId", action="store", help=".Well ID (e.g. D04)", metavar="STRING")
parser.add_option("-L", "--logTag", dest="logTag", help="log file tag", metavar="STRING")
parser.add_option("-o", "--output", dest="outputTag", help="output file tag", metavar="STRING")
parser.add_option("-r", "--dynal-results", dest="dynalResults", help="Previous dynal results to compare results to", metavar="FILE")
(options,args) = parser.parse_args()


run_analysis(referenceName=options.referenceName, inputPath=options.inputPath, wellId=options.wellId, logTag=options.logTag, outputTag=options.outputTag, dynalResults=options.dynalResults)


