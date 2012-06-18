'''
Created on 29 May 2012

@author: sh695
'''

import logging

class ContaminantInfo(object):
    '''
    classdocs
    '''
    _rawAlignmentResults = {}
    #raw alignment data gets split into contaminants and non-contaminants
    _contamResults = {}
    _realResults = {}
    def __init__(self, rawResults):
        '''
        Constructor
        '''
        self._rawAlignmentResults = rawResults 
        
        
    def calculate_contamination_levels(self):
        logging.info('Determining levels of contamination')
        if len(self._rawAlignmentResults) == 0:
            logging.info('Failed to determine contamination levels as no alignments were reported')
            return
        
        #determine number of amplicons being used in this exp (e.g. if 2_2_F, 2_4_F, 2_2_R, 2_5_R: 4 amplicons
        numAmplicons = -1
        for entry in self._rawAlignmentResults.items():
            numAmplicons = len(entry[1])
            break
        
        contamTotals = [0]*numAmplicons
        realTotals = [0]*numAmplicons
        overallTotals = [0]*numAmplicons
        
        for entry in self._rawAlignmentResults.items():
            refName = entry[0]
            refScores = entry[1]
            
            #split into contam/real depending on name
            if "DRB1" in refName:
                self._realResults[refName]=refScores
                for r in range(0,numAmplicons):
                    realTotals[r]+=refScores[r]
                    
            else:
                self._contamResults[refName]=refScores
                for r in range(0,numAmplicons):
                    contamTotals[r]+=refScores[r]
                    
            for r in range(0,numAmplicons):
                overallTotals[r]+=refScores[r]
                
        contamProps = [0.0]*numAmplicons
        realProps = [0.0]*numAmplicons
        
        for r in range(0,numAmplicons):
            if overallTotals[r]==0:
                contamProps[r]=0
                realTotals[r]=0
            else:
                contamProps[r] = (contamTotals[r]*1.0)/(overallTotals[r]*1.0)
                realProps[r] = (realTotals[r]*1.0)/(overallTotals[r]*1.0)
        
        print "*** Contamination report ***"
        print "Contam counts\t%s"%('\t'.join(map(str,contamTotals)))
        print "Contam props\t%s"%('\t'.join(map(str,contamProps)))
        print "DRB1 counts\t%s"%('\t'.join(map(str,realTotals)))
        print "DRB1 props\t%s"%('\t'.join(map(str,realProps)))
        print "Total counts\t%s"%('\t'.join(map(str,overallTotals)))
        print
        print "Overall Contamination\t%.2f%%"%((sum(contamTotals)*1.0/sum(overallTotals)*1.0)*100)
        print
        print "Contaminant\tCounts"
        for c in self._contamResults.items():
            print "%s\t%s"%(c[0],'\t'.join(map(str,c[1])))
            
        return sum(contamTotals)/sum(overallTotals)
 
                
            
                
            
                
            
        
        