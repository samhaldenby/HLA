'''
Created on 29 May 2012

@author: sh695
'''

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
        contamTotals = [0]*4
        realTotals = [0]*4
        overallTotals = [0]*4
        
        for entry in self._rawAlignmentResults.items():
            refName = entry[0]
            refScores = entry[1]
            
            #split into contam/real depending on name
            if "DRB1" in refName:
                self._realResults[refName]=refScores
                for r in range(0,4):
                    realTotals[r]+=refScores[r]
            else:
                self._contamResults[refName]=refScores
                for r in range(0,4):
                    contamTotals[r]+=refScores[r]
                    
            for r in range(0,4):
                overallTotals[r]+=refScores[r]
                
        contamProps = [0.0]*4
        realProps = [0.0]*4
        
        for r in range(0,4):
            contamProps[r] = contamTotals[r]/overallTotals[r]
            realProps[r] = realTotals[r]/overallTotals[r]
        
        print "CONTAM: %s"%contamTotals
        print "      : %s"%contamProps
        print "REAL  : %s"%realTotals
        print "      : %s"%realProps
        print "TOTAL : %s"%overallTotals
      
        
        for c in self._contamResults.items():
            print c
            
        return sum(contamTotals)/sum(overallTotals)
 
                
            
                
            
                
            
        
        