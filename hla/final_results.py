'''
Created on 20 Jun 2012

@author: sh695
'''

from hla.allele_status import Status


class FinalResultBundle(object):
    '''
    Contains all final results information, i.e. several FinalResult objects plus wellId etc
    '''
    
    results = {}    #map of {resultScore : correspondingFinalResult}
    resultsSortedList = [] #sorted list by score of FinalResult instances
    wellId = ""
    comments = ""   
    calledAlleles = ""
    accepted_alleles = []
    
    
    
    def result_string(self):
        
        #sort alleles in to order
        self.accepted_alleles = []
        count = 0


        for entry in sorted(self.results.iteritems(), key=lambda (x,y): float(x), reverse=True):
            score = entry[0]
            result = entry[1]
            self.resultsSortedList.append(result)

            
        #prepare comments and top alleles
        self._calculate_called_alleles()
        self._generate_comment()
        
        #build string    
        finalString = "%s"%self.wellId
        finalString = "%s\t%s\t%s"%(finalString, '\t'.join(map(str,self.accepted_alleles)), self.comments)
        for entry in self.resultsSortedList:
            print "Adding: ",entry.get_result_string()
            finalString = "%s\t%s"%(finalString, entry.get_result_string())
            
        
        return finalString
        
        
        
    def _calculate_called_alleles(self):
        
        count = 0
        
        for i in range(0,2):
            result = self.resultsSortedList[i]
            if result.status == Status.Pass():
                #carve down allele
                modifiedName = "%s%s"%(result.targetName[5:7],result.targetName[8:])
                self.accepted_alleles.append(modifiedName)

        if len(self.accepted_alleles) == 0:
            self.accepted_alleles.append("*")
            self.accepted_alleles.append("*")
        if len(self.accepted_alleles) == 1:
            self.accepted_alleles.append(self.accepted_alleles[0])
            
            
    
    def _generate_comment(self):
        #Have both alleles passed?
        if "*" not in self.accepted_alleles:
            #Homozygous?
            if self.accepted_alleles[0] == self.accepted_alleles[1]:
                self.comments = "%sHomozygous; "%(self.comments)
            #third allele close to top two?
            if len(self.resultsSortedList)>2 and self.resultsSortedList[2].status == Status.Warn_not_top_two():
                    self.comments = "%sSecond score close to third; "%(self.comments)
            #second allele failed only because of score differential (e.g. score 1 = 4.5 and score 2 = 3.0?
            if self.resultsSortedList[1].status == Status.Fail_score_diff():
                self.comments = "%sSecond allele score could be acceptable, but a lot lower than first; "%(self.comments)
                
        else:
            #two fails
            if self.accepted_alleles[0]=="*" and self.accepted_alleles[1]=="*":
                self.comments = "%sBoth calls failed; "%(self.comments)
            #only one failed
            else:     
                self.comments = "%sOne call failed; "%(self.comments)
            
       
                    
                        
    
            
class FinalResult(object):
    '''
    Contains all information relevant to result
    '''
    
    targetName = ""
    status = ""
    rawCounts = []
    modifiedCounts = []
    score = 0.0
    
        
        
    def get_result_string(self):
        resultString = "%s\t%s\t%.2f\t%s"%(self.targetName , self.rawCounts , self.score, self.status)
        return resultString
        
        