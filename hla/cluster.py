import logging
import re
from hla.utils import *

'''
Created on 17 May 2012

@author: sh695
'''
class Clusters(object):
    '''
    loads sequence data and groups identical reads together in clusters
    This is more optimal than aligning each read individually as, e.g. a cluster of 500
    identical reads could be aligned with one calculation, rather than 500
    
    Use one of these class objects per file of sequences
    
    NB: This is a base class which DnaClusters and AaClusters derive from. Do not use directly
    '''
    
    _clusters = {} #{read:count}
    _firstBase = -1 #trim all before this position
    _lastBase = -1 #trim all after this position
    
    totalReads = 0
    totalClusters = 0
    clustergram = [] #histogram: bin number is number of reads, value in bin is number of clusters containing this many reads
    name = ""
    
    
    def set_read_region(self, first, last):
        '''
        sets region of interest from reads. Basically, trimming
        
        args:
            first    Start position
            last     End position
        '''
        
        logging.info('Setting read region: %d-%d',first,last)
        if first >= 0 and last > first:
            self._firstBase = first
            self._lastBase = last
        else:
            logging.warn('Region out of range. Not set')
    
    
    def get_clusters(self):
        '''
        gets clusters
        
        return:
            clusters: dict of format {sequence:count}
        '''
        
        return self._clusters
    
    
    def load_from_fq(self, fileName):
        '''
        loads and clusters all reads within a fastq file
        
        args:
            fileName    fastq file name
        '''
        raise NotImplementedError("Subclasses are responsible for creating this method")
 
                

    def stats(self):
        '''
        prints basic stats regarding clusters to the console
        '''
        
        print "%d reads across %d _clusters" % (self.totalReads, self.totalClusters)
        
        #find largest cluster
        largestClusterSize = 0
        for cluster in self._clusters.items():
            clusterSize = cluster[1]
            if clusterSize > largestClusterSize:
                largestClusterSize = clusterSize
        
        print "Largest cluster contains %d reads" % (largestClusterSize) 
        
        #create histogram
        self.clustergram = [0]*(largestClusterSize+1)
        for cluster in self._clusters.items():
            clusterSize = cluster[1]
            self.clustergram[clusterSize]+=1
        
        print "%d _clusters contain only one reads" % (self.clustergram[1])
        
        for x in range(1,largestClusterSize+1):
            if self.clustergram[x] >0:
                print "%d\t%d" %(x, self.clustergram[x])
    
    def export_clusters(self, outFileName, minClusterSize = 1):
        outFile = open(outFileName, "w")
        for cluster in self._clusters.items():
            sequence = cluster[0]
            count = cluster[1]
            if count >= minClusterSize:
                outFile.write(">%d\n%s\n"%(count,sequence))
                
        outFile.close()
    
    
    
    
    

class DnaClusters(Clusters):
    '''
    loads sequence data and groups identical reads together in clusters
    This is more optimal than aligning each read individually as, e.g. a cluster of 500
    identical reads could be aligned with one calculation, rather than 500
    
    Use one of these class objects per file of sequences
    '''
    
    
    def load_from_fq(self, fileName):
        '''
        loads and clusters all reads within a fastq file
        
        args:
            fileName    fastq file name
        '''
        
        logging.info('Generating DNA sequence clusters from fastq file: %s', fileName)
        
        #prepare variables
        self.name = fileName.split("/")[-1].split(".")[0] #e.g. /path/to/myFastq.fq > myFastq
        self._clusters = {}
        self.clustergram = []*0
        self.totalReads = 0
        self.totalClusters = 0
        trimming = self._firstBase >= 0 and self._lastBase > self._firstBase
        
        #open file
        try:
            fqFile = open(fileName)
        except IOError:
            logging.error('Could not open %s', fileName)
        
        #parse file
        dnaRegEx = re.compile('^[ACGTN]+$')
        for line in fqFile:
            isSequenceLine = dnaRegEx.match(line) != None
            if isSequenceLine is True:                
                self.totalReads+=1
                sequence = line.strip()
                if trimming:
                    sequence = line[self._firstBase : self._lastBase]
                if sequence not in self._clusters:
                    self._clusters[sequence] = 1
                else:
                    self._clusters[sequence] +=1
       
        
        #clean-up
        fqFile.close()
        self.totalClusters = len(self._clusters)
        logging.info('Generated %d clusters from %d sequences', self.totalClusters, self.totalReads)
   
                
                
                
                
class AaClusters(Clusters):
    '''
    loads sequence data, translates in 3 frames on top strand and groups identical sequences together in clusters
    This is more optimal than aligning each read individually as, e.g. a cluster of 500
    identical reads could be aligned with one calculation, rather than 500
    
    Use one of these class objects per file of sequences
    '''
    
    
    def load_from_fq(self, fileName):
        '''
        loads and clusters all reads within a fastq file
        
        args:
            fileName    fastq file name
        '''
        
        logging.info('Generating translated sequence clusters from fastq file: %s', fileName)
        
        #prepare variables
        self.name = fileName.split("/")[-1].split(".")[0] #e.g. /path/to/myFastq.fq > myFastq
        self._clusters = {}
        self.clustergram = []*0
        self.totalReads = 0
        self.totalClusters = 0
        trimming = self._firstBase >= 0 and self._lastBase > self._firstBase
        
        #open file
        try:
            fqFile = open(fileName)
        except IOError:
            logging.error('Could not open %s', fileName)
        
        #parse file
        dnaRegEx = re.compile('^[ACGTN]+$')
        for line in fqFile:
            isSequenceLine = dnaRegEx.match(line) != None
            if isSequenceLine is True:                
                self.totalReads+=1
                sequence = line.strip()
                if trimming:
                    sequence = line[self._firstBase : self._lastBase]
                #add all 3 translations (top strand) to clusters)
                translations = translate_3_frames(sequence)
                for translation in translations:
                    if translation not in self._clusters:
                        self._clusters[translation] = 1
                    else:
                        self._clusters[translation] +=1
       
        
        #clean-up
        fqFile.close()
        self.totalClusters = len(self._clusters)
        logging.info('Generated %d clusters from %d sequences', self.totalClusters, self.totalReads)
   
                


    
    
                

        
                        
                    
                
                
        

        