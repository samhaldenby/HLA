import math
import logging
from hla.utils import *

'''
Created on 17 May 2012

@author: sh695
'''

class Aligner(object):
    '''
    for carrying out alignments and storing subsequent data.
    
    NB: This is an abstract base class. Use either DnaAligner or AaAligner
    '''
    _counts = {}#{clusterName : {targetName:hitCount} }
    _results = {}#{targetName:[hits,hits.....hits]} modified by number of ambiguous hits
    #_rawResults = {}#unmodified hit counts, in same format as results
    _totalReads = {}#{clusterName : numOfReadsInTotal}
    _mappedReadCounts = {}#{clusterName : numOfMappedReads}
    
    def align(self,clusters,refs, omitHitsTo = set()):
        '''
        derives hits from clusters to references
        
        args:
            clusters    Clusters class object
            references  Reference class object
            omitHitsTo  A set of geneNames (abbreviated!). Any reads that map to these genes (regardless of where else they map) will not be considered
        
        return:
            Nothing
        '''
        raise NotImplementedError("Subclasses are responsible for creating this method")
    
    
    
    def compile_results(self):
        '''
        once all alignments have been carried out, do this.
        It will create a dict of hits for each gene e.g. a row
        might look like:
        {DRB1*01:01 , [1,100,20,4]}
        '''
        
        logging.info('Compiling results')
        
        #abbreviate gene names and consolidate e.g. DRB1*01:01:02:04 -> DRB1*01:01
        for read in self._counts.items():
            readName = read[0]
            hits = read[1]
            abbreviatedHits = {}
            for hit in hits.items():
                geneName = hit[0]
                geneScore = hit[1]
                
#                #trim down geneName - this won't matter in AaAlignment as will have already been trimmed during loading of reference
#                shorterName = shorten_name(geneName)
                
#                if shorterName not in abbreviatedHits:
#                    abbreviatedHits[shorterName]=geneScore
#                else:
#                    abbreviatedHits[shorterName]+=geneScore  
            #replace old scores with new ones
#            self._counts[readName] = abbreviatedHits
        
        #get list of name of all genes that were hit
        hitGenes = set()
        for dataSet in self._counts.items():
            for entry in dataSet[1].items():
                hitGenes.add(entry[0])
                
        #print "Found %d genes"%len(hitGenes)
        
        self._results ={}
        numReads = len(self._counts)
        for gene in hitGenes:
            self._results[gene]=[0] * numReads
            counter = 0
            for dataSet in self._counts.items():
                if gene in dataSet[1]:
                    self._results[gene][counter] = dataSet[1][gene]
                #else:
                #    self._results[gene][counter]
                counter+=1
            #print "%s\t%s\t%.2f" %(gene, self._results[gene],sum(self._results[gene]))
            

    def report_basic_stats(self):
        '''
        basic alignment stats can be displayed with this
        '''
        
        print 
        print "*** Basic stats ***"
        for entry in self._totalReads.items():
            clusterName = entry[0]
            totalReads = entry[1]
            #now count how many reads mapped
            totalHits =0 
            for hitEntry in self._counts[clusterName].items():
                hitCount = hitEntry[1]
                totalHits +=hitCount
        
            print "%s\t%s\t%s\t%.2f%%"%(clusterName,totalHits,totalReads,(totalHits*100.0)/(totalReads*1.0))
            
        print
        
    def get_results(self):
        '''
        gets results. Be sure to have run compile_results() first
        
        return:
            results: dict of format e.g., {DRB1*01:01 , [1,100,20,4]}
        '''
        return self._results
    
    
    
    
    
    
    
class DnaAligner(Aligner):
    '''
    for carrying out nucleotide alignments and storing subsequent data.
    '''
    
    def report_basic_stats(self):
        '''
        basic alignment stats can be displayed with this
        '''
        
        print ">>"
        print ">> *** Basic stats ***"
        print ">> Cluster\tTotalHits\tTotalMappedReads\tTotalReads\t%Hits/Total\t$Mapped/Total\tAve Hits Per Mapped"
        for entry in self._totalReads.items():
            clusterName = entry[0]
            totalReads = entry[1]
            #now count how many reads mapped
            totalHits =0 
            for hitEntry in self._counts[clusterName].items():
                hitCount = hitEntry[1]
                totalHits +=hitCount
        
            print ">> %s\t%d\t%d\t%s\t%.2f%%\t%.2f%%\t%.2f"%(clusterName,totalHits,self._mappedReadCounts[clusterName],totalReads,(totalHits*100.0)/(totalReads*1.0),(self._mappedReadCounts[clusterName]*100.0)/(totalReads*1.0), totalHits/self._mappedReadCounts[clusterName])
            
        print ">>"
        
        
    def _calc_hits_against_ref(self, query, references):
        '''
        takes a query sequence and Reference class and calculates the targets that the query sequence hits
        
        args:
            query    DNA sequence string
            references    filled instance of Reference class
            
        returns:
            set of reference sequence names that the query maps to
        '''
        
        hitTargets = set()
        for ref in references.get_references().items():
            refName = ref[0]
            refSeqs = ref[1]
                
            for refSeq in refSeqs:
                if query in refSeq:
                    hitTargets.add(refName)
                    
        return hitTargets
    
    
    
    
    def align(self,clusters,refs,omitHitsTo = set()):
        '''
        derives hits from clusters to references
        
        args:
            clusters    Clusters class object
            references  Reference class object
            omitHitsTo  A set of geneNames (abbreviated!). Any reads that map to these genes (regardless of where else they map) will not be considered
        '''

        
        logging.info('Aligning %d clusters from %s to %d references, omitting reads hitting %s',len(clusters.get_clusters()), clusters.name, len(refs.get_references()), omitHitsTo)
        
        
        
        #prepare variables
        referenceHitCounts = {} #{referenceSequenceName : hitCount}
        queryClusters = clusters.get_clusters()
        totalCount = 0
        mappedReads = 0
        #For each cluster of reads, see if they match any reference sequences
        for cluster in queryClusters.items():
            querySeq = cluster[0]
            queryCount = cluster[1]
            totalCount+= queryCount
            
            
            
            #get a set of all targets that this cluster hits
            hitTargets = self._calc_hits_against_ref(query=querySeq, references=refs)
            
            
            if len(hitTargets) > 0:
                mappedReads += queryCount
                #check that it doesn't map to any of the omitHitTo set. If it does, chuck the cluster away
                keepCluster = True
                for target in hitTargets:
                    if target in omitHitsTo:
                        keepCluster = False
                        break
                
                #if all is fine, add hits to each hit target   
                if keepCluster:
                    score = queryCount #1 hit per read in cluster
                    for target in hitTargets:
                        if target not in referenceHitCounts:
                            referenceHitCounts[target] =score
                        else:
                            referenceHitCounts[target]+=score
        
        #add to _counts
        self._counts[clusters.name] = referenceHitCounts
        self._totalReads[clusters.name] = totalCount
        self._mappedReadCounts[clusters.name] = mappedReads
        
        
    def modify_scores_for_read_imbalances(self):
        #calculate read imbalances
        modifiers = []
        total = 0
        for entry in self._mappedReadCounts.items():
            total+=entry[1]
        
        for entry in self._mappedReadCounts.items():
            ampliconName = entry[0]
            count = entry[1]
            modifiers.append((total*1.0)/(count*1.0))
            
        #Now modify all results
        modifiedResults = {}
        for entry in self._results.items():
            targetName = entry[0]
            allHits = entry[1] #e.g. [100,500,200,140]
            counter = 0
            modifiedHits = [0.0]*len(modifiers)
            for i in range(0,len(modifiers)):
                allHits[i] = (allHits[i]*1.0) * modifiers[i] * (1.0/len(modifiers)*1.0)
            #modifiedResults[targetName] = modifiedHits
            
        return modifiedResults
        
class AaAligner(Aligner):
    '''
    for carrying out alignments of translated data and storing subsequent data.
    '''
    
    
    def align(self,clusters,refs,  omitHitsTo = set()):
        '''
        derives hits from clusters to references
        
        args:
            clusters    Clusters class object
            references  Reference class object
        
        return:
            Nothing
        '''
        
        logging.info('Aligning %d translated clusters from %s to %d translated references',len(clusters.get_clusters()), clusters.name, len(refs.get_references()))
        
        #prepare variables
        geneHitsForCluster = {}
        
        queryClusters = clusters.get_clusters()
        #For each cluster of reads, see if they match any reference sequences
        

        for cluster in queryClusters.items():
            querySeq = cluster[0]
            queryCount = cluster[1]
            hitTargets = []
            for ref in refs.get_references().items():
                refName = ref[0]
                refTranslations = ref[1]
                
                #keep track of references hit by this cluster
                for refTranslation in refTranslations:
                    
                    if querySeq in refTranslation:
                        hitTargets.append(refName)
                    
            #now recalibrate score based on number of genes hit
            #If a cluster only hits one target, that has more weight than if a cluster hits 50 reference targets
            if len(hitTargets) > 0:
                #check that it doesn't map to any of the omitHitTo set. If it does, chuck the cluster away
                keepCluster = True
                for target in hitTargets:
                    if target in omitHitsTo:
                        keepCluster = False
                        break
                
                #if all ok, add hits to each hit target   
                if keepCluster:
                    score = queryCount
                   # score = 1.0/(len(hitTargets)*1.0)
                    #score *= queryCount #also weight by number of reads in cluster
                    for target in hitTargets:
                        if target not in geneHitsForCluster:
                            geneHitsForCluster[target] =score
                        else:
                            geneHitsForCluster[target]+=score
        
        #add to _counts
        self._counts[clusters.name] = geneHitsForCluster

    

        
        
        
        
                
            
            