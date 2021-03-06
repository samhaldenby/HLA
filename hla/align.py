import math
import logging
from hla.utils import *
import commands

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
    
    Order = [""] * 4
    
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
            clusterName = dataSet[0]
            print "<<< CLUSTERNAME: ", clusterName
            hits = dataSet[1]
            for entry in hits.items():
                targetName = entry[0]
                hitGenes.add(targetName)
                
        #print "Found %d genes"%len(hitGenes)
        
        self._results ={}
        numReads = len(self._counts)
        for gene in hitGenes:
            self._results[gene]=[0] * numReads
            counter = 0
            for dataSet in self._counts.items():
                clusterName = dataSet[0]
                hits = dataSet[1]
                Aligner.Order[counter]=clusterName
                if gene in hits:
                    self._results[gene][counter] = hits[gene]
                #else:
                #    self._results[gene][counter]
                counter+=1
        print Aligner.Order
            #print "%s\t%s\t%.2f" %(gene, self._results[gene],sum(self._results[gene]))
            

    def report_basic_stats(self):
        '''
        basic alignment stats can be displayed with this
        '''
        
        pass
    
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
                
            percentHitsByTotal = 0.0
            percentMappedByTotal = 0.0
            if totalReads > 0:
                percentHitsByTotal = (totalHits*100.0)/(totalReads*1.0)
                percentMappedByTotal = (self._mappedReadCounts[clusterName]*100.0)/(totalReads*1.0)
                
            aveHitsPerMapped = 0.0
            if  self._mappedReadCounts[clusterName] > 0:
                aveHitsPerMapped = totalHits/self._mappedReadCounts[clusterName]
            print ">> %s\t%d\t%d\t%s\t%.2f%%\t%.2f%%\t%.2f"%(clusterName,totalHits,self._mappedReadCounts[clusterName],totalReads,percentHitsByTotal,percentMappedByTotal, aveHitsPerMapped)
            
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
        
        #print "Hit targets are: %s",hitTargets
        return hitTargets
    
    
    
    
    def align(self,clusters,refs,omitHitsTo = set(), unmappedFileName = None):
        '''
        derives hits from clusters to references
        
        args:
            clusters    Clusters class object
            references  Reference class object
            omitHitsTo  A set of geneNames (abbreviated!). Any reads that map to these genes (regardless of where else they map) will not be considered
        '''

        
        logging.info('Aligning %d clusters from %s to %d references, omitting reads hitting %s',len(clusters.get_clusters()), clusters.name, len(refs.get_references()), omitHitsTo)
        print ">>>> ORDER: ",clusters.name
        
###        print "*************"
###        print "**** %s"%clusters.name
###        print "*************"
        #prepare variables
        referenceHitCounts = {} #{referenceSequenceName : hitCount}
        queryClusters = clusters.get_clusters()
        totalCount = 0
        mappedReads = 0
        
        #if writing unmapped to a file, open that file
        unmappedFile = None
        if unmappedFileName != None:
            unmappedFile = open(unmappedFileName,"w")
            
        #variables for generating stats histograms (more of a debugging technique at the mo, but could be good in final reports)
        unmappedClusterHistogram = {}
        mappedClusterHistogram = {}
        #For each cluster of reads, see if they match any reference sequences
        for cluster in queryClusters.items():
            
            clusterSequence = cluster[0]
            numReadsInCluster = cluster[1]
            totalCount+= numReadsInCluster
            
            
            
            #get a set of all targets that this cluster hits
            hitTargets = self._calc_hits_against_ref(query=clusterSequence, references=refs)
            
            
            if len(hitTargets) > 0:
                mappedReads += numReadsInCluster
###                print "MappedReads = %d (this cluster size: %d)\t(%s)"%(mappedReads, numReadsInCluster, hitTargets)
                #check that it doesn't map to any of the omitHitTo set. If it does, chuck the cluster away
                keepCluster = True
                for target in hitTargets:
                    if target in omitHitsTo or "DRB3" in target or "DRB5" in target:
                        keepCluster = False
                        break
                
                #if all is fine, add hits to each hit target   
                if keepCluster:
                    #add to mappedClusterHistogram
                    if numReadsInCluster not in mappedClusterHistogram:
                        mappedClusterHistogram[numReadsInCluster] = 1
                    else:
                        mappedClusterHistogram[numReadsInCluster] +=1
                    score = numReadsInCluster #1 hit per read in cluster
                    for target in hitTargets:
                        if target not in referenceHitCounts:
                            referenceHitCounts[target] =score
                        else:
                            referenceHitCounts[target]+=score
            #if doesn't map, write to unmapped reads file
            elif unmappedFile != None:
                if numReadsInCluster not in unmappedClusterHistogram:
                    unmappedClusterHistogram[numReadsInCluster] = 1
                else:
                    unmappedClusterHistogram[numReadsInCluster] +=1
                
###                if numReadsInCluster > 10:
###                    print "Cluster of size %d not mapped!:\t%s"%(numReadsInCluster,clusterSequence)
                unmappedFile.write("%d\t%s\n"%(numReadsInCluster, clusterSequence))

        #add to _counts
        self._counts[clusters.name] = referenceHitCounts
        self._totalReads[clusters.name] = totalCount
        self._mappedReadCounts[clusters.name] = mappedReads
###        print ">>> MappedReads[%s] = %d (this cluster size: %d)\t(%s)"%(clusters.name,mappedReads, numReadsInCluster, hitTargets)
        
        if unmappedFile != None:
            unmappedFile.close()
            
###       print
###       print "*** Unmapped"
###        print "ClusterSize\tNumClusters"
###        for entry in unmappedClusterHistogram.items():
###            binId = entry[0]
###            frequency = entry[1]
###            print "%d\t%d"%(binId,frequency)
###        print
###        print "*** Mapped"
###        print "ClusterSize\tNumClusters"
###        for entry in mappedClusterHistogram.items():
###            binId = entry[0]
###            frequency = entry[1]
###            print "%d\t%d"%(binId,frequency)
        
        
    def align_by_blast(self,clusters,refs,omitHitsTo = set(), unmappedFileName = None):
        '''
        derives hits from clusters to references but uses blast - experimental so don't use for production!
        
        args:
            clusters    Clusters class object
            references  Reference class object
            omitHitsTo  A set of geneNames (abbreviated!). Any reads that map to these genes (regardless of where else they map) will not be considered
        '''

        
        logging.info('Aligning %d clusters from %s to %d references with BLAST, omitting reads hitting %s',len(clusters.get_clusters()), clusters.name, len(refs.get_references()), omitHitsTo)
        
        
###        print "*************"
###        print "**** %s"%clusters.name
###        print "*************"
        #prepare variables
        referenceHitCounts = {} #{referenceSequenceName : hitCount}
        queryClusters = clusters.get_clusters()
        totalCount = 0
        mappedReads = 0
        
        #if writing unmapped to a file, open that file
        unmappedFile = None
        if unmappedFileName != None:
            unmappedFile = open(unmappedFileName,"w")
            
        #variables for generating stats histograms (more of a debugging technique at the mo, but could be good in final reports)
        unmappedClusterHistogram = {}
        mappedClusterHistogram = {}
        #For each cluster of reads, see if they match any reference sequences
        for cluster in queryClusters.items():
            
            clusterSequence = cluster[0]
            numReadsInCluster = cluster[1]
            totalCount+= numReadsInCluster
            #create fasta entry
            
            if numReadsInCluster >1:
                cmd = "echo %s | blastn -query - -db /scratch/sh695/HLA/subject.fa -ungapped -outfmt 6 "%clusterSequence.strip()
                #print cmd
             #   blastResults = subprocess.check_output(cmd,shell=True)
                status, blastResults = commands.getstatusoutput(cmd)
                #get top hits
                blastScoreList = []
                blastHitList = []
                #print blastResults
                hitTargets = set()
                blastTopScore = 0
                firstHit=True
                #grab hits with top score
                for blastLine in blastResults.split("\n"):
                    splitByTab = blastLine.split("\t")
                    if len(splitByTab) ==12:
                        blastHitName = splitByTab[1]
                        thisBlastScore = float(splitByTab[11])
                        if firstHit:
                            blastTopScore = thisBlastScore
                            firstHit = False
                        if thisBlastScore == blastTopScore:#
                            hitTargets.add(blastHitName)
                        else:
                            break
                    
                #print hitTargets   
                    
                        
               #     print splitByTab[1],splitByTab[11]
                #sort for top hit
                keydict = dict(zip(blastHitList, blastScoreList))
                blastHitList.sort(key=keydict.get)
               # print blastHitList
    
    
                    
                #resultsSplitByTab = blastResults.split("\t")
                
                #print resultsSplitByLine
                
                
                
                
                #get a set of all targets that this cluster hits
                #hitTargets = self._calc_hits_against_ref_with_blast(query=clusterSequence, references=refs)
                
                
                if len(hitTargets) > 0:
                    mappedReads += numReadsInCluster
    ###                print "MappedReads = %d (this cluster size: %d)\t(%s)"%(mappedReads, numReadsInCluster, hitTargets)
                    #check that it doesn't map to any of the omitHitTo set. If it does, chuck the cluster away
                    keepCluster = True
                    for target in hitTargets:
                        if target in omitHitsTo:
                            keepCluster = False
                            break
                    
                    #if all is fine, add hits to each hit target   
                    if keepCluster:
                        #add to mappedClusterHistogram
                        if numReadsInCluster not in mappedClusterHistogram:
                            mappedClusterHistogram[numReadsInCluster] = 1
                        else:
                            mappedClusterHistogram[numReadsInCluster] +=1
                        score = numReadsInCluster #1 hit per read in cluster
                        for target in hitTargets:
                            if target not in referenceHitCounts:
                                referenceHitCounts[target] =score
                            else:
                                referenceHitCounts[target]+=score
                #if doesn't map, write to unmapped reads file
                elif unmappedFile != None:
                    if numReadsInCluster not in unmappedClusterHistogram:
                        unmappedClusterHistogram[numReadsInCluster] = 1
                    else:
                        unmappedClusterHistogram[numReadsInCluster] +=1
                    
    ###                if numReadsInCluster > 10:
    ###                    print "Cluster of size %d not mapped!:\t%s"%(numReadsInCluster,clusterSequence)
                    unmappedFile.write("%d\t%s\n"%(numReadsInCluster, clusterSequence))

        #add to _counts
        self._counts[clusters.name] = referenceHitCounts
        self._totalReads[clusters.name] = totalCount
        self._mappedReadCounts[clusters.name] = mappedReads
###        print ">>> MappedReads[%s] = %d (this cluster size: %d)\t(%s)"%(clusters.name,mappedReads, numReadsInCluster, hitTargets)
        
        if unmappedFile != None:
            unmappedFile.close()
            

    def align_by_blat(self,clusters,refs,omitHitsTo = set(), unmappedFileName = None):
        '''
        derives hits from clusters to references but uses blat - experimental so don't use for production!
        
        args:
            clusters    Clusters class object
            references  Reference class object
            omitHitsTo  A set of geneNames (abbreviated!). Any reads that map to these genes (regardless of where else they map) will not be considered
        '''

        
        logging.info('Aligning %d clusters from %s to %d references with BLAST, omitting reads hitting %s',len(clusters.get_clusters()), clusters.name, len(refs.get_references()), omitHitsTo)
        
        
###        print "*************"
###        print "**** %s"%clusters.name
###        print "*************"
        #prepare variables
        referenceHitCounts = {} #{referenceSequenceName : hitCount}
        queryClusters = clusters.get_clusters()
        totalCount = 0
        mappedReads = 0
        
        #if writing unmapped to a file, open that file
        unmappedFile = None
        if unmappedFileName != None:
            unmappedFile = open(unmappedFileName,"w")
            
        #variables for generating stats histograms (more of a debugging technique at the mo, but could be good in final reports)
        unmappedClusterHistogram = {}
        mappedClusterHistogram = {}
        #For each cluster of reads, see if they match any reference sequences
        for cluster in queryClusters.items():
            
            clusterSequence = cluster[0]
            numReadsInCluster = cluster[1]
            totalCount+= numReadsInCluster
            #create fasta entry
            
            if numReadsInCluster >1:
                #cmd = blat DRB_exon2_3.8.0.fa  A05_clusters.fa A05_noHead.psl -noHead -fastMap -minScore=130 -minIdentity=100
                cmd = "echo %s | blastn -query - -db /scratch/sh695/HLA/subject.fa -ungapped -outfmt 6 "%clusterSequence.strip()
                #print cmd
             #   blastResults = subprocess.check_output(cmd,shell=True)
                status, blastResults = commands.getstatusoutput(cmd)
                #get top hits
                blastScoreList = []
                blastHitList = []
                #print blastResults
                hitTargets = set()
                blastTopScore = 0
                firstHit=True
                #grab hits with top score
                for blastLine in blastResults.split("\n"):
                    splitByTab = blastLine.split("\t")
                    if len(splitByTab) ==12:
                        blastHitName = splitByTab[1]
                        thisBlastScore = float(splitByTab[11])
                        if firstHit:
                            blastTopScore = thisBlastScore
                            firstHit = False
                        if thisBlastScore == blastTopScore:#
                            hitTargets.add(blastHitName)
                        else:
                            break
                    
                #print hitTargets   
                    
                        
               #     print splitByTab[1],splitByTab[11]
                #sort for top hit
                keydict = dict(zip(blastHitList, blastScoreList))
                blastHitList.sort(key=keydict.get)
               # print blastHitList
    
    
                    
                #resultsSplitByTab = blastResults.split("\t")
                
                #print resultsSplitByLine
                
                
                
                
                #get a set of all targets that this cluster hits
                #hitTargets = self._calc_hits_against_ref_with_blast(query=clusterSequence, references=refs)
                
                
                if len(hitTargets) > 0:
                    mappedReads += numReadsInCluster
    ###                print "MappedReads = %d (this cluster size: %d)\t(%s)"%(mappedReads, numReadsInCluster, hitTargets)
                    #check that it doesn't map to any of the omitHitTo set. If it does, chuck the cluster away
                    keepCluster = True
                    for target in hitTargets:
                        if target in omitHitsTo:
                            keepCluster = False
                            break
                    
                    #if all is fine, add hits to each hit target   
                    if keepCluster:
                        #add to mappedClusterHistogram
                        if numReadsInCluster not in mappedClusterHistogram:
                            mappedClusterHistogram[numReadsInCluster] = 1
                        else:
                            mappedClusterHistogram[numReadsInCluster] +=1
                        score = numReadsInCluster #1 hit per read in cluster
                        for target in hitTargets:
                            if target not in referenceHitCounts:
                                referenceHitCounts[target] =score
                            else:
                                referenceHitCounts[target]+=score
                #if doesn't map, write to unmapped reads file
                elif unmappedFile != None:
                    if numReadsInCluster not in unmappedClusterHistogram:
                        unmappedClusterHistogram[numReadsInCluster] = 1
                    else:
                        unmappedClusterHistogram[numReadsInCluster] +=1
                    
    ###                if numReadsInCluster > 10:
    ###                    print "Cluster of size %d not mapped!:\t%s"%(numReadsInCluster,clusterSequence)
                    unmappedFile.write("%d\t%s\n"%(numReadsInCluster, clusterSequence))

        #add to _counts
        self._counts[clusters.name] = referenceHitCounts
        self._totalReads[clusters.name] = totalCount
        self._mappedReadCounts[clusters.name] = mappedReads
###        print ">>> MappedReads[%s] = %d (this cluster size: %d)\t(%s)"%(clusters.name,mappedReads, numReadsInCluster, hitTargets)
        
        if unmappedFile != None:
            unmappedFile.close()
            
            
            
        
    def modify_scores_for_read_imbalances(self):
        #calculate read imbalances
        modifiers = []
        total = 0
        for entry in self._mappedReadCounts.items():
            total+=entry[1]
        
        for entry in self._mappedReadCounts.items():
            ampliconName = entry[0]
            count = entry[1]
            if count > 0:
                modifiers.append((total*1.0)/(count*1.0))
            else:
                modifiers.append(0.0)   #TODO: Is this correct?!
            
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

    

        
        
        
        
                
            
            