import math
import logging

'''
Created on 17 May 2012

@author: sh695
'''

class Aligner(object):
    '''
    for carrying out alignments and storing subsequent data.
    '''
    _counts = {}#{clusterName : {targetName:hitCount} }
    _results = {}#{targetName:[hits,hits.....hits]}
    
    
    def align(self,clusters,refs):
        '''
        derives hits from clusters to references
        
        args:
            clusters    Clusters class object
            references  Reference class object
        
        return:
            Nothing
        '''
        
        logging.info('Aligning %d clusters from %s to %d references',len(clusters.get_clusters()), clusters.name, len(refs.get_references()))
        
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
                refSeq = ref[1]
                
                #keep track of references hit by this cluster
                if querySeq in refSeq:
                    hitTargets.append(refName)
                    
            #now recalibrate score based on number of genes hit
            #If a cluster only hits one target, that has more weight than if a cluster hits 50 reference targets
            if len(hitTargets) > 0:
                #score = 1.0/math.log10(len(hitTargets)+1) TODO : Restore this!
                score = 1.0/(len(hitTargets)*1.0)
                score *= queryCount #also weight by number of reads in cluster
                for target in hitTargets:
                    if target not in geneHitsForCluster:
                        geneHitsForCluster[target] =score
                    else:
                        geneHitsForCluster[target]+=score
        
        #add to _counts
        self._counts[clusters.name] = geneHitsForCluster
    
    
    
    def calculate_alignment_cross_talk(self,clusters,refs):
        '''
        generates a grid of geneName x geneName and shows number of reads that map to each
        e.g.
        
            g1    g2    g3    g4
        g1  x     0     10    0
        g2  0     x     25    5
        g3  10    25    x     0
        g4  0     5     0     x
        
        args:
            clusters    Clusters class object
            references  Reference class object
        
        return:
            grid
        '''
        logging.info('Calculating cross-talk')
        
        #prepare variables
        geneHitsForCluster = {}
        
        queryClusters = clusters.get_clusters()
        allHitTargets = set()
        
        #find out what genes are mapped to
        for cluster in queryClusters.items():
            querySeq = cluster[0]
            queryCount = cluster[1]
            for ref in refs.get_references().items():
                refName = ref[0]
                refSeq = ref[1]
                
                #keep track of references hit by this cluster
                if querySeq in refSeq:
                    allHitTargets.add(refName)
                    
        
        
        #assign a keyId to each gene name
        idToNameDict = {}
        nameToIdDict = {}
        
        currId = 0
        print allHitTargets
        
        for ref in allHitTargets:
            print ref,currId
            idToNameDict[currId] = ref
            nameToIdDict[ref] = currId
            currId+=1

        
        #generate a grid
        grid = []
        for y in range(0,len(allHitTargets)):
            grid.append([])
            for x in range(0,len(allHitTargets)):
                grid[y].append(0.0)
    #    grid =[[0.0]*len(allHitTargets)]*len(allHitTargets)
        print grid
        print "Number of targets hit by these clusters: %d"%len(allHitTargets)
        
        totScore =0
        #For each cluster of reads, see if they match any reference sequences
        for cluster in queryClusters.items():
            querySeq = cluster[0]
            queryCount = cluster[1]
            hitTargets = []
            for ref in refs.get_references().items():
                refName = ref[0]
                refSeq = ref[1]
                
                #keep track of references hit by this cluster
                if querySeq in refSeq:
                    hitTargets.append(nameToIdDict[refName])
        
            #add data to grid
            
            for first in range (0,len(hitTargets)):
                for second in range(0,len(hitTargets)):
                    if first!=second:
                        grid[hitTargets[first]][hitTargets[second]]+=queryCount
                        grid[hitTargets[second]][hitTargets[first]]+=queryCount
                        print "[%d][%d]+=%.2f"%(hitTargets[first], hitTargets[second],queryCount)
            if len(hitTargets)==1:
                grid[hitTargets[0]][hitTargets[0]]+=1000*queryCount
                print "SAME:%d"%hitTargets[first]
        #    for hit1 in hitTargets:
        #        for hit2 in hitTargets:
        #            if hit1!=hit2:
        #                grid[hit1][hit2]+=queryCount
            
            totScore+=queryCount
            if len(hitTargets) >1:
                raw_input("round complete")
            print totScore
            
        print grid
        
        outFile = open("grid2.txt","w")
        #print header
        outFile.write("\t")
        for head in range(0,len(allHitTargets)):
            outFile.write("%s\t"%idToNameDict[head])
        outFile.write("\n")
        
        for y in range(0,len(allHitTargets)):
            outFile.write("%s\t"%idToNameDict[y])
            for x in range(0,len(allHitTargets)):
                outFile.write("%d\t"%grid[y][x])
            outFile.write("\n")
        outFile.close()
        raw_input("Written grid to file")
                    
        
        

        for cluster in queryClusters.items():
            querySeq = cluster[0]
            queryCount = cluster[1]
            hitTargets = []
            for ref in refs.get_references().items():
                refName = ref[0]
                refSeq = ref[1]
                
                #keep track of references hit by this cluster
                if querySeq in refSeq:
                    hitTargets.append(refName)
                    
            #now recalibrate score based on number of genes hit
            #If a cluster only hits one target, that has more weight than if a cluster hits 50 reference targets
            if len(hitTargets) > 0:
                score = 1.0/math.log((len(hitTargets)+1) , 10)
                score *= queryCount #also weight by number of reads in cluster
                for target in hitTargets:
                    if target not in geneHitsForCluster:
                        geneHitsForCluster[target] =score
                    else:
                        geneHitsForCluster[target]+=score
        
        #add to _counts
        self._counts[clusters.name] = geneHitsForCluster
        
        
        
        
    def align_unique_only(self,clusters,refs):
        '''
        derives hits from clusters to references, BUT only accepting hits which only hit one target
        It's total crap. So few reads are actually uniquely mapping to one target there's no point
        It biases heavily against all but a few HLAs
        
        args:
            clusters    Clusters class object
            references  Reference class object
        
        return:
            Nothing
        '''
        
        logging.info('Aligning %d clusters from %s to %d references',len(clusters.get_clusters()), clusters.name, len(refs.get_references()))
        
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
                refSeq = ref[1]
                
                #keep track of references hit by this cluster
                if querySeq in refSeq:
                    hitTargets.append(refName)
                    
            #only accept it unique
            if len(hitTargets) == 1:
                score = queryCount 
                for target in hitTargets:
                    if target not in geneHitsForCluster:
                        geneHitsForCluster[target] =score
                    else:
                        geneHitsForCluster[target]+=score
        
        #add to _counts
        self._counts[clusters.name] = geneHitsForCluster
   
        
    
    def print_counts(self):
        '''
        print results out. Use only after compile results
        '''
        for dataSet in self._counts.items():
            readName = dataSet[0]
            for entry in self._counts[readName].items():
                print readName,entry[0],entry[1]
                
                
                
                
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
                
                #trim down geneName
                colonSplit = geneName.split(":")
                shorterName = "%s:%s"%(colonSplit[0],colonSplit[1]) 
                if shorterName not in abbreviatedHits:
                    abbreviatedHits[shorterName]=geneScore
                else:
                    abbreviatedHits[shorterName]+=geneScore  
            #replace old scores with new ones
            self._counts[readName] = abbreviatedHits
        
        #get list of name of all genes that were hit
        hitGenes = set()
        for dataSet in self._counts.items():
            for entry in dataSet[1].items():
                hitGenes.add(entry[0])
                
        #print "Found %d genes"%len(hitGenes)
        
        self._results ={}
        numReads = len(self._counts)
        for gene in hitGenes:
            self._results[gene]=[0.0] * numReads
            counter = 0
            for dataSet in self._counts.items():
                if gene in dataSet[1]:
                    self._results[gene][counter] = dataSet[1][gene]
                #else:
                #    self._results[gene][counter]
                counter+=1
            #print "%s\t%s\t%.2f" %(gene, self._results[gene],sum(self._results[gene]))
            

    
    def get_results(self):
        '''
        gets results. Be sure to have run compile_results() first
        
        return:
            results: dict of format e.g., {DRB1*01:01 , [1,100,20,4]}
        '''
        return self._results
    
    
    
    
    

        
        
        
        
                
            
            