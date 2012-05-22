import logging
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

'''
Created on 17 May 2012

@author: sh695
'''

class Reference(object):
    '''
    classdocs
    '''


    _refs = {}#{name:sequence}
    

    
    
    
    def load_reference_tabbed(self, refName):
        '''
        load references from a tab delimited file. First column is gene name, second column is sequence
        '''
        
        logging.info('Loading references from tabbed file: %s',refName)
        
        
        #open file
        try:
            refFile = open(refName)
        except IOError:
            logging.error("Unable to open %s",refName)
        
        #parse file
        for entry in refFile:
            entry = entry.strip()
            entrySplit = entry.split("\t")
            if len(entrySplit)==2:
                name = entrySplit[0]
                sequence = entrySplit[1]
                self._refs[name] = sequence
        
        #clean up
        logging.info('Loaded %d reference sequences',len(self._refs))
        refFile.close()
    
                
            
            
            
    def num_references(self):
        '''
        returns number of loaded references
        '''
        return len(self._refs)
    
    
    
    def get_references(self):
        return self._refs
    
    
    def get_translated_references(self):
        return self._translated_refs


        
                