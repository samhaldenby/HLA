import logging
from hla.utils import *
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
    

    _refs = {}#{name:sequence} or {name:[translation1, translation2, translation3]}
    
    
    def load_reference_tabbed(self, refName):
        '''
        load references from a tab delimited file. First column is gene name, second column is sequence
        '''
        raise NotImplementedError("Subclasses are responsible for creating this method")
        
        
        
    def num_references(self):
        '''
        returns number of loaded references
        '''
        return len(self._refs)
    
    
    
    def get_references(self):
        return self._refs
        
        
        
        
        
        
class DnaReference(Reference):
    '''
    classdocs
    '''
    
    
    def load_reference_tabbed(self, refName):
        '''
        load references from a tab delimited file. First column is gene name, second column is sequence
        '''
        
        logging.info('Loading and translating references from tabbed file: %s',refName)
        
        
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
        
        
        
        
class AaReference(Reference):
    '''
    classdocs
    '''


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
                #translate in 3 frames and add list to dict
                name = entrySplit[0]
                sequence = entrySplit[1]
                translations = translate_3_frames(sequence)
                self._refs[name] = translations
        
        #clean up
        logging.info('Loaded %d reference sequences',len(self._refs))
        refFile.close()


        
                