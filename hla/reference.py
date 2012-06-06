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
        
        saves the reference data in dict of type {sequence : count}
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
                
                #abbreviate name and add 
                shorterName = shorten_name(name)
                
                if shorterName not in self._refs:
                    self._refs[shorterName] = []
                self._refs[shorterName].append(sequence)
                
                    
                     
        
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
                
                #shorten name
                shorterName = shorten_name(name)
                #debugging code
                if shorterName in self._refs:
                    #compare translations to ensure they're the same
                    matches = 0
                    for q in translations:
                        for r in self._refs[shorterName]:
                            if q in r or r in q:
                                matches+=1
                   # print "%s -> %s already exists! %d translations match"%(name, shorterName, matches)
                    
                self._refs[shorterName] = translations
        
        #clean up
        logging.info('Loaded %d reference sequences',len(self._refs))
        refFile.close()


        
                