from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

'''
Created on 23 May 2012

@author: sh695
'''

def translate_3_frames(seq):
    translations =[]
    for p in range(0,3):
        dna=Seq(seq[p:].strip(),generic_dna)
        translations.append(dna.translate().tostring())
    return translations



def shorten_name(name): 
    #trim down geneName
    colonSplit = name.split(":")
    shorterName = "%s:%s"%(colonSplit[0],colonSplit[1]) 
    return shorterName