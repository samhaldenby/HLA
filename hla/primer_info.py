'''
Created on 25 Jul 2012

@author: sh695
'''

#Info on how reads should be trimmed
class PrimerInfo:
    name = ""
    trimFrom = None
    trimTo = None
    
    def __init__(self,name,trimFrom,trimTo):
        self.name = name
        self.trimFrom = trimFrom
        self.trimTo = trimTo