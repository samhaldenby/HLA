'''
Created on 20 Jun 2012

@author: sh695
'''

class Status(object):
    '''
    Holds all statuses
    '''

    @staticmethod
    def Pass():
        return "P1"
    
    @staticmethod
    def Fail_low_score():
        return "F1"
    
    @staticmethod
    def Fail_score_diff():
        return "F2"
    
    @staticmethod
    def Warn_not_top_two():
        return "W1"
        