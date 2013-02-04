# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/09/12
content:    Implement some globals used by many parsers that interact with a
            reference.
'''
# Modules
import numpy as np
import Bio.AlignIO as AlignIO
from itertools import permutations

from longitudinal.classes.helper import find_with_gaps
from longitudinal.classes.parser_generic import alpha



# Globals
# these are derived from each other using MSA.corresponding_pattern
# (and ultimately by hand from the horrible XLS file for HXB2)
seq_patterns = {'SHAPE':{'V1': ['TGCACTGATTTGAA', 'TGCTCTTTCAATATCAGC'],
                         'V2': ['TGCTCTTTCAATATCAGC', 'TGTAACACCTCAGTCAT'],
                         'V3': ['TGTACAAGACCCAACAACA', 'TGTAACATTAGTAGA'],
                         'V4': ['TGTAATTCAACACAACTGTTTAA', 'TGCAGAATAAAACAATT'],
                         'V5': ['TAATAACAACAATGGGTCCGAGAT', 'CCTGGAGGA'],
                         'env_signal': ['ATGAGAGTGAAGGAGAA', 'ACAGAAAAATTGTGGGTC']},
                'NL4-3':{'V1': ['TGCACTGATTTGAA', 'TGCTCTTTCAATATCAGC'],
                         'V2': ['TGCTCTTTCAATATCAGC', 'TGTAACACCTCAGTCAT'],
                         'V3': ['TGTACAAGACCCAACAACA', 'TGTAACATTAGTAGA'],
                         'V4': ['TGTAATTCAACACAACTGTTTAA', 'TGCAGAATAAAACAATT'],
                         'V5': ['TAATAACAACAATGGGTCCGAGAT', 'CCTGGAGGA']},
                'HXB2': {'V1': ['TGCACTGATTTGAA', 'TGCTCTTTCAATATCAGC'],
                         'V2': ['TGCTCTTTCAATATCAGC', 'TGTAACACCTCAGTCAT'],
                         'V3': ['TGTACAAGACCCAACAACA', 'TGTAACATTAGTAGA'],
                         'V4': ['TGTAATTCAACACAACTGTTTAA', 'TGCAGAATAAAACAAAT'],
                         'V5': ['TAATAGCAACAATGAGTCCGAGAT', 'CCTGGAGGA'],
                         'env_signal': ['ATGAGAGTGAAGGAGAA', 'ACAGAAAAATTGTGGGTC']}
               }

genes = {'SHAPE': {'gag': ['ATGGGTGCGAGAGCGTCGGTA', 'AGATAGGGGGGCAATTAA'],
                   'pol': ['CCTCAGATCACTCTTTGGCAGCGACCCC', 'CACATGGAAAAGATTAGTA'],
                   'env': ['ATGAGAGTGAAGGAG', 'GATGGGTGGCAAGTG'],
                   'nef': ['ATGGGTGGCAAGTG', 'CATCGAGCTTGCTACAAGGG']},
         'HXB2': {'gag': ['ATGGGTGCGAGAGCGTCAGTA', 'AGATAGGGGGGCAACTAA'],
                  'pol': ['CCTCAGGTCACTCTTTGGCAACGACCCC', 'AACATGGAAAAGTTTAGTA'],
                  'env': ['ATGAGAGTGAAGGAG', 'GATGGGTGGCAAGTG'],
                  'nef': ['ATGGGTGGCAAGTG', 'CATCGAGCTTGCTACAAG']},
        }

datadir = '/home/fabio/university/phd/longitudinal/data/reference'



class MSA(object):
    '''An MSA of the given references.

    Note: the MSA must already exist, we are not calling MUSCLE here for reasons
    of clarity (although we could).
    '''

    def __init__(self, references):
        '''Initialize'''
        for tmplist in permutations(references):
    
            datafile = datadir+'/'+'_'.join(tmplist)
            if len(references) > 1:
                datafile = datafile+'_aligned'
            datafile = datafile+'.fasta'
    
            try:
                self.MSA = AlignIO.read(datafile, format='fasta')
                break
            except IOError:
                pass
    
        if not hasattr(self, 'MSA'):
            raise IOError('Reference alignment not found')


    def __repr__(self):
        return self.MSA.__repr__()


    def __str__(self):
        return self.MSA.__str__()


    def corresponding_pattern(self, pattern, origin, target):
        L = len(pattern)
        fo = lambda x: origin in self.MSA[x].name
        ft = lambda x: target in self.MSA[x].name
        iseqo = filter(fo, xrange(len(self.MSA)))[0]
        iseqt = filter(ft, xrange(len(self.MSA)))[0]

        index = find_with_gaps(str(self.MSA[iseqo].seq), pattern)
        if index == (-1):
            raise ValueError('Pattern not found in origin')

        tmplist = []
        while (len(tmplist) < L) and (index < self.MSA.get_alignment_length()):
            base = self.MSA[iseqt, index]
            if base != '-':
                tmplist.append(base)
            index += 1

        return ''.join(tmplist)


    
