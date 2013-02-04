# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/07/12
content:    Parser for Bunnik2008 (subparser of LANL)
'''
# Standard modules
from operator import *
import numpy as np

from Bio.Alphabet.IUPAC import *

from longitudinal.classes.parser_LANL import patient as patientLANL
from longitudinal.classes.parser_LANL import parse_sequences as parse_sequencesLANL
from longitudinal.classes.helper import translate, find_with_gaps
from longitudinal.classes.parser_reference import seq_patterns


# Globals
# Second rev exon in env (tat is included)
rev2 = {'ACH19999': ['ACCCACCTCCCA', 'TGCTGTTAGCTTG'],
        'ACH18969': ['ACCCGCCTCCCA', 'TGCTGTTAGCTTG'],
        'ACH19659': ['ACCCGCCTCCCA', 'TGCTGTTAGCTTG']}




# Classes
class patient(patientLANL):
    '''Specialized version of patient from parser_LANL.'''

    aaS = ['EEE', 'RPGGG']
    '''Shankarappa region'''

    def filter_for_longitudinal_genetics(self):
        '''Filter only time points with sequences.'''
        self.filter_seqs_with_attribute('days_from_seroconversion')


    def only_Shanka(self):
        '''Denote the Shankarappa region.
 
        Returns:
            index: boolean index with the Shankarappa region marked by 'True's.

        *Note*: The region is marked by 'EEE' at the beginning and 'RPGGG' at the end.
        '''

        cons = self.get_consensus_initial()
        is_Shanka = np.zeros(len(cons), bool)
        consaa = ''.join(translate(cons))
        is_Shanka[consaa.find(self.aaS[0]) * 3: consaa.find(self.aaS[1]) * 3] = True
        return is_Shanka


    def only_rev2(self):
        '''Denote the second rev exon.

        Returns:
            index: boolean index with the Shankarappa region marked by 'True's.
        '''
        cons = ''.join(self.get_consensus_initial())
        is_rev2 = np.zeros(len(cons), bool)
        is_rev2[cons.find(rev2[self.name][0]): cons.find(rev2[self.name][1])] = True
        return is_rev2


    def find_V_region(self, name):
        '''Find a variable region based on the reference'''
        if name not in self.seq_patterns:
            raise ValueError(str(name)+" is not in self.seq_patterns.")
        
        patterns = self.seq_patterns[name]
        seq = str(self.reference_seq.seq)
        start = find_with_gaps(seq, patterns[0])
        end = find_with_gaps(seq[start:], patterns[1])
        if end != -1:
            end += start 
        return [start, end]


    def show_V_region(self, name):
        '''Show a variable region'''
        indices = self.find_V_region(name)
        return self.reference_seq[indices[0]: indices[1]]


    def only_V_regions(self):
        '''Denote the hypervariable regions'''

        cons = self.get_consensus_initial()
        is_V = np.zeros(len(cons), bool)
        for name in self.seq_patterns:
            indices = self.find_V_region(name)
            is_V[indices[0]:indices[1]] = True
        return is_V



# Functions
def parse_sequences(*args, **kwargs):
    '''Parse refined sequences.'''
    filename = 493
    patients = parse_sequencesLANL(filename, patientclass=patient, *args, **kwargs)
    for p in patients:
        if hasattr(p, 'reference_name') and p.reference_name in seq_patterns:
            p.seq_patterns = seq_patterns[p.reference_name]
            '''Variable regions in env'''
    return patients

