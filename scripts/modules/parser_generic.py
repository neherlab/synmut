# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/10/12
content:    Base class for patients and such with a few methods shared by all
            parsers.
'''
# Standard modules
import os
import re
import glob
from operator import *
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as ppl
import matplotlib.cm as cm

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import *



# Globals
alpha = np.array(['A', 'C', 'G', 'T', '-'])
alphal = list(alpha)



# Classes
class patient_generic(object):
    '''Generic class for a patient'''

    def __repr__(self):
        return self.name
    

    def L(self):
        return len(self.seqs[0])

   
    def allele_frequency_trajectories(self, first_visit=None, last_visit=None,
                                      alpha=alpha):
        '''Allele frequency trajectories.
        
        Parameters:
        -----------------
        first_visit: the first visit to track. Default is the first visit found.
        last_visit: the last visit to track. Default is the last visit found.
        alpha: dictionary for possible alleles. Default is usual DNA code.

        Returns:
        -----------------
        A tensor of size La x L x n_visits.
        '''
        if first_visit == None:
            first_visit = self.visit[0]
        if last_visit == None:
            last_visit = self.visit[-1]

        vs = self.visit[(self.visit >= first_visit) & (self.visit <= last_visit)]
        n = len(vs)
        L = self.L()
        La = len(alpha)
        aft = np.zeros((n, La, L))
        for i, v in enumerate(vs):
            aft[i] = self.allele_frequencies(seqs=self.seqs_from_visit(v),
                                             alpha=alpha)
        return aft.swapaxes(0,1).swapaxes(1,2)


    def allele_frequencies(self, seqs=None, alpha=alpha):
        '''Allele freqs in the choses seqs.

        Parameters:
        -----------------
        seqs: which seqs to calculate from. Default is all seqs of the patient.
        alpha: dictionary for possible alleles. Default is usual DNA code.
        '''
        if seqs is None:
            seqs = self.seqs
        smat = np.asarray(seqs)
        L = smat.shape[1]
        La = len(alpha)

        af = np.zeros((La,L))
        for i,a in enumerate(alpha):
            af[i] = (smat == a).mean(axis=0)
        return af


    def allele_count_trajectories(self, first_visit=None, last_visit=None,
                                      alpha=alpha):
        '''Allele count trajectories.
        
        Parameters:
        -----------------
        first_visit: the first visit to track. Default is the first visit found.
        last_visit: the last visit to track. Default is the last visit found.
        alpha: dictionary for possible alleles. Default is usual DNA code.

        Returns:
        -----------------
        A tensor of size La x L x n_visits.
        '''
        if first_visit == None:
            first_visit = self.visit[0]
        if last_visit == None:
            last_visit = self.visit[-1]

        vs = self.visit[(self.visit >= first_visit) & (self.visit <= last_visit)]
        n = len(vs)
        L = self.L()
        La = len(alpha)
        act = np.zeros((n, La, L), int)
        for i, v in enumerate(vs):
            act[i] = self.allele_counts(seqs=self.seqs_from_visit(v),
                                        alpha=alpha)
        return act.swapaxes(0,1).swapaxes(1,2)


    def allele_counts(self, seqs=None, alpha=alpha):
        '''Get the allele counts (useful for spotting singletons)'''
        if seqs is None:
            seqs = self.seqs
        smat = np.asarray(seqs)
        L = smat.shape[1]
        La = len(alpha)

        counts = np.zeros((La,L), int)
        for i,a in enumerate(alpha):
            counts[i] = (smat == a).sum(axis=0)
        return counts


    def set_consensus_initial(self):
        '''Set the consensus at the first time point.'''
        af = self.allele_frequencies(seqs=self.seqs_from_visit(self.visit[0]))
        self.consensus = alpha[af.argmax(axis=0)]

    
    def get_consensus_initial(self):
        '''Get the consensus at the first time point.'''
        if not hasattr(self, 'consensus'):
            self.set_consensus_initial()
        return self.consensus

