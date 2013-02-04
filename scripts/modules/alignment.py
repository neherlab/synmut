# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/11/12
content:    Tools for managing sequences and alignments of general nature.
'''
# Modules
import os
import glob
import numpy as np
from Bio import SeqIO, AlignIO
import general.classes.alphabet as abc


# Globals
datadir = '../data/'
genes = ['gag', 'pol', 'env', 'nef', 'vpu', 'vpr', 'vif', 'tat']



# Classes
class MSA(object):
    '''Multiple sequence alignment with general properties'''
    def __repr__(self):
        s = 'MSA(kind="'+self.kind+'", gene="'+self.gene+'"'
        if hasattr(self, 'reference'):
            s += ', reference="'+self.reference+'"'
        s += ')'
        return s

    def __init__(self, kind='subtypeB', reference=None, gene=None):

        kinddir = (datadir + kind).rstrip('/') + '/'
        self.kind = kind.rstrip('/')
        if reference is not None:
            kinddir = (kinddir + 'to_reference/')
            self.reference = reference
        if not os.path.isdir(kinddir):
            raise IOError('Directory not found')

        found = False
        for fn in glob.glob(kinddir+'*.fasta'):
            if gene in fn:
                MSA = AlignIO.read(fn, 'fasta')
                if reference is None:
                    self.MSA = MSA
                else:
                    self.reference_seq = MSA[0]
                    if reference not in self.reference_seq.id:
                        raise ValueError('Reference name not found in first sequence of the MSA')
                    self.MSA = MSA[1:]
                self.gene = gene
                found = True
                break
        if not found:
            raise IOError('Alignment file not found')


    def allele_frequencies(self, alpha=abc.alpha):
        '''Get the allele frequencies for an alphabet
        
        Parameters:
            - alpha: an alphabet
        '''
        alpha = list(alpha)
        L = len(self.MSA[0])
        afs = np.zeros((len(alpha), L))
        mat = np.array(self.MSA)
        for i, a in enumerate(alpha):
            afs[i] = (mat == a).sum(axis=0)
        afs /= afs.sum(axis=0)

        return afs




# Script
if __name__ == '__main__':

    msa = MSA('subtypeB', gene='nef')

