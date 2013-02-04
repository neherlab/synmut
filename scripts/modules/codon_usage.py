# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/10/12
content:    Read codon usage tables
'''
# Modules
import re
from longitudinal.classes.parser_generic import alpha
from Bio.Seq import translate



# Globals
datadir = '/home/fabio/university/phd/general/data/codon_usage/'
alpha = list(alpha)



# Classes
class codon_table(object):

    def __init__(self, organism='homo sapiens'):

        self.organism = organism
        self.datafile = datadir + re.sub(' ', '_', organism) +'.dat'
        
        table_within = {}
        table_abs = {}
        with open(self.datafile, 'r') as f:
            for line in f:
                if line[0] != '#':
                    fields = line.rstrip('\n').split('\t')
                    key = fields[0]
                    aa = fields[1]
                    freqwithin = float(fields[2])
                    freqabs = 1e-3 * float(fields[3])

                    table_within[key] = (aa, freqwithin)
                    table_abs[key] = freqabs

        self.within = table_within
        self.absolute = table_abs



# Functions
def all_codons():
    '''Generate a list with all codons'''
    return ('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC',
            'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT',
            'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC',
            'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
            'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
            'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT',
            'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
            'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
            'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC',
            'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT',
            'TTA', 'TTC', 'TTG', 'TTT')


def all_amino_acids():
    '''All amino acids'''
    return ('A', 'C', 'E', 'D', 'G', 'F',
            'I', 'H', 'K', '*', 'M', 'L',
            'N', 'Q', 'P', 'S', 'R', 'T',
            'W', 'V', 'Y')


def number_of_codons(amino_acid):
    '''Get the number of codons per amino acid'''
    d = {'L': 6, 'S': 6, 'R': 6, 'A': 4,
         'G': 4, 'P': 4, 'T': 4, 'V': 4,
         'I': 3, '*': 3, 'C': 2, 'E': 2,
         'D': 2, 'F': 2, 'H': 2, 'K': 2,
         'N': 2, 'Q': 2, 'Y': 2, 'M': 1,
         'W': 1}
    return d[amino_acid]


def volatility(codon):
    '''Number of nonsilent single mutants'''
    codon = ''.join(codon)
    if codon not in all_codons():
        raise ValueError('Codon not recognized. Remember we are in a DNA world.')

    aa = translate(codon)
    v = 0
    for pos in xrange(3):
        for m in alpha:
            if m not in ('-', codon[pos]):
                tmp = list(codon)
                tmp[pos] = m
                tmp = ''.join(tmp)
                v += translate(tmp) != aa
    return v


def effective_number_of_codons(codonlist):
    '''Calculate the ENC from a list of codons'''
    def F(n, codonlist):
        '''The probaibility that two codons are the same in the list'''
        import numpy as np
        from collections import Counter

        Ps = []
        for aa in all_amino_acids():
            if number_of_codons(aa) == n:
                cods = filter(lambda x: translate(x) == aa, codonlist)
                if cods:
                    P = 1.0 * np.array(Counter(cods).values()) / len(cods)
                    Ps.append(np.dot(P, P))
        if Ps:
            return np.mean(Ps)
        else:
            return None

    codonlist = map(''.join, codonlist)
    F2 = F(2, codonlist)
    F3 = F(3, codonlist)
    F4 = F(4, codonlist)
    F6 = F(6, codonlist)
    ENC = 2
    if F2 is not None:
        ENC += 9 / F2
    if F3 is not None:
        ENC += 1 / F3
    if F4 is not None:
        ENC += 5 / F4
    if F6 is not None:
        ENC += 3 / F6

    return ENC


                


# Script
if __name__ == '__main__':

    table = codon_table()
