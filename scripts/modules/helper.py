# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/03/12
content:    Helper functions for all parsers
'''
# Standard modules
import numpy as np
import matplotlib.pyplot as ppl
import Bio.Seq
import re


# Functions
def calculate_distance(smat0, smat1):
    dM = np.zeros((smat0.shape[0], smat1.shape[0]), int)
    for i0, s0 in enumerate(smat0):
        dM[i0] = (s0 != smat1).sum(axis=1)
    return dM


def plot_distance(dM, ax=None):
    if ax is not None:
        ppl.sca(ax)
    ppl.imshow(dM, interpolation='nearest')
    ppl.ylabel('First seqs')
    ppl.xlabel('Second seqs')
    ppl.title('Hamming Distance', fontsize=16)
    
    # set ticks to integers
    dMmax = dM.max()
    tickwidth = ((dMmax - 1) // 10) + 1
    tickmax = (dMmax // tickwidth + 1) * (tickwidth)
    ppl.colorbar(ticks=np.arange(tickmax))

    ppl.ion()
    ppl.show()


def translate(seqdna):
    '''Translate one or more sequences with gaps.'''
    if np.rank(seqdna) > 1:
        aa = []
        for seq in seqdna:
            aa.append(translate(seq))
        return np.asarray(aa)
    else:
        tmp = seqdna.reshape((len(seqdna) / 3, 3))
        aa = np.zeros(tmp.shape[0], str)
        for i1 in xrange(aa.shape[0]):
            if not (tmp[i1] == '-').any():
                aa[i1] = Bio.Seq.translate(''.join(tmp[i1]))
        return aa


def is_syn_table(sequence, *args, **kwargs):
    '''Find what loci would lead to synonymous mutations.'''
    return -is_nonsyn_table(sequence, *args, **kwargs)


def is_nonsyn_table(sequence, alpha=np.array(['A', 'C', 'G', 'T', '-'])):
    '''Find what loci would lead to nonsynonymous mutations.'''
    sequence = np.asarray(sequence)
    codons = sequence.reshape(sequence.shape[0] / 3, 3)
    aa = translate(sequence)

    # All possible single mutants
    is_nonsyn = np.ones((sequence.shape[0], 5), bool)
    for i in xrange(sequence.shape[0]):
        # Note: a gap is always nonsynonymous
        for a in xrange(4):
            tmpcod = codons[i // 3].copy()
            tmpcod[i % 3] = alpha[a]
            is_nonsyn[i,a] = aa[i // 3] != translate(tmpcod)[0]
    return is_nonsyn


def expand_gaps(is_nongap):
    '''Expand an array of indices for gaps to whole codons.'''
    is_nongap_exp = np.array(is_nongap, bool)
    i = 0
    while i < len(is_nongap):
        if not is_nongap[i]:
            is_nongap[i - i % 3: i - i % 3 + 3] = False 
            i += 3 - i % 3
        else:
            i += 1
    return is_nongap_exp


def find_with_gaps(seq, pattern):
    '''Find a pattern in a gapped sequence string'''
    sug = re.sub('-','', seq)
    index = sug.find(pattern)
    if index == -1:
        return index

    # Recover the right index in the gapped sequence
    iu = 0
    for ig, base in enumerate(seq):
        if iu == index:
            return ig
        if base != '-':
            iu += 1

    raise RuntimeError('I have trouble with finding the pattern')


def allele_entropy(af):
    '''Calculate the allele entropy'''
    ind = af > 0
    return np.dot(af[ind], -np.log(af[ind]))


