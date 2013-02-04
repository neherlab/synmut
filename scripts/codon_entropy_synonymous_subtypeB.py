# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/12/12
content:    Measure the level of conservation of synonymous mutations in subtype
            B.
'''
# Modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Custom modules
import sys
sys.path.insert(0, '.')
import modules.alignment as alignment
from modules.alphabet import alpha
from modules.codon_usage import number_of_codons
from modules.helper import translate
from modules.parser_reference import seq_patterns, find_with_gaps



# Globals
genes = ['gag', 'pol', 'env', 'nef']
reference = 'SHAPE'
seq_patterns = seq_patterns[reference]



# Functions
def smoothen(x, f):
    l = len(x)
    x = x[:l - l % f]
    x = x.reshape((len(x) / f, f)).mean(axis=1)
    return x



# Script
if __name__ == '__main__':

    # Load data for a certain gene
    gene = genes[0]
    ali = alignment.MSA(gene=gene, reference=reference)
    msa = np.array(ali.MSA)
    afs = ali.allele_frequencies()

    # Filter only part of pol
    if gene == 'pol':
        msa = msa[:, :2718]
        afs = afs[:, :2718]

    # Find consensus sequence
    consind = afs.argmax(axis=0)
    consensus = alpha[consind]

    # Some seqs are not viable and include frameshifts that mess up the
    # translation, hence we restrict to positions for which gaps are a minority
    is_gap = consensus == '-'
    # Exclude full codons
    tmp = np.unique(is_gap.nonzero()[0] / 3)
    is_gap[np.concatenate([tmp * 3, tmp * 3 +1, tmp * 3 + 2])] = True

    # Exclude stop codons
    is_stop = np.zeros_like(is_gap)
    tmp = (translate(consensus) == '*').nonzero()[0]
    is_stop[np.concatenate([tmp * 3, tmp * 3 +1, tmp * 3 + 2])] = True

    ## Plot base prevalence
    #for i in xrange(4):
    #    plt.plot(np.arange(len(consensus)), afs[i],
    #             lw=1.5, alpha=0.5)
    #plt.xlim(0, 2600)
    #plt.ylim(-0.05, 1.25)
    #plt.xlabel('position in '+gene)
    #plt.ylabel('allele frequency')
    #plt.legend(alpha, loc=9)

    # Good are codons with no gaps and no stops
    is_good = (-is_gap) & (-is_stop)

    # For each codon, calculate the entropy
    msa = msa[:, is_good]
    consaa = translate(consensus[is_good])
    entropy = np.zeros(len(consaa))
    from collections import Counter
    for i, aa in enumerate(consaa):
        tmp = msa[:, i * 3: (i + 1) * 3]
        count = Counter(map(''.join, tmp))
        abus = []
        for (cod, abu) in count.iteritems():
            if translate(np.array(list(cod))) == aa:
                abus.append(abu)
        abus = np.array(abus)
        freqcod = 1.0 * abus / abus.sum()
        entropy[i] = -np.sum(freqcod * np.log(freqcod))

    # Maximal entropy per codon
    ncod = np.array(map(number_of_codons, consaa))
    entropymax = np.log(ncod)
    entropyfrac = entropy / (1e-5 + entropymax)

    # Plot the entropy
    plt.figure()
    plt.plot(np.arange(len(entropy)), entropyfrac, c='b', lw=1)
    plt.xlabel('Position in '+gene)
    plt.ylabel('Codon entropy / maximal entropy')
    plt.ylim(0, 1)

    # Plot smoothed curve
    smoothrange = 15
    xsmooth = smoothen(np.arange(len(entropy)), smoothrange)
    ysmooth = smoothen(entropyfrac, smoothrange)
    plt.plot(xsmooth, ysmooth, c='r', lw=2, ls='-')
    plt.title(gene, fontsize=18)

    # Mark V loops
    if gene == 'env':
        refstr = str(ali.reference_seq.seq)
        tmp = (-is_good).nonzero()[0]
        Vind = {}
        for key in ['V1', 'V2', 'V3', 'V4', 'V5']:
            i0 = find_with_gaps(refstr, seq_patterns[key][0])
            i0 -= (tmp < i0).sum()
            i0 //= 3
            i1 = find_with_gaps(refstr, seq_patterns[key][1])
            i1 -= (tmp < i1).sum()
            i1 //= 3
            Vind[key] = (i0, i1)
            h = Rectangle((i0, 0), i1 - i0, 1, color='yellow', alpha=0.5)
            plt.gca().add_patch(h)
            plt.text(i0 + 0.2 * (i1 - i0), 0.95, key, fontsize=12)
        

    plt.ion()
    plt.show()

    
