# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/12/12
content:    Measure the level of conservation of synonymous mutations in subtype
            B.
'''
# Modules
import re
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

    # Store data
    data = {}

    # Make figure for all four genes
    fig, axs = plt.subplots(2, 2)
    axs = axs.flatten()

    for k, gene in enumerate(genes):

        print gene

        # Select the axis
        ax = axs[k]
        plt.sca(ax)

        # Load data for a certain gene
        ali = alignment.MSA(gene=gene, reference=reference)
        if gene in data:
            msa = data[gene]
        else:
            msa = np.array(ali.MSA)
            data[gene] = msa
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
            abundances = []
            for (cod, abundance) in count.iteritems():
                if bool(re.findall('[ACGT]{3}', cod)) and (translate(np.array(list(cod))) == aa):
                    abundances.append(abundance)
            abundances = np.array(abundances)
            freqcod = 1.0 * abundances / abundances.sum()
            entropy[i] = -np.sum(freqcod * np.log(freqcod))
    
        # Maximal entropy per codon
        ncod = np.array(map(number_of_codons, consaa))
        entropymax = np.log(ncod)
        entropyfrac = entropy / (1e-5 + entropymax)
    
        # Plot the entropy
        plt.plot(np.arange(len(entropy)), entropyfrac, c='b', lw=1)
        plt.xlabel('Position in '+gene)
        plt.ylim(0, 1)
        if k in [0, 2]:
            plt.ylabel('Codon entropy / maximal entropy')
        else:
            ax.set_yticklabels([''])
    
        # Plot smoothed curve
        smoothrange = 15
        xsmooth = smoothen(np.arange(len(entropy)), smoothrange)
        ysmooth = smoothen(entropyfrac, smoothrange)
        plt.plot(xsmooth, ysmooth, c='r', lw=2, ls='-')
        plt.title(gene, fontsize=18)
    
        # Mark V loops and rev 2nd exon
        if gene == 'env':
            # V loops
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

            # Rev exon
            s = ''.join(consensus[is_good])
            i0 = s.find('ACCCGCCTCCCA') // 3
            i1 = s.find('TGCTGTTAGCTTG') // 3
            h = Rectangle((i0, 0), i1 - i0, 1, color='magenta', alpha=0.5)
            plt.gca().add_patch(h)
            t = plt.text(i1 + 0.15 * (i1 - i0), 0.85, 'rev\n2nd\nexon', fontsize=12)

        

    plt.ion()
    plt.show()

    
