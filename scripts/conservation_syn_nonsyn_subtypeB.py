# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/12/12
content:    Measure the level of conservation of synonymous and nonsynonymous 
            mutations in subtype B. Nonsynonymous mutations are expected to be
            more broadly conserved, because they impair protein function.
'''
# Modules
import sys
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Custom modules
import modules.alignment as alignment
from modules.alphabet import alpha, alphal_ng
from modules.codon_usage import all_codons
from modules.codon_usage import translate as transstrict
from modules.helper import translate
from modules.parser_reference import seq_patterns, find_with_gaps



# Globals
gene = 'gag'
reference = 'SHAPE'
seq_patterns = seq_patterns[reference]



# Functions
def smooth(x, f):
    l = len(x)
    x = x[:l - l % f]
    x = x.reshape((len(x) / f, f)).mean(axis=1)
    return x


def codon_single_mutants_synnonsyn(codon):
    codon = list(codon)
    aa = transstrict(''.join(codon))
    syn = []
    nonsyn = []
    for pos in xrange(3):
        for mut in alphal_ng:
            if mut != codon[pos]:
                tmp = list(codon)
                tmp[pos] = mut
                tmp = ''.join(tmp)
                if transstrict(tmp) == aa:
                    syn.append(tmp)
                else:
                    nonsyn.append(tmp)
    return {'syn': syn, 'nonsyn': nonsyn}




# Script
if __name__ == '__main__':

    # Load data for a certain gene
    ali = alignment.MSA(gene=gene, reference=reference)
    msaraw = np.array(ali.MSA)
    afs = ali.allele_frequencies()

    # Filter only part of pol
    if gene == 'pol':
        msaraw = msaraw[:, :2718]
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

    # Get all syn and nonsyn single mutants of all codons
    sms = {cod: codon_single_mutants_synnonsyn(cod) for cod in all_codons()}

    # For each codon, calculate the number of syn and nonsyn changes from the
    # consensus
    msa = msaraw[:, is_good]
    conscods = consensus[is_good].reshape((is_good.sum() / 3, 3))
    # Count the number of sequences that enter the statistics for each codon, to
    # counterbalance for ambiguous sequences
    nseqs = msa.shape[0] * np.ones(len(conscods))
    # Prepare the structure to save the counts (columns: not changed, syn, nonsyn)
    counts3 = np.zeros((len(conscods), 3), int)
    
    # Run over consensus codons and count
    for i, conscod in enumerate(conscods):
        conscod = ''.join(conscod)
        smscod = sms[conscod]
        codsall = map(''.join, msa[:, i * 3: (i + 1) * 3])
        counts = Counter(codsall)
        for (cod, abundance) in counts.iteritems():
            if cod == conscod:
                counts3[i, 0] += abundance
            elif cod in smscod['syn']:
                counts3[i, 1] += abundance
            elif cod in smscod['nonsyn']:
                counts3[i, 2] += abundance

    # Calculate the number of syn and nonsyn for each codon along the consensus
    nsmscons = {'syn': [],
                'nonsyn': []}
    for cod in conscods:
        tmp = sms[''.join(cod)]
        for key in ['syn', 'nonsyn']:
            nsmscons[key].append(len(tmp[key]))
    for key in ['syn', 'nonsyn']:
        nsmscons[key] = np.array(nsmscons[key])

    # Take e.g. four-fold degenerate codons (as SMs only) and calculate frequencies
    is_deg = nsmscons['syn'] == 3
    nsmscons2 = np.vstack([nsmscons['syn'], nsmscons['nonsyn']]).T
    freqs2 = 1.0 * counts3[:, 1:][is_deg] / nsmscons2[is_deg] / msa.shape[0]

    ## Plot a histogram of conservation levels
    #plt.figure()
    #plt.hist(freqs2,color=('r', 'b'), bins=np.linspace(0, 0.04, 20),
    #         normed=True, label=('syn', 'nonsyn'))
    #plt.xlabel('diversity')
    #plt.ylabel('histogram')
    #plt.title('diversity at 4-fold degen codons, '+gene)

    # Plot cumulative
    freqsyn = np.sort(freqs2[:,0])
    freqnonsyn = np.sort(freqs2[:,1])
    fig = plt.figure(figsize=[6,4.6])
    plt.plot(freqsyn, np.linspace(0, 1, len(freqsyn)), c='r', lw=2, label='syn')
    plt.plot(freqnonsyn, np.linspace(0, 1, len(freqnonsyn)), c='b', lw=2, label='nonsyn')
    plt.xlabel('diversity')
    plt.ylabel('cumulative')
    plt.xlim(0, 0.04)
    plt.legend(loc=4)
    plt.title('diversity at 4-fold degen codons, '+gene)
    line = plt.plot([0.003] * 2, [0.08, 0.60], c='k', lw=1.5, ls='--')
    txt = plt.text(0.0035, 0.5, r'$\Delta \approx 0.52$', fontsize=14)

    plt.ion()
    plt.show()
