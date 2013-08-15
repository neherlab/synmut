# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/03/2012
content:    Trajectories of allele frequencies.
'''
# Standard modules
from operator import *
import numpy as np
import matplotlib.cm as cm
from Bio.Seq import translate
import matplotlib.pyplot as plt


# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
from modules.helper import is_nonsyn_table
from modules.alphabet import alpha
from conservation_syn_nonsyn_subtypeB import codon_single_mutants_synnonsyn


def get_is_mutation(consensus):
    is_mutation = np.ones((len(alpha), len(consensus)), bool)
    alphal = list(alpha)
    for i, a in enumerate(consensus):
        is_mutation[alphal.index(a), i] = False
    return is_mutation.T



# Script
if __name__ == '__main__':

    # Define the patients
    patients = pS.parse_sequences(exclude=['p4', 'p7', 'p8', 'p11'])

    # Iterate over patients
    for k,p in enumerate(patients[:1]):
        p.filter_only_sequenced()

        # Measure the allele frequencies (of all alleles at all positions)
        paf = p.allele_frequencies
        afs = np.asarray([paf(seqs=p.seqs_from_visit(v)) for v in p.visit])

        # Reshape so that we get the site as first axis, the nucleotide as
        # second, the time as third
        afs = afs.swapaxes(0,2)

        # Eliminate gaps (whole codons are excluded to keep translation possible)
        # first and translate then!
        # Note: this has been implemented at the root
        is_nongap = (afs[:,4,:] == 0).all(axis=1)

        # Divide into syn/nonsyn to spot hitchhiking
        consensus = alpha[afs[:,:,0].argmax(axis=1)]
        conscods = map(''.join, consensus.reshape((len(consensus) / 3, 3)))
        is_mut = get_is_mutation(consensus)
        is_nonsyn = is_mut & is_nonsyn_table(consensus)
        is_syn = is_mut & (-is_nonsyn)

        # Times (in months?)
        ts = np.array([p.visit_to_time(v) for v in p.visit]) * 30.5

        # Divergence
        divergence = {key: np.zeros(afs.shape[-1]) for key in ['syn', 'nonsyn']}
        diversity = {key: np.zeros(afs.shape[-1]) for key in ['syn', 'nonsyn']}
        for i in xrange(len(ts)):
            v = p.visit[i]
            seqs = np.array(p.seqs_from_visit(v))
           
            counts = {key: 0 for key in ['syn', 'nonsyn']}
            chances = {key: 0 for key in ['syn', 'nonsyn']}
            counts2 = {key: 0 for key in ['syn', 'nonsyn']}
            chances2 = {key: 0 for key in ['syn', 'nonsyn']}
            # Proceed codon by codon
            for j in xrange(afs.shape[0] // 3):
                # Exclude gaps
                if not is_nongap[3 * j: 3* (j+1)].all():
                    continue

                conscod = conscods[j]
                cods = map(''.join, seqs[:, 3 * j: 3* (j+1)])               
                cods_sm = codon_single_mutants_synnonsyn(conscod)

                # Divergence
                # Restrict to 4-fold degenerate for normalization issues?
                for key in divergence.iterkeys():
                    counts[key] += sum([cods.count(c) for c in cods_sm[key]])
                    chances[key] += len(seqs) * len(cods_sm[key])

                # Diversity
                for ci, cod in enumerate(cods):
                    cods_sm = codon_single_mutants_synnonsyn(cod)
                    for key in diversity.iterkeys():
                        counts2[key] += len([c for c in cods[:ci] if c in cods_sm[key]])
                        chances2[key] += len(cods_sm[key]) * ci

            # Normalize
            for key in divergence.iterkeys():
                divergence[key][i] = 1.0 * counts[key] / chances[key]
                diversity[key][i] = 1.0 * counts2[key] / chances2[key]


        # Plot
        plt.figure()
        plt.plot(ts, divergence['syn'] * 9.0 / 3, c='r', lw=2, label='divergence syn')
        plt.plot(ts, divergence['nonsyn'] * 9.0 / 3, c='b', lw=2, label='divergence nonsyn')
        plt.plot(ts, (divergence['syn'] + divergence['nonsyn']) * 9.0 / 3,
                 c='y', lw=2, label='divergence all')
        plt.plot(ts, diversity['syn'] * 9.0 / 3, c='g', lw=2, label='diversity syn')
        plt.plot(ts, diversity['nonsyn'] * 9.0 / 3, c='k', lw=2, label='diversity nonsyn')
        plt.plot(ts, (diversity['syn'] + diversity['nonsyn']) * 9.0 / 3,
                 c='orange', lw=2, label='diversity all')
        plt.xlabel('Time [days after SC]')
        plt.ylabel('div/div')
        plt.title(str(p), fontsize=20)
        plt.legend(loc=2)

        # Save figure
        plt.savefig('../figures/Shankarappa_divergence_diversity_'+str(p)+'.pdf')

    plt.ion()
    plt.show()

