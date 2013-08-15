# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/03/2012
content:    Trajectories of allele frequencies.
'''
# Standard modules
from operator import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Bio.Seq import translate


# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
from modules.helper import is_nonsyn_table
from modules.alphabet import alpha
from conservation_syn_nonsyn_subtypeB import codon_single_mutants_synnonsyn



# Script
if __name__ == '__main__':

    # Define the patients
    patients = pS.parse_sequences(exclude=['p4', 'p7', 'p8', 'p11'])

    # Iterate over patients
    for k,p in enumerate(patients):
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
        is_nongap_cod = np.array(map(lambda x: '-' not in x, conscods))

        # Times (in months?)
        ts = np.array([p.visit_to_time(v) for v in p.visit]) * 30.5

        # Scan the genome codon-wise
        counts = np.zeros((len(ts), afs.shape[0] // 3), int)
        chances = np.zeros(afs.shape[0] // 3, int)
        # Proceed codon by codon
        for i in xrange(len(ts)):
            v = p.visit[i]
            seqs = np.array(p.seqs_from_visit(v))

            for j in xrange(afs.shape[0] // 3):
                # Exclude gaps
                if not is_nongap_cod[j]:
                    continue
    
                conscod = conscods[j]
                cods = map(''.join, seqs[:, 3 * j: 3* (j+1)])
                cods_sm = codon_single_mutants_synnonsyn(conscod)
    
                counts[i, j] += len([c for c in cods_sm['syn']
                                     if 0.25 <= 1.0 * cods.count(c) / len(cods) <= 0.75])
                if i == 0:
                    chances[j] += 1#len(cods_sm['syn'])

            
        is_good = is_nongap_cod & (chances > 0)
        p.diversity_25_75 = (1.0 * counts[:, is_good] / chances[is_good]).mean(axis=1)

        # Plot
        plt.figure()
        plt.plot(ts, p.diversity_25_75 * 3.0 / 3, c='r', lw=2, label='diversity 25-75')
        plt.xlabel('Time [days after SC]')
        plt.title(str(p), fontsize=20)

    # Poll results
    diversity_25_75 = np.mean([p.diversity_25_75.mean() * 3.0 / 3 for p in patients])
    print diversity_25_75
    

