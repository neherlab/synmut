# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/03/2012
content:    Trajectories of allele frequencies.
'''
# Standard modules
from operator import *
import numpy as np
import matplotlib.pyplot as ppl
import matplotlib.cm as cm
from Bio.Seq import translate


# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
from modules.helper import is_nonsyn_table
from modules.alphabet import alpha



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
        is_nonsyn = is_nonsyn_table(consensus)

        # Scan the genome codon-wise
        counts = np.zeros(afs.shape[0], int)
        chances = np.zeros(afs.shape[0], int)
        for i in xrange(len(afs)):
            non_double_hit = ((afs[i - i % 3: i - i % 3 + 3] > 0).sum(axis=0).sum(axis=0) <= 4).all()
            if (not is_nongap[i]):# or (not non_double_hit):
                continue
            af = afs[i]
            for j, a in enumerate(alpha):
                # Gate for nonsynonymous changes
                if ((not is_nonsyn[i, j]) and       # Not nonsyn
                    (consensus[i] != a)):       # An actual mutation

                    chances[i] += 1
                    counts[i] += ((af[j] >= 0.25) & (af[j] <= 0.75)).sum()

        is_good = is_nongap & (chances > 0)
        p.diversity_25_75 = 1.0 * (counts[is_good]).sum() / (chances).sum() / (afs.shape[2])

        print p, p.diversity_25_75

    # Poll results
    print np.mean(map(attrgetter('diversity_25_75'), patients))
