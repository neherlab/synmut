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


# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
from modules.helper import is_nonsyn_table



# Tables
alpha = pS.alpha



# Functions
def expanded_isnongap(afs):
    '''Restrict the sequence to non-gapped CODONS.'''
    is_gap = (afs[:,4,:].any(axis=1)).nonzero()[0]
    is_gap2 = set()
    for i in is_gap:
        is_gap2 |= set(range(i - i % 3, i - i % 3 + 3))
    is_nongap = np.ones(afs.shape[0], bool)
    is_nongap[list(is_gap2)] = False
    return is_nongap



# Script
if __name__ == '__main__':

    # Define the patients
    patients = pS.parse_sequences(reference='HXB2')

    for k,p in enumerate(patients[:1]):
        print p

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
        #is_nongap = expanded_isnongap(afs)
        is_nongap = (afs[:,4,:] == 0).all(axis=1)
        is_V_region = p.only_V_regions()

        # Divide into syn/nonsyn to spot hitchhiking
        consensus = alpha[afs[:,:,0].argmax(axis=1)]
        is_nonsyn = is_nonsyn_table(consensus)

        # Plot trajectories (no gaps, no conserved)
        figsyn = ppl.figure(figsize=[9.35, 6.325])
        fignonsyn = ppl.figure(figsize=[9.35, 6.325])
        axsyn = figsyn.add_subplot(1,1,1)
        axnonsyn = fignonsyn.add_subplot(1,1,1)
        x = p.mo_to_SC
        n_cols = afs.shape[0]
        cols = cm.jet([int(255.0 * i / n_cols) for i in xrange(n_cols)])

        for i in xrange(n_cols):
            # Exclude gaps and restrict to regions if requested
            if is_nongap[i]:# and not is_V_region[i]:
    
                # Do not plot consensus allele
                indmaj = afs[i,:,0].argmax()
                for a in xrange(4):
                    # Exclude initial consensus allele and singletons
                    if (a != indmaj) and ((afs[i,a] > 0).sum() > 1):
                        # Approximate months-to-days conversion
                        xx = 30.5 * (x + 1.5 / n_cols * i)
                        if is_nonsyn[i,a]:
                            axnonsyn.plot(xx, afs[i,a], lw=1.2, color=cols[i])
                        else:
                            axsyn.plot(xx, afs[i,a], lw=1.2, color=cols[i])

        for ax in [axsyn, axnonsyn]:
            ax.set_xlabel('Time from seroconversion [days]', fontsize=14)
            ax.set_ylabel(r'$\nu$', fontsize=18)
            ax.set_ylim(-0.05,1.05)

        axsyn.set_title('Derived synonymous allele frequency trajectories, '+str(p), fontsize=18)
        axnonsyn.set_title('Derived nonsynonymous allele frequency trajectories, '+str(p), fontsize=18)

    ppl.ion()
    ppl.show()
        

