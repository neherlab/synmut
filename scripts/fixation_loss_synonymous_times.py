# vim: fdm=indent
'''
author:     Fabio Zanini
date:       15/06/12
content:    Calculate how many of the synonymous alleles at intermediate nu0
            reach either boundary after a certain time.
'''
# Standard modules
from sys import stderr
import re
from operator import *
import numpy as np
import matplotlib.pyplot as ppl
import matplotlib.cm as cm
import matplotlib.patches as patches

# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
from modules.helper import translate
from modules.alphabet import alpha

# Tables
# Verbosity level (0-3)
VERBOSE = 1

# Bad patients (selected manually in a conservative way looking e.g. at PCA)
bad_patients = ['p4', 'p7', 'p8', 'p11']

# Initial frequency bins
nu0ss = [[0.15, 0.25],
         [0.25, 0.4],
         [0.4, 0.6]]

# Colors for the plot
colors = ['b', 'g', 'r']

# Maximal time to plot
maxt = 40 # Months

# Approximate months-to-days conversion
mo_to_gen = 30.5



# Script
if __name__ == '__main__':

    fixed = []
    lost = []
    fixlost = [lost, fixed]

    # Counts and times (overall)
    # The first list of each nu0s is lost, the second fixed, the third floating
    counts_all = [[0, 0, 0] for nu0s in nu0ss]
    # The first list of each nu0s is lost, the second fixed
    times_all = [[[], []] for nu0s in nu0ss]
    
    # Define the patients
    patients = pS.parse_sequences(exclude=bad_patients)
    
    # Aggregate information from all patients
    for k, p in enumerate(patients):

        if VERBOSE >= 1:
            stderr.write(str(p)+'\n')

        # Filter the time points to sequenced times
        p.filter_only_sequenced()
        n = len(p.visit)
    
        # Filter away gaps (conservative), keep reading frame
        is_nongap = ((np.array(p.seqs) == '-').sum(axis=0) == 0)
        for i in xrange(len(is_nongap)):
            if not is_nongap[i]:
                is_nongap[i - i % 3: i - i % 3 + 3] = False
    
        # Collect allele frequencies
        afs = p.allele_frequency_trajectories().swapaxes(1,2).swapaxes(0,1)
        # Consensus at the initial time point to decide on
        # mutation and synonymity
        consensus = p.get_consensus_initial()
    
        # Iterate over starting frequency windows
        for iji, nu0s in enumerate(nu0ss):
            if VERBOSE >= 2:
                stderr.write('nu0s = ['+str(nu0s[0])+', '+str(nu0s[1])+']\n')
    
            # Filter away already used alleles (we do it within each frequency
            # bin, hence it is not maximally conservative in terms of confidence
            # levels -- but we do not calculate any of those, and a small degree
            # of double counting is fine with our scarse data points, since it
            # should not bias the results)
            is_nonused = np.ones(p.L(), bool)
    
            # Define the counts (lost, fixed, neither)
            counts = counts_all[iji]
            delta_ts = times_all[iji]
    
            # Initial time point can be everything up to the second-last visit
            for i in xrange(n-1):

                if VERBOSE >= 2:
                    stderr.write('v'+str(p.visit[i])+': ')
    
                # Get allele frequencies
                af0 = afs[i]
    
                # Filter only alleles in the nu0 range and have never gaps
                ind_nu0 = (af0 > nu0s[0]) & (af0 <= nu0s[1]) & is_nongap
                alleles = np.transpose(ind_nu0.nonzero())
    
                # Filter only actual mutations (at t0) from consensus
                if len(alleles):
                    is_mut = np.array([alpha[al[0]] != consensus[al[1]] for al in alleles], bool)
                    alleles = alleles[is_mut]
    
                # filter only synonymous changes
                if len(alleles):
                    is_syn = np.zeros(len(alleles), bool)
                    for j, al in enumerate(alleles):
                        pos = al[1]
                        mut = alpha[al[0]]
                        codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                        cod = codcons.copy()
                        cod[pos % 3] = mut
                        is_syn[j] = (translate(codcons) == translate(cod))
                    alleles = alleles[is_syn]
    
                # Test more stringently for synonymity
                # (we want to avoid double-hits in one single codon)
                if len(alleles):
                    seqs = np.array(p.seqs_from_visit(p.visit[i]))
                    is_single = np.zeros(len(alleles), bool)
                    for j, al in enumerate(alleles):
                        pos = al[1]
                        mut = alpha[al[0]]
                        # Check whether sequences have double-hits
                        # If a double mutant is *ever* observed, discard the allele.
                        # Note: This is a very conservative measure and must
                        # be avoided when estimating densities.
                        codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                        seqal = seqs[seqs[:,pos] == mut]
                        codseq = seqal[:, pos - pos % 3: pos - pos % 3 + 3]
                        is_single[j] = ((codcons != codseq).sum(axis=1) == 1).all()
                        if (VERBOSE > 2) and (is_single[j]):
                            stderr.write(''.join(codcons)+' --> '+''.join(codseq[0])+'\n')
    
                    alleles = alleles[is_single]
                    if VERBOSE > 1:
                        stderr.write(str(len(alleles))+' single mutant alleles in the right range\n')
            
                # Check whether they fix/are lost in at least one time point
                if len(alleles):
                    for al in alleles:
                        alf_sel0 = af0[al[0], al[1]]
    
                        for l in xrange(i+1, n):
                            alf_sel1 = afs[l][al[0], al[1]]
            
                            # Increment the counts of lost/fixed
                            for j in xrange(2):
                                # Test whether the locus has been used already
                                if (alf_sel1 == j) and (is_nonused[al[1]]):

                                    # Finally: set the counts
                                    counts[j] += 1
                                    delta_ts[j].append(p.mo_to_SC[l]-p.mo_to_SC[i])
                                    is_nonused[al[1]] = False

                                    # Add the alleles to the list of dicts
                                    fixlost[j].append({'patient': str(p),
                                                       'pos': al[1],
                                                       'all': alpha[al[0]],
                                                       'all_cons': consensus[al[1]],
                                                       })
    
                                    # Print verobsely
                                    if VERBOSE > 1:
                                        stderr.write(str(p))
                                        stderr.write(', delta vi = '+str(l-i))
                                        stderr.write(', delta t = '+str(p.mo_to_SC[l]-p.mo_to_SC[i]))
                                        stderr.write(', nu0 = ')
                                        stderr.write('{:1.2f}'.format(alf_sel0))
                                        stderr.write(', nu1 = '+str(j)+'\n')
    
                        # If the allele has neither fixed not vanished by the end,
                        # count it in the third bucket
                        # Note: exclude late ones which have not been around for
                        # long enough
                        if is_nonused[al[1]] and (p.mo_to_SC[-1]-p.mo_to_SC[i] > maxt):
                            if VERBOSE >= 3:
                                stderr.write(str(p)+', pos '+str(al[1])+'\n')
                            counts[2] += 1
                            is_nonused[al[1]] = False
    
    
                if VERBOSE and (i == n-2):
                    stderr.write('\n')
        if VERBOSE >= 1:
            stderr.write('\n')
    
    # Plot fixation/extinction times
    fig = ppl.figure()
    ax = fig.add_subplot(1,1,1) 
    for iji, nu0s in enumerate(nu0ss):
    
        # Define the counts (lost, fixed, neither)
        counts = counts_all[iji]
        delta_ts = times_all[iji]
        P0 = 1.0 * counts[0] / sum(counts)
        P1 = 1.0 * counts[1] / sum(counts)
    
        # Plot the fixation/loss probability together as a function of time,
        # by simply counting how many alleles (in percentage) had fixed/been lost after time dt
        totcounts = sum(counts)
        dt_fix = np.sort(delta_ts[1]) * mo_to_gen
        dt_loss = np.sort(delta_ts[0]) * mo_to_gen
        nu0m = sum(nu0s) / len(nu0s)
    
        # Merge the extinction/fixation curves
        x = np.concatenate([dt_fix, dt_loss[::-1]])
        y = np.concatenate([np.linspace(0,P1,len(dt_fix)),
                            (1 - np.linspace(0,P0,len(dt_loss)))[::-1]])
        ax.plot(x, y,
                lw=2, c=colors[iji],
                label=r'$\nu_0 = ['+str(nu0s[0])+','+str(nu0s[1])+']$')
    
        # Draw rectangle
        rec = patches.Rectangle((0,nu0s[0]), 2000, nu0s[1] - nu0s[0],
                                color=colors[iji], alpha=0.3)
        ax.add_patch(rec)

    ax.set_xlabel(r'$t_\text{boundary} - t_i$ [days]', fontsize=16)
    ax.set_ylabel(r'$\nu$', fontsize=18)
    ax.set_title('Fixation/extinction times')
    ax.legend(loc=1)
    ax.set_xlim(0, 2000)

    ppl.ion()
    ppl.show()
