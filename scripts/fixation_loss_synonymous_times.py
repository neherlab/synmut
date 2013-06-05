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
from matplotlib import rcParams
params = {'backend': 'ps',
        'axes.labelsize': 20, 
        'text.fontsize': 20,
    'legend.fontsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'text.usetex': True}
    
rcParams.update(params)
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

    fixlost = {'syn': [[], []], 'nonsyn': [[], []]}

    # Counts and times (overall)
    # The first list of each nu0s is lost, the second fixed, the third floating
    counts_all = {key: [[0, 0, 0] for nu0s in nu0ss] for key in ['syn', 'nonsyn']}
    # The first list of each nu0s is lost, the second fixed
    times_all = {key: [[[], []] for nu0s in nu0ss] for key in ['syn', 'nonsyn']}
    
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

                if not len(alleles):
                    continue
    
                ###############################################################
                # Filter only synonymous/nonsynonymous changes
                ###############################################################
                # First test
                is_syn = np.zeros(len(alleles), bool)
                for j, al in enumerate(alleles):
                    pos = al[1]
                    mut = alpha[al[0]]
                    codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                    cod = codcons.copy()
                    cod[pos % 3] = mut
                    is_syn[j] = (translate(codcons) == translate(cod))

                alleles_syn = alleles[is_syn]
                alleles_nonsyn = alleles[-is_syn]
    
                # Test more stringently (avoid double-hits in one single codon)
                if len(alleles_syn):
                    seqs = np.array(p.seqs_from_visit(p.visit[i]))
                    is_single = np.zeros(len(alleles_syn), bool)
                    for j, al in enumerate(alleles_syn):
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

                    alleles_syn = alleles_syn[is_single]

                # Print some info
                if VERBOSE > 1:
                    stderr.write(str(len(alleles_syn))+' single mutant syn alleles and ')
                    stderr.write(str(len(alleles_nonsyn))+' single mutant nonsyn alleles in the right range\n')

                alleles_all = {'syn': alleles_syn, 'nonsyn': alleles_nonsyn}
                ###############################################################
                # End of synonymous test
                ###############################################################
            
                # Check whether they fix/are lost in at least one later time point
                for (key, alleles) in alleles_all.iteritems():
                    # Define the counts (lost, fixed, neither)
                    counts = counts_all[key][iji]
                    delta_ts = times_all[key][iji]

                    for al in alleles:
                        # Initial allele frequency at t0
                        alf_sel0 = af0[al[0], al[1]]
    
                        # Look at later time points
                        for l in xrange(i+1, n):

                            # Later allele frequency at t1
                            alf_sel1 = afs[l][al[0], al[1]]
            
                            # Look at LOST (j=0) and FIXED (j=1)
                            for j in xrange(2):
                                # Test whether the locus has been used already
                                if (np.abs(alf_sel1 - j) < 1e-2) and (is_nonused[al[1]]):

                                    # Finally: set the counts
                                    counts[j] += 1
                                    delta_ts[j].append(p.mo_to_SC[l]-p.mo_to_SC[i])

                                    # The position is used now
                                    is_nonused[al[1]] = False

                                    # Add the alleles to the list of dicts
                                    newitem = {'patient': str(p),
                                               'pos': al[1],
                                               'all': alpha[al[0]],
                                               'all_cons': consensus[al[1]],
                                               'cod_cons': ''.join(consensus[al[1] - al[1] % 3:
                                                                             al[1] - al[1] % 3 + 3]),
                                               'codpos': al[1] % 3,
                                               'time': p.mo_to_SC[l]-p.mo_to_SC[i],
                                               'type': key,
                                               'nu0': alf_sel0,
                                               'nu1': alf_sel1,
                                               't0': p.mo_to_SC[i],
                                               't1': p.mo_to_SC[l],
                                               'nu0_index': iji,
                                                }
                                    fixlost[key][j].append(newitem)
    
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
    
    
                if (VERBOSE > 1) and (i == n-2):
                    stderr.write('\n')
        if VERBOSE > 1:
            stderr.write('\n')
    
    # Plot fraction of surviving mutations over time
    fig, axs = ppl.subplots(1, 2, figsize=(14, 6))
    for i, key in enumerate(['syn', 'nonsyn']):
        ax = axs[i]
        for iji, nu0s in enumerate(nu0ss):
 
            # Define the counts (lost, fixed, neither)
            counts = counts_all[key][iji]
            delta_ts = times_all[key][iji]
            P0 = 1.0 * counts[0] / sum(counts)
            P1 = 1.0 * counts[1] / sum(counts)
        
            # Plot the fixation/loss probability together as a function of time,
            # by simply counting how many alleles (in percentage) had fixed/been lost after time dt
            dt_fix = np.sort(delta_ts[1]) * mo_to_gen
            dt_loss = np.sort(delta_ts[0]) * mo_to_gen

            ## Equivalent but slower
            #[dt_loss2, dt_fix2] = [np.sort(map(itemgetter('time'),
            #                              filter(lambda x: x['nu0_index'] == iji,
            #                                     fixlost[key][y]))) * mo_to_gen
            #                       for y in xrange(2)]

            #if (dt_fix != dt_fix2).all() or (dt_loss != dt_loss2).all():
            #    raise ValueError('Not the same!')

            ax.plot(dt_loss, 1.0-np.linspace(0,P0,len(dt_loss)),
                    lw=2, c=colors[iji], ls = '-',
                    label=r'$\nu_0 \in ['+str(nu0s[0])+','+str(nu0s[1])+']$')    
            #ax.plot(dt_fix, np.linspace(0, P1, len(dt_fix)),
            #        lw=2, c=colors[iji], ls = '-')    
            ax.plot([0,2000], np.mean(nu0s)*np.ones(2), lw=2, ls='--', c=colors[iji])

        ax.set_xlabel(r'time interval $\Delta t$ [days]', fontsize=20)
        ax.set_ylabel(r'fraction of surviving mutations', fontsize=20)
        ax.set_title(key.capitalize()+'onymous')
        if i == 0:
            ax.legend(loc=1)
        ax.set_xlim(0, 2000)

    # Replot fraction of surviving mutations over time normalized together
    fig, ax = ppl.subplots(1, 1, figsize=(9, 9))
    lss = {'syn': '--', 'nonsyn': '-'}
    dt_losses = {'syn': [], 'nonsyn': []}
    for iji, nu0s in enumerate(nu0ss):
        for i, key in enumerate(['syn', 'nonsyn']):
 
            # Define the counts (lost, fixed, neither)
            counts = counts_all[key][iji]
            delta_ts = times_all[key][iji]
            P0 = 1.0 * counts[0] / sum(counts)
            P1 = 1.0 * counts[1] / sum(counts)
        
            # Plot the fixation/loss probability together as a function of time,
            # by simply counting how many alleles (in percentage) had fixed/been lost after time dt
            dt_fix = np.sort(delta_ts[1]) * mo_to_gen
            dt_loss = np.sort(delta_ts[0]) * mo_to_gen

            dt_losses[key].append(dt_loss)
            
            ax.plot(dt_loss, 1.0-np.linspace(0, 1,len(dt_loss)),
                    lw=2, c=colors[iji], ls=lss[key],
                    label=key+r' $\nu_0 \in ['+str(nu0s[0])+','+str(nu0s[1])+']$')

    ax.set_xlabel(r'time interval $\Delta t$ [days]', fontsize=20)
    ax.set_ylabel(r'fraction of surviving alleles (normalized)', fontsize=20)
    ax.legend(loc=1)
    ax.set_xlim(0, 2000)


    ppl.ion()
    ppl.show()
