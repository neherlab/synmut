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
import Bio

# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
import modules.parser_Bunnik2008 as pB
from modules.helper import translate
from modules.alphabet import alpha



# Globals
# Verbosity level (0-3)
VERBOSE = 1

# Initial frequency bins
nu0ss = [[0.05, 0.15],
         [0.15, 0.25],
         [0.25, 0.4],
         [0.4, 0.6],
         [0.6, 0.8]]

# Bad patients (selected manually in a conservative way looking e.g. at PCA)
bad_patients = {'Bunnik': ['ACH19542', 'ACH19768'],
                'Shankarappa': ['p4', 'p7', 'p8', 'p11']}

# Classes of mutations to be shown
classes = ['syn_Shanka',
           'syn_nonShanka',
           'nonsyn']

# Classes of mutations including script intermediates
classes_all = ['nonsyn',
               'syn',
               'syn_Shanka',
               'syn_nonShanka']

# Colors for the plot
colors = {'syn_Shanka': 'red',
          'syn_nonShanka': 'green',
          'nonsyn': 'blue'}

# Labels for the plot
labels = {'syn_Shanka': 'synonymous, C2-V5',
          'syn_nonShanka': 'synonymous, other',
          'nonsyn': 'nonsynonymous'}

# Distance from last time point to be considered likely to reach either boundary
# before the sequencing stops (later sequences will probably stay floating).
# Note: this is not used in the plot (only for testing).
maxt = {'Shankarappa': 40,  # Months
        'Bunnik': 1200}     # Days



# Script
if __name__ == '__main__':

    # Define the patients (excluding problematic ones)
    patientsB = pB.parse_sequences(reference='SHAPE', exclude=bad_patients['Bunnik'])
    patientsS = pS.parse_sequences(reference='SHAPE', exclude=bad_patients['Shankarappa'])
    patients = patientsB + patientsS

    # Counts (overall)
    # The first list of each nu0s is lost, the second fixed, the third floating
    # Moreover, record the patient number they came from
    counts = {x: np.zeros((len(nu0ss), 3, len(patients)), int) for x in classes}

    # Aggregate information from all patients
    for k, p in enumerate(patients):

        if VERBOSE >= 1:
            stderr.write(str(p)+'\n')

        # Filter the time points to sequenced times
        p.filter_for_longitudinal_genetics()
        n = len(p.visit)

        # Filter away gaps (conservative), keep reading frame
        is_nongap = ((np.array(p.seqs) == '-').sum(axis=0) == 0)
        for i in xrange(len(is_nongap)):
            if not is_nongap[i]:
                is_nongap[i - i % 3: i - i % 3 + 3] = False
    
        # Collect allele frequencies
        p.afs = afs = p.allele_frequency_trajectories().swapaxes(1,2).swapaxes(0,1)

        # Consensus at the initial time point to decide on mut and synonymity
        consensus = alpha[afs[0].argmax(axis=0)]

        # For Bunnik patients, keep track of the Shankarappa region and the
        # second exon of rev (including tat)
        if p.name[:3] == 'ACH':
            is_Shanka = p.only_Shanka()
            is_rev2 = p.only_rev2()

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

                if VERBOSE >= 3:
                    stderr.write('v'+str(p.visit[i])+'\n')

                # Time limit for floating/settling
                if p.name[:3] == 'ACH':
                    in_time = (p.days_from_seroconversion[-1] - p.days_from_seroconversion[i] > maxt['Bunnik'])
                else:
                    in_time = (p.mo_to_SC[-1]-p.mo_to_SC[i] > maxt['Shankarappa'])
    
                # Get allele frequencies
                af0 = afs[i]
    
                # Keep only af in the nu0 range and have never gaps
                ind_nu0 = (af0 > nu0s[0]) & (af0 < nu0s[1]) & is_nongap
                alleles = np.transpose(ind_nu0.nonzero())
    
                # Keep only actual mutations from consensus
                if len(alleles):
                    is_mut = np.array([alpha[al[0]] != consensus[al[1]] for al in alleles], bool)
                    alleles = alleles[is_mut]
    
                # Divide alleles by class
                all_cla = {x: [] for x in classes_all}
                if len(alleles):
                    is_syn = np.zeros(len(alleles), bool)
                    for j, al in enumerate(alleles):
                        pos = al[1]
                        mut = alpha[al[0]]
                        codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                        cod = codcons.copy()
                        cod[pos % 3] = mut
                        aacons = translate(codcons)
                        aa = translate(cod)
                        is_syn[j] = (aacons == aa)

                    all_cla['nonsyn'] = alleles[-is_syn]
                    all_cla['syn'] = alleles[is_syn]

                # Test more stringently for synonymity
                # (we want to avoid double-hits in one single codon)
                if len(all_cla['syn']):
                    seqs = np.array(p.seqs_from_visit(p.visit[i]))
                    is_single = np.zeros(len(all_cla['syn']), bool)
                    for j, al in enumerate(all_cla['syn']):
                        pos = al[1]
                        mut = alpha[al[0]]
                        # Check whether sequences have double-hits
                        codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                        seqal = seqs[seqs[:,pos] == mut]
                        codseq = seqal[:, pos - pos % 3: pos - pos % 3 + 3]
                        # If at least one sequence has a double hit, ignore the
                        # allele (very conservative)
                        is_single[j] = ((codcons != codseq).sum(axis=1) == 1).all()
                    all_cla['syn'] = all_cla['syn'][is_single]

                    # Divide by Shanka and outside
                    if p.name[:3] != 'ACH':
                        all_cla['syn_Shanka'] = all_cla['syn']
                        all_cla['syn_nonShanka'] = []
                    else:
                        all_cla['syn_Shanka'] = [a for a in all_cla['syn'] if is_Shanka[a[1]]]
                        all_cla['syn_nonShanka'] = [a for a in all_cla['syn'] if
                                ((not is_Shanka[a[1]]) and (not is_rev2[a[1]]))]
            
                # Check whether they fix/are lost in at least one time point
                # (the possible fates are:
                # - fixed
                # - lost
                # - floating since long
                # - floating because of lack of time
                # Note: the latter two fates are not analyzed further, only for
                # testing purposes
                for cla in classes:
                    alleles = all_cla[cla]
                    for al in alleles:
                        alf_sel0 = af0[al[0], al[1]]
    
                        # Check all later time points for fixation/loss
                        for l in xrange(i+1, n):
                            alf_sel1 = afs[l][al[0], al[1]]
            
                            # Increment the counts of lost/fixed
                            for j in xrange(2):
                                # Test whether the locus has been used already
                                if (alf_sel1 == j) and (is_nonused[al[1]]):
                                    counts[cla][iji, j, k] += 1
                                    is_nonused[al[1]] = False 
    
                        # If the allele has neither fixed not vanished by the end,
                        # and if there is enough time to consider it 'floating',
                        # count it in the third bucket
                        if is_nonused[al[1]] and in_time:
                            counts[cla][iji, 2, k] += 1
                            is_nonused[al[1]] = False
    
    if VERBOSE >= 1:
        stderr.write('\n')

    # Calculate average fixation probabilities starting from a certain nu0
    Pf = {x: np.zeros((len(nu0ss), 2)) for x in classes}
    for cla in classes:
        for i in xrange(counts[cla].shape[0]):
            tot = counts[cla][i].sum()
            Pf[cla][i, 1] = 0.5 * counts[cla][i, 2].sum() / tot
            Pf[cla][i, 0] = 1.0 * counts[cla][i, 1].sum() / tot + Pf[cla][i, 1]

    # Bootstrap to estimate errorbars
    nboot = 100
    Pfboot = {x: np.ma.masked_all((len(nu0ss), 2, nboot)) for x in classes}
    for cla in classes:
        for i in xrange(counts[cla].shape[0]):
            for k in xrange(nboot):
                if cla == 'syn_nonShanka':
                    is_good = np.random.randint(len(patientsB), size=len(patientsB))
                else:
                    is_good = np.random.randint(len(patients), size=len(patients))
                tmp = counts[cla][i,:, is_good]
                tot = tmp.sum()
                if tot == 0:
                    print cla
                    print nu0ss[i]
                    print is_good
                    print tmp
                    continue
                Pfboot[cla][i, 1, k] = 0.5 * tmp[2].sum() / tot
                Pfboot[cla][i, 0, k] = 1.0 * tmp[1].sum() / tot + Pfboot[cla][i, 1, k]

    # Calculate stddev of bootstrapping
    Pfstd = {x: np.zeros((len(nu0ss), 2)) for x in classes}
    for cla in classes:
        for i in xrange(counts[cla].shape[0]):
            Pfstd[cla][i] = np.std(Pfboot[cla][i], axis=1)

    # Calculate areas
    # Note: skip the last point (via i_final < len(x) - 1) because many simulations
    # have nothing there (not enough time or very noisy)
    nu0mid = [0.5 * sum(nu0s) for nu0s in nu0ss]
    x = [0] + list(nu0mid) + [1]
    ysyn = [0] + list(Pf['syn_Shanka'][:, 0]) + [1]
    i_final_offset = 2              # Skip last bin
    i_final = len(x) - 1 - i_final_offset
    Atot = 0.5 * x[i_final]**2
    Asyn = - Atot + sum([0.5 * (ysyn[i+1] + ysyn[i]) * (x[i+1] - x[i]) for i in xrange(i_final)])

    # Prepare figure and plot results
    nu0err = [0.5 * (nu0s[1] - nu0s[0]) for nu0s in nu0ss]
    fig, ax = ppl.subplots(1,1)
    for cla in ['nonsyn', 'syn_nonShanka', 'syn_Shanka']:
        ax.errorbar(nu0mid, Pf[cla][:, 0],
                    xerr=nu0err,
                    yerr=Pfstd[cla][:, 0],
                    lw=2,
                    label=labels[cla],
                    c=colors[cla])
    ax.plot([0, 1], [0, 1], lw=2, ls='--', c='k')
    ax.set_xlabel(r'$\nu$', fontsize=16)
    ax.set_ylabel(r'$P_\text{fix}(\nu)$', fontsize=16)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc=2)

    ppl.ion()
    ppl.show()
