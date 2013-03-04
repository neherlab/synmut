# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/07/12
content:    Correlate the fixation/loss of an allele with the amount of
            secondary structure and codon usage bias. Use three categories:
                1. polymorphic above 15% but finally extinct
                2. polymorphic above 15% and fixed
                3. never observed.
            All three refer to synonymous, derived alleles.

            Of course, the 3rd category includes alleles that are deleterious
            and other that were simply not hit by mutation + hitchhiking.
'''
# Standard modules
from sys import stderr
from operator import *
import numpy as np
import matplotlib.pyplot as ppl
from scipy import stats

# Custom modules
import sys
sys.path.insert(0, '.')
import modules.parser_Shankarappa as pS
import modules.parser_Bunnik2008 as pB
import modules.parser_LANL as pL
import modules.secondary as sec
from modules.helper import translate
from modules.alphabet import alpha


# Globals
# Verbosity level (0-3)
VERBOSE = 1

# Bad patients (selected manually in a conservative way looking e.g. at PCA)
bad_patients = {'Bunnik': ['ACH19542', 'ACH19768'],
                'Shankarappa': ['p4', 'p7', 'p8', 'p11']}

# Liu LANL data set code
Liu_code = 2432

# Thresholding method to spot polymorphisms ('counts' or 'frequency' are
# possible: only non-singletons are filtered by the first, only derived alleles
# in a certain frequency window by the second -- same results).
tmethod = 'counts'
# Details of the thresholding methods
thresholdcount = 2
thresholds = [0.15, 0.8]

# Level of conservativeness in defining 'synonymous' (one of 'normal', 'most',
# 'least'). Too strict means not many data, too loose means mixing up classes.
conservative_level = 'normal'

# Region of the genome to analyze
regions = ['V1', 'V2', 'V3', 'V4', 'V5']
only_regions = 1    # 1: only V regions and flanks, 0: all, -1: only nonV
flank = 100     # length in nucleotides of the flanks around V loops

# Keep only alleles that are really found in the SHAPE array (i.e. mutation from
# the NL4-3-like allele to a different allele).
only_from_SHAPE = True



# Script
if __name__ == '__main__':

    # Fixed/lost alleles
    fixed = []
    lost = []
    never = []
    fixlost = [lost, fixed, never]

    # Define the patients (excluding problematic ones)
    patientsB = pB.parse_sequences(reference='SHAPE', exclude=bad_patients['Bunnik'])
    patientsS = pS.parse_sequences(reference='SHAPE', exclude=bad_patients['Shankarappa'])
    patientsL = pL.parse_sequences(Liu_code, reference='SHAPE')
    patients = patientsS + patientsB + patientsL

    # Filter time points on which we have sequences
    for p in patientsS:
        p.filter_only_sequenced()
    for p in patientsB:
        p.filter_seqs_with_attribute('days_from_seroconversion')
    for p in patientsL:
        p.filter_seqs_with_attribute('days_from_seroconversion')
    
    # Aggregate information from all patients
    for k, p in enumerate(patients):

        if VERBOSE >= 1:
            stderr.write(str(p)+'\n')
    
        n = len(p.visit)
    
        # Filter away gaps (conservative), keep reading frame
        is_nongap = ((np.array(p.seqs) == '-').sum(axis=0) == 0)
        for i in xrange(len(is_nongap)):
            if not is_nongap[i]:
                is_nongap[i - i % 3: i - i % 3 + 3] = False

        # Restrict to V regions and surroundings. Fortunately, each
        # region is in only one chopped patient (approximately,
        # regardless for flanking regions)
        is_region = np.zeros(len(is_nongap), bool)
        start = end = None
        for region in regions:
            try:
                (start, end) = p.only_region(region)
                # Extend to flanking regions
                start = max(start - flank, 0)
                end = min(end + flank, len(is_nongap))
                is_region[start:end] = True
            except ValueError:
                pass
        if only_regions == 1:
            is_nongap &= is_region
        elif only_regions == -1:
            is_nongap &= -is_region
    
        # Collect allele frequencies and counts
        afs = p.allele_frequency_trajectories().swapaxes(1,2).swapaxes(0,1)
        acs = p.allele_count_trajectories().swapaxes(1,2).swapaxes(0,1)

        # Consensus at the initial time point to decide on
        # mutation and synonymity
        consind = afs[0].argmax(axis=0)
        consensus = alpha[consind]

        # Check what derived alleles are considered real polymorphisms, then get
        # fixed/extinct, then check for synonymity, then add them to the list
        if tmethod == 'counts':
            ind = (acs.sum(axis=0) >= thresholdcount)
        else:
            ind = ((afs > thresholds[0]) & (afs <= thresholds[1])).any(axis=0)
        ind &= is_nongap
        alleles = np.transpose(ind.nonzero())

        # Keep only actual mutations from consensus
        if len(alleles):
            is_mut = np.array([alpha[al[0]] != consensus[al[1]] for al in alleles], bool)
            alleles = alleles[is_mut]

        # Keep only synonymous changes
        if len(alleles):
            tmp = []
            for j, al in enumerate(alleles):
                pos = al[1]
                mut = alpha[al[0]]
                codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                cod = codcons.copy()
                cod[pos % 3] = mut
                aacons = translate(codcons)
                aa = translate(cod)
                if aacons == aa:
                    tmp.append(al)
                    if (VERBOSE >= 3):
                        print ''.join(codcons), '-->', ''.join(cod)
            alleles = tmp

        # Test more stringently for synonymity
        # (we want to avoid double-hits in one single codon)
        if conservative_level in ['most', 'normal']:
            if len(alleles):
                is_single = np.zeros(len(alleles), bool)
                for j, al in enumerate(alleles):
                    # Find the first polymorphic time point
                    if tmethod == 'counts':
                        tmpac = acs[:, al[0], al[1]]
                        i = (tmpac >= 1).nonzero()[0][0]
                    else:
                        tmpaf = afs[:, al[0], al[1]]
                        i = ((tmpaf > thresholds[0]) & (tmpaf <= thresholds[1])).nonzero()[0][0]
                    seqs = np.array(p.seqs_from_visit(p.visit[i]))
                    pos = al[1]
                    mut = alpha[al[0]]
                    # Check whether sequences have double-hits
                    codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                    seqal = seqs[seqs[:,pos] == mut]
                    codseq = seqal[:, pos - pos % 3: pos - pos % 3 + 3]
                    # If at least one sequence has a double hit, ignore the
                    # allele (very conservative)
                    is_single[j] = ((codcons != codseq).sum(axis=1) == 1).all()
                alleles = np.array(alleles)[is_single]

        # Most conservative method for calling synonymity (many data discarded)
        # filter out any case where there were more than one change in the codon
        if conservative_level in ['most']:
            if len(alleles):
                tmp = []
                for j, al in enumerate(alleles):
                    # Find the first polymorphic time point
                    if tmethod == 'counts':
                        tmpac = acs[:, al[0], al[1]]
                        i = (tmpac >= 1).nonzero()[0][0]
                    else:
                        tmpaf = afs[:, al[0], al[1]]
                        i = ((tmpaf > thresholds[0]) & (tmpaf <= thresholds[1])).nonzero()[0][0]
                    pos = al[1]
                    for codpos in xrange(3):
                        if (codpos != (pos % 3)):
                            xtmp = pos - pos % 3 + codpos
                            tmpaf = afs[i:, consind[xtmp] , xtmp]
                            if not (tmpaf == 1).all():
                                break
                    else:
                        tmp.append(al)
                alleles = tmp
        
        # CLASSES 1 + 2: fixed or lost
        # Check whether they fix/are lost in at least one time point
        if len(alleles):
            for al in alleles:
                tmpaf = afs[:, al[0], al[1]]
                # Find the first polymorphic time point
                if tmethod == 'counts':
                    tmpac = acs[:, al[0], al[1]]
                    i = (tmpac >= 1).nonzero()[0][0]
                else:
                    i = ((tmpaf > thresholds[0]) & (tmpaf <= thresholds[1])).nonzero()[0][0]
    
                # Check all later time points for fixation/loss
                for l in xrange(i+1, n):
                    # Increment the counts of lost
                    if tmpaf[l] == 0:
                        lost.append({'patient': str(p),
                                     'pos': al[1],
                                     'pp': p.SHAPE_p[al[1]],
                                     'r': p.SHAPE_r[al[1]],
                                     'all': alpha[al[0]],
                                     'all_cons': consensus[al[1]],
                                     'all_SHAPE': p.reference_seq[al[1]],
                                     'fate': 'lost',
                                     'codon_cons': ''.join(consensus[al[1] - 2: al[1] + 1])
                                     })
                        break

                    # Increament the counts of fixed
                    elif tmpaf[l] == 1:
                        fixed.append({'patient': str(p),
                                      'pos': al[1],
                                      'pp': p.SHAPE_p[al[1]],
                                      'r': p.SHAPE_r[al[1]],
                                      'all': alpha[al[0]],
                                      'all_cons': consensus[al[1]],
                                      'all_SHAPE': p.reference_seq[al[1]],
                                      'fate': 'fixed',
                                      'codon_cons': ''.join(consensus[al[1] - 2: al[1] + 1])
                                      })
                        break

        if (VERBOSE >= 2) and (i == n-2):
            stderr.write('\n')

        # CLASS 3: never observed
        # Conserved sites: second-codon positions are always nonsynonymous
        is_conserved = (afs == 1).all(axis=0).any(axis=0)
        is_nonsecond = np.arange(afs.shape[2]) % 3 != 1
        sites = (is_conserved & is_nonsecond & is_nongap).nonzero()[0]
        for pos in sites:
            for j, mut in enumerate(alpha):
                # Exclude consensus alleles and gaps
                if (mut != '-') and (mut != consensus[pos]):
                    # Check for synonymity (easier for conserved sites)
                    codcons = consensus[pos - pos % 3: pos - pos % 3 + 3]
                    cod = codcons.copy()
                    cod[pos % 3] = mut
                    aacons = translate(codcons)
                    aa = translate(cod)
                    if aacons == aa:
                        never.append({'patient': str(p),
                                      'pos': pos,
                                      'pp': p.SHAPE_p[pos],
                                      'r': p.SHAPE_r[pos],
                                      'all': mut,
                                      'all_cons': consensus[pos],
                                      'all_SHAPE': p.reference_seq[pos],
                                      'fate': 'never',
                                      'codon_cons': ''.join(consensus[pos - 2: pos + 1])
                                      })

                        # Do not double count conserved sites (a alleles will be
                        # overrepresented, but this is a technical issue)!
                        break

    ##############################################
    ##          RNA SECONDARY STRUCTURE         ##
    ##############################################
    # Filter only alleles away from NL4-3 = consensus
    if only_from_SHAPE:
        fixed_fromSHAPE = [a for a in fixed if (a['all'] != a['all_SHAPE']) and (a['all_cons'] == a['all_SHAPE'])]
        lost_fromSHAPE = [a for a in lost if (a['all'] != a['all_SHAPE']) and (a['all_cons'] == a['all_SHAPE'])]
        never_fromSHAPE = [a for a in never if (a['all'] != a['all_SHAPE']) and (a['all_cons'] == a['all_SHAPE'])]
        fixlost = [lost_fromSHAPE, fixed_fromSHAPE, never_fromSHAPE]
    
    # Get reactivities
    rs = [np.array(map(itemgetter('r'), a)) for a in fixlost]
    rs = map(lambda x: x[-np.isnan(x)], rs)
    rs = map(np.sort, rs)

    # Plot cumulative histogram for SHAPE reactivities VS fixed/lost
    fig, ax = ppl.subplots(1,1)
    ax.plot(rs[0], np.linspace(0, 1, len(rs[0])), lw=2, color='b',
            label='lost, '+str(len(rs[0])))
    ax.plot(rs[1], np.linspace(0, 1, len(rs[1])), lw=2, color='r',
            label='fixed, '+str(len(rs[1])))
    ax.plot(rs[2], np.linspace(0, 1, len(rs[2])), lw=2, color='g',
            label='never, '+str(len(rs[2])))

    ax.legend(loc=4)
    ax.set_title('Cumulative reactivity (syn)')
    ax.set_xlabel('SHAPE Reactivity')
    ax.set_xlim(0, 1.8)
 
    # Level of significance via Kolmogorov-Smirnov
    Prks = stats.ks_2samp(rs[0], rs[1])
    ax.text(1.2, 0.3, 'KS P = '+str(np.round(Prks[1], 4)), fontsize=14)

    # Plot some guides to the eye
    ax.plot([0,1.8], [0.35] * 2, c='k', ls='--')
    ax.plot([0,1.8], [0.78] * 2, c='k', ls='--')
    ax.plot([0.65] * 2, [0,1], c='k', ls='--')

    # Plot histograms of reactivities
    fig, ax = ppl.subplots(1,1, figsize=[ 6.5875,  6.075 ])
    ax.hist([rs[0], rs[2], rs[1]], bins=[0,0.33,0.67,1.0], normed=True,
            color=('b', 'g', 'r'),
            label=('lost', 'never', 'fixed'))
    ax.set_ylabel('normed histogram', fontsize=16)
    ax.set_xlabel('SHAPE Reactivity', fontsize=16)
    ax.set_title('SHAPE reactivity (syn)')
    ax.set_ylim(0, 2)
    ax.legend(loc=2)

    ppl.ion()
    ppl.show()
