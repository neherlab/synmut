#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       28/07/13
content:    Perform the simulations for the depression in fixation probability
            and analyze them at once, to avoid occupying lots of disk space.
            This script uses a short genome model (1000 sites).

'''
# Standard modules
import os
import sys
from operator import *
import numpy as np
from itertools import combinations
import time

# Env vars
JOBSSCRIPT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSSCRIPT)+'/'
JOBMODULES = JOBDIR+'modules'    
# Set the modules path for the cluster nodes (as opposed to my laptop)
on_laptop = any(map(lambda x: '/home/fabio' in x, os.sys.path))

#if (not on_laptop) and (JOBMODULES not in sys.path):
#    sys.path.insert(0, JOBMODULES)

# Custom modules (import local FFPopSim)
import FFPopSim as h



# Globals
VERBOSE = 2
L = 1000
crossover = 1e-3
epitope_length = 9

# Simulation length
dt = 10
n_time_points = 600
dt_sample = 20

# Analysis parameters
t_max = 3000
IS_ZERO = 0.000001
max_counts = 10000
nirands = 10
sfs_bins = np.logspace(-3, -0.1, 30)
savefig = True
sep = '\t'



# Classes
class haploid_highd_gradual(h.haploid_highd):
    '''Population with dynamic fitness landscape'''
    

    def __init__(self, *args, **kwargs):
        if 'nescepi' not in kwargs:
            raise ValueError('Please add the "nescepi" argument.')
        newkwargs = {key: value for key, value in kwargs.iteritems() if key != 'nescepi'}
        super(haploid_highd_gradual, self).__init__(*args, **newkwargs)
        self.nescepi = kwargs['nescepi']
        self.epitopes = []
        self.epitope_times = []
        self.epitope_fitness = []


    def add_epitope(self, value):

        # Complex epitope creation
        if self.nescepi > 1:
            # Choose a random good position
            pos = np.random.randint(self.L - epitope_length)
            epitope = [p for p in range(pos, pos + epitope_length) if (p + 1) % 3]
            np.random.shuffle(epitope)
            epitope = epitope[:self.nescepi]

            add, epis = create_epitope(epitope, value)
            fitness = self.get_trait_additive()
            for i, v in add.iteritems():
                fitness[i] = v
            self.set_trait_additive(fitness)
            for loci, v in epis.iteritems():
                self.add_trait_coefficient(v, loci)

            # Store changes
            self.epitopes.append(epitope)
            self.epitope_times.append(self.generation)
            self.epitope_fitness.append(value)

        # Single-locus epitope (simpler)
        else:
            # Choose a random good position
            pos = 2
            while (((pos + 1) % 3) == 0) or (pos in self.epitopes):
                pos = np.random.randint(self.L)
    
            # Change fitness coefficient at that site
            fitness = self.get_trait_additive()
            fitness[pos] = value
            self.set_trait_additive(fitness)
    
            # Store changes
            self.epitopes.append([pos])
            self.epitope_times.append(self.generation)
            self.epitope_fitness.append(value)




# Functions
def create_epitope(epitope, adaptive_effect):
    '''Assign fitness coefficients to the epitope
    
    We want the single mutant to have a fitness advantage over the wildtype
    that is randomly distributed around adaptive_effect (times two). Hence,
    we must fiddle a bit around with the coefficients.

    Parameters:
        - epitope: sites that are part of the epitope and get beneficial
                   effects if flipped. Synonymous sites should NOT be part of
                   this list;
        - adaptive_effect: the magnitude of a single mutation in the epitope,
                   i.e. the difference between wildtype and single mutant.

    Returns:
        - add: dictionary of additive effects
        - epis: dictionary of epistatic effects

    See also:
    - /home/fabio/university/phd/notes/epitopes/epitopes.html

    '''
    L = len(epitope)
    # This makes sure that that locus is actually more beneficial than one
    # hitchhiker (to save simulation time)
    a = adaptive_effect / (1 << (L - 1))
    add = {i: a for i in epitope}
    epis = {}
    for order in xrange(2, L + 1):
        if order % 2:
            sign = +1
        else:
            sign = -1
        for c in combinations(epitope, order):
                epis[c] = sign * a
    return add, epis



# Script
if __name__ == '__main__':

    t0_sim = time.time()

    ###########################################################################
    # INPUT
    ###########################################################################
    # Take folder from input
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please specify a folder.')
    # Output dir
    results_dir = args[1].rstrip('/')+'/'

    # Read number of iterations to perform within this script
    if len(args) < 4:
        raise ValueError('Please specify how many times to iterate in one '+
                         'script.')
    simsn = int(args[2])

    # Read parameters from input
    if len(args) < 18:
        raise ValueError('Please specify population size, coinfection rate, \
                          deleterious effect, adaptive effect, \
                          deleterious fraction, adaptive density, \
                          escape mutations per epitope, \
                          generation time, mutation rate, \
                          deleterious effect of nonsyn non-escape mutations, \
                          time dependent selection coefficient, \
                          and bin edges for the analysis (in this order).')
    N = int(float(args[3]))
    coi = float(args[4])
    del_effect = float(args[5])
    adaptive_effect = float(args[6])
    del_fraction = float(args[7])
    adaptive_density = float(args[8])
    escape_muts_per_epitope = int(args[9])
    lag_time = 1.0 / adaptive_density
    generation_time = float(args[10])
    mu = float(args[11])
    nonsyn_dead_eff = float(args[12])
    time_dep_selection_coeff = float(args[13])
    nu0ss = np.array(args[14:], float).reshape(len(args[14:]) / 2, 2)

    # Note: The input parameters are in the folder name!

    #######################################################################
    # Prepare structures for results
    #######################################################################
    # Fixation/loss counts
    counts_syn = [{'fixed': 0, 'lost': 0, 'never': 0} for i in nu0ss]
    counts_nonsyn = [{'fixed': 0, 'lost': 0, 'never': 0} for i in nu0ss]

    # Rate of substitution
    rate_subs = {'syn': np.ma.masked_all(simsn),
                 'nonsyn': np.ma.masked_all(simsn)}

    # Counts between 25% and 75%
    poly_25to75s = {'syn': np.ma.masked_all(simsn),
                    'nonsyn': np.ma.masked_all(simsn)}

    # Fitness lanscape parameters
    del_effs = np.ma.masked_all(simsn)
    ada_effs = np.ma.masked_all(simsn)
    del_fracs = np.ma.masked_all(simsn)
    ada_dens = np.ma.masked_all(simsn)
    trashed = 0
    afs_all = []

    # Count how many simulations contribute to the counts for each bin
    contribs = simsn * np.ones(len(nu0ss), int)

    # Repeat everything simsn number of times
    for i in xrange(simsn):

        #######################################################################
        # PREPARE SIMULATION
        #######################################################################
        if VERBOSE >= 1:
            sys.stdout.write('Iteration n.'+str(i + 1)+' of '+str(simsn)+'...')
        if VERBOSE >= 2:
            sys.stdout.write('Preparing...')

        # Create population
        pop = haploid_highd_gradual(L, nescepi=escape_muts_per_epitope)
        pop.set_wildtype(N)
        pop.mutation_rate = mu * generation_time
        pop.outcrossing_rate = coi * generation_time
        pop.crossover_rate = crossover
    
        # Set replication (fitness) landscape
        is_third = (np.arange(L) % 3) == 2
        is_fstsnd = -is_third
        fitness = -del_effect * np.ones(L) * generation_time
    
        # Add neutral loci
        if del_fraction < (1 - 1.0 / L):
            is_neutral = np.random.random(L) > del_fraction
            fitness[is_neutral] = 0

        # Kill nonsynonymous sites: they are very deleterious
        fitness[-is_third] = - nonsyn_dead_eff * generation_time

        ### THE INITIAL LANDSCAPE IS READY HERE! ###
        pop.set_trait_additive(fitness, 0)
        if VERBOSE >= 2:
            sys.stdout.write('DONE. Evolving...')
    
        # Prepare output structures
        afs = []
        gs = []

        #######################################################################
        # EVOLVE
        #######################################################################
        pop.evolve(30)          # Burn-in
        gs.append(pop.generation)
        afs.append(pop.get_allele_frequencies())
        # Store data a la Shankarappa if on laptop (for supplementary figure)
        if on_laptop: seqs = [pop.random_genomes(10)]

        if escape_muts_per_epitope:
            # Do not add initial epitope
            #pop.add_epitope(np.random.exponential(adaptive_effect))
            last_epitope_time = pop.generation
            next_epitope_time = last_epitope_time + np.random.exponential(lag_time) * generation_time
    
        last_percent = -1
        for j in xrange(n_time_points - 1):
            pop.evolve(dt)
            g = pop.generation
    
            # Sample allele frequencies
            if (g - gs[-1]) >= dt_sample:
                gs.append(g)
                afs.append(pop.get_allele_frequencies())
                if on_laptop: seqs.append(pop.random_genomes(10))
    
            # Change the fitness landscape introducing new epitopes
            # Note: now this works by exponential times
            if escape_muts_per_epitope:
                if g >= next_epitope_time:
                    for ii in xrange(int((g - last_epitope_time) / lag_time)):
                        pop.add_epitope(np.random.exponential(adaptive_effect) * generation_time)
                    last_epitope_time = g
                    next_epitope_time = last_epitope_time + np.random.exponential(lag_time) * generation_time

                # If time-dependent selection is involved (only simple epitopes),
                # check whether or not it has been caught already
                if time_dep_selection_coeff > 0:
                    neutralize_ind = []
                    for ie in xrange(len(pop.epitopes)):
                        af = pop.get_allele_frequency(pop.epitopes[ie][0])
                        is_flip = np.random.rand() > (1 - time_dep_selection_coeff * dt * af)
                        if is_flip:
                            neutralize_ind.append(ie)

                    # If any epitope is to be neutralized, do so
                    if len(neutralize_ind):
                        fitness = pop.get_trait_additive()
                        for ie in neutralize_ind:
                            # They turn dead again
                            fitness[pop.epitopes[ie][0]] = - nonsyn_dead_eff * generation_time
                        pop.set_trait_additive(fitness)
    
            # Output percent
            if VERBOSE >= 3:
                new_percent = 100 * j / (n_time_points - 1)
                if last_percent != new_percent:
                    sys.stdout.write('\n'+str(i)+': '+str(new_percent)+'% completed')
                    last_percent = new_percent
    
        # Transform into ndarrays and ignore the first 1000 generations
        gs = np.array(gs)
        ind_relaxed = gs > 1000
        gs = gs[ind_relaxed]
        afs = np.array(afs)[ind_relaxed]
        if on_laptop: seqs = np.array(seqs)[ind_relaxed]
    
        if VERBOSE >= 3:
            sys.stdout.write('\n'+str(i)+': ')
        if VERBOSE >= 2:
            sys.stdout.write('DONE. Analyzing...')
    
        ########################################################################
        # ANALYSIS
        ########################################################################
        # Get fitness landscape
        fitness = pop.get_trait_additive()

        # Adaptive effects have to be read from fitness changes
        changes = np.asarray(pop.epitope_times)

        # Rate of epitope creation
        # There are a couple of corner cases here: if everything happens at one
        # time point, then there's little way of determining the rate of new epitopes
        if not len(changes):
            # if no sweeps at all, the density is zero
            ada_density = 0
        elif len(changes) == 1:
            # ex: in case of a single switch put it negative
            # Note: this changes 15/08/2013
            ada_density = 1.0 / n_time_points / dt_sample
        else:
            # in case of several switches, it depends whether all happen at the same time
            dt_epi = np.diff(changes).mean()
            if dt_epi > 0:
                ada_density = 1.0 / dt_epi
            else:
                # ex: in case of a single sweep with several switches put it negative
                # Note: this changes 15/08/2013
                ada_density = 1.0 * len(changes) / n_time_points / dt_sample

        # Some simulations have no sweeps at all, those have no adaptive effects
        if ((fitness > 0) & is_fstsnd).any():
            ada_eff = fitness[(fitness > 0) & is_fstsnd].mean()
        else:
            ada_eff = 0

        # Assign fitness coefficients
        del_fracs[i] = (fitness < 0)[is_third].mean()
        del_effs[i] = -fitness[(fitness < 0) & is_third].mean()
        ada_dens[i] = ada_density
        ada_effs[i] = ada_eff

        # sim parameters
        L = afs.shape[1]
        n = len(gs)

        # Count the number of substitutions
        af0 = afs[0]
        af1 = afs[-1]
        consensus = (af0 > 0.5)
        consensus1 = (af1 > 0.5)
        is_changed = consensus1 != consensus
        tt = gs[-1] - gs[0]
        rate_subs['syn'][i] = 1.0 * (is_changed & is_third).sum() / tt / is_third.sum()
        rate_subs['nonsyn'][i] = 1.0 * (is_changed & is_fstsnd).sum() / tt / is_fstsnd.sum()

        # Measure the diversity per time point at random time points
        aftmp = afs[np.random.randint(len(gs), size=nirands)]
        is_poly = (aftmp > 0.25) & (aftmp < 0.75)
        poly_25to75s['syn'][i] = 1.0 * (is_poly & is_third).sum() / nirands / is_third.sum()
        poly_25to75s['nonsyn'][i] = 1.0 * (is_poly & is_fstsnd).sum() / nirands / is_fstsnd.sum()

        # Establish last accepted t0, because simulations run of out steam after
        # a while (faster than real infections)
        il = (gs <= t_max).nonzero()[0].max()
   
        # Loop over initial frequencies
        for j, nu0s in enumerate(nu0ss):
            # Check whether we already have enough counts for this bin
            if ((sum(counts_syn[j].values()) >= max_counts) and
                (sum(counts_nonsyn[j].values()) >= max_counts)):
                contribs[j] -= 1
                continue

            # Alleles in this frequency bin
            is_nu0 = ((afs[:il+1] > nu0s[0]) & (afs[:il+1] <= nu0s[1])).any(axis=0)
   
            # Look at their fate
            candidates_syn = (is_nu0 & is_third).nonzero()[0]
            candidates_nonsyn = (is_nu0 & is_fstsnd).nonzero()[0]
            candidates_all = {'syn': candidates_syn,
                              'nonsyn': candidates_nonsyn}
            for key, candidates in candidates_all.iteritems():
                if key == 'syn':
                    counts = counts_syn[j]
                else:
                    counts = counts_nonsyn[j]
                
                # Store them
                for a in candidates:
                    traj = afs[:, a]
                    it0 = ((traj > nu0s[0]) & (traj <= nu0s[1])).nonzero()[0][0]
                    it1 = it0 + 1
                    while (it1 < n):
                        if traj[it1] > 0.95:
                            counts['fixed'] += 1
                            break
                        elif traj[it1] < 0.05:
                            counts['lost'] += 1
                            break
                        else:
                            it1 += 1
                    if it1 == n:
                        counts['never'] += 1

            if i < 30:
                afs_all.append(afs)
        
        if VERBOSE >= 1:
            sys.stdout.write('DONE.\n')


    #######################################################################
    # Calculate parameters as averages over all simulations
    #######################################################################
    if VERBOSE >= 1:
        sys.stdout.write('Averaging results over iterations...')

    # Fitness landscape
    del_eff = del_effs.mean()
    ada_eff = ada_effs.mean()
    del_frac = del_fracs.mean()
    ada_den = np.abs(ada_dens).mean()

    # Calculate probabilities as averages over all simulations
    Psyn = np.ma.masked_all(len(nu0ss))
    Pnonsyn = np.ma.masked_all(len(nu0ss))
    for (P, counts) in [(Psyn, counts_syn),
                        (Pnonsyn, counts_nonsyn)]:
        for j in xrange(len(nu0ss)):
            dic = counts[j]
            su = sum(dic.values())
            if su > 0:
                P[j] = 1.0 * dic['fixed'] / su

    # Calculate the average number of counts for each frequency bin
    counsn_syn = 1.0 / np.array([sum(dic.values()) for dic in counts_syn]) / contribs
    counsn_nonsyn = 1.0 / np.array([sum(dic.values()) for dic in counts_nonsyn]) / contribs

    # Calculate the substitution rates as an average over all simulations
    rates = {'syn': np.mean(rate_subs['syn']),
             'nonsyn': np.mean(rate_subs['nonsyn'])}

    # Calculate the synonymous diversity as an average over all simulations
    poly_25to75 = {'syn': np.mean(poly_25to75s['syn']),
                   'nonsyn': np.mean(poly_25to75s['nonsyn'])}

    # Calculate allele frequency spectra
    afs_all = np.array(afs)
    af_spectrum = np.histogram(afs.ravel(), bins=sfs_bins)

    if VERBOSE >= 1:
        sys.stdout.write('DONE. Saving...')

    ###########################################################################
    # OUTPUT
    ###########################################################################
    # Write parameters
    with open(results_dir+'parameters.dat', 'w') as f:
        f.write('# '+sep.join(['del_effect', 'adaptive_effect',
                                'del_fraction', 'adaptive_rate',
                                'N', 'coi', 'mu',
                                'trashed_simulations',
                                'good_simulations'])+'\n')
        f.write(sep.join(map(str, [del_eff, ada_eff,
                                   del_frac, ada_den,
                                   N, coi, mu,
                                   trashed, simsn - trashed]))+'\n')

    # Write fixation probabilities
    with open(results_dir+'Pfix.dat', 'w') as f:
        f.write('# '+sep.join(['nu0_start', 'nu0_end',
                               'Psyn', 'Pnonsyn',
                               'counts_syn', 'counts_nonsyn'])+'\n')
        for j, nu0s in enumerate(nu0ss):
            f.write(sep.join(map(str, [nu0s[0], nu0s[1],
                                       Psyn[j], Pnonsyn[j],
                                       counsn_syn[j], counsn_nonsyn[j]]))+'\n')

    # Write substitution rates
    with open(results_dir+'substitutions.dat', 'w') as f:
        f.write('# '+sep.join(['syn', 'nonsyn'])+'\n')
        f.write(sep.join(map(str, [rates['syn'], rates['nonsyn']]))+'\n')

    # Write polymorphic counts
    with open(results_dir+'poly_25to75.dat', 'w') as f:
        f.write('# '+sep.join(['syn', 'nonsyn'])+'\n')
        f.write(sep.join(map(str, [poly_25to75['syn'], poly_25to75['nonsyn']]))+'\n')
    
    # Write allele frequency spectrum
    with open(results_dir+'sfs.dat', 'w') as f:
        f.write('# bins\tnumber of alleles')
        for j in xrange(len(sfs_bins) - 1):
            f.write('\t'.join(map(str, [sfs_bins[j], af_spectrum[0][j]]))+'\n')
        f.write(str(sfs_bins[-1])+'\n')

    if VERBOSE >= 1:
        sys.stdout.write('DONE.\n')

    # Plot divergence and diversity if on laptop
    if on_laptop:
        divergence = seqs.mean(axis=1).mean(axis=1)
        diversity = []
        for i in xrange(len(gs)):
            s0 = seqs[i, np.random.randint(seqs.shape[1], size=seqs.shape[1])]
            s1 = seqs[i, np.random.randint(seqs.shape[1], size=seqs.shape[1])]
            diversity.append((s0 != s1).mean())
        diversity = np.asarray(diversity)

        from matplotlib import rcParams
        rcParams.update({'axes.labelsize': 17, 
                         'text.fontsize': 17,
                         'legend.fontsize': 17,
                         'xtick.labelsize': 17,
                         'ytick.labelsize': 17,
                         'text.usetex': True})
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(gs, divergence, c='b', lw=2, label='Divergence')
        plt.plot(gs, diversity, c='r', lw=2, label='Diversity')
        plt.xlabel('Time [generations]')
        plt.legend(loc=2)
        plt.tight_layout()
        if savefig:
            for ext in ['svg', 'pdf']:
                plt.savefig('/home/fabio/publications/synmut/figures/simulations_diversity_divergence_example.'+ext)

        plt.ion()
        plt.show()


    t1_sim = time.time()
    if VERBOSE >= 1:
        sys.stdout.write('The simulation took '+str(t1_sim - t0_sim)+' seconds.\n')

