#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/08/12
content:    Study fixation and loss of synonymous/nonsynonymous changes by a
            short-genome model (1000 sites). This script is the one eventually
            used to produce the simulations for the paper.
'''
# Standard modules
import os
import sys
from operator import *
import numpy as np
from itertools import combinations

# Env vars
JOBSSCRIPT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSSCRIPT)+'/'
JOBMODULES = JOBDIR+'modules'    
if JOBMODULES not in sys.path: sys.path.insert(0,JOBMODULES)

# Custom modules (import local FFPopSim)
import FFPopSim as h



# Globals
VERBOSE = 0
L = 1000
mu = 2e-5
crossover = 1e-3
epitope_length = 9

# Simulation length
dt = 10
n_time_points = 600
dt_sample = 20



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

    ###########################################################################
    # INPUT
    ###########################################################################
    # Take folder from input
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please specify a folder.')

    # Read number of iterations to perform within this script
    if len(args) < 4:
        raise ValueError('Please specify how many times to iterate in one '+
                         'script.')
    quot = int(args[2])

    # Read parameters from input
    if len(args) < 10:
        raise ValueError('Please specify population size, coinfection rate, \
                          deleterious effect, adaptive effect, \
                          deleterious fraction, adaptive density, \
                          and escape mutations per epitope (in this order).')
    N = int(float(args[3]))
    coi = float(args[4])
    del_effect = float(args[5])
    adaptive_effect = float(args[6])
    del_fraction = float(args[7])
    adaptive_density = float(args[8])
    escape_muts_per_epitope = int(args[9])
    lag_time = 1.0 / adaptive_density

    # Check for time-dependent selection (only in some sims)
    if len(args) > 10:
        time_dep_selection_coeff = float(args[10])
    else:
        time_dep_selection_coeff = 0


    # Repeat everything quot number of times (this is only an efficiency issue)
    for iquot in xrange(quot):

        #######################################################################
        # PREPARE
        #######################################################################
        sys.stdout.write('Preparing...\n')
        sys.stdout.write('Iteration n.'+str(iquot + 1)+'...\n')
        sys.stdout.flush()

        SGEnum = (int(os.getenv('SGE_TASK_ID')) - 1) * quot + iquot + 1
        SGEstr = str(SGEnum)

        # Set folder name
        if os.getenv('SGE_TASK_ID'):
            results_dir = args[1].rstrip('/')+'/'+SGEstr+'/'
        else:
            results_dir = args[1].rstrip('/')+'/test/'
    
        # Check existence of the folder and create it
        if os.path.isdir(results_dir):
            sys.stderr.write('Folder present! Skipping.\n')
        else:
            os.mkdir(results_dir)
            if VERBOSE: sys.stdout.write('Folder created: '+results_dir+'\n')
    
        # Create population
        pop = haploid_highd_gradual(L, nescepi=escape_muts_per_epitope)
        pop.set_wildtype(N)
        pop.mutation_rate = mu
        pop.outcrossing_rate = coi
        pop.crossover_rate = crossover
    
        # Set replication (fitness) landscape
        is_third = (np.arange(L) % 3) == 2
        fitness = -del_effect * np.ones(L)
    
        sys.stdout.write('Adding neutral loci...')
        sys.stdout.flush()
    
        # Add neutral loci
        if del_fraction < (1 - 1.0 / L):
            is_neutral = np.random.random(L) > del_fraction
            fitness[is_neutral] = 0

        # Kill nonsynonymous sites: they are very deleterious
        fitness[-is_third] = -0.01

        ### THE LANDSCAPE IS READY HERE! ###

        # Set landscape
        pop.set_trait_additive(fitness, 0)
        
        sys.stdout.write('DONE\nEvolving...\n')
        sys.stdout.flush()
    
        # Prepare output structures
        afs = []
        gs = []
    
        #######################################################################
        # EVOLVE
        #######################################################################
        pop.evolve(30)          # Burn-in
        gs.append(pop.generation)
        afs.append(pop.get_allele_frequencies())
        pop.add_epitope(np.random.exponential(adaptive_effect))
        last_epitope_time = 0
    
        for i in xrange(n_time_points - 1):
            pop.evolve(dt)
            g = pop.generation
    
            # Sample allele frequencies
            if (g - gs[-1]) >= dt_sample:
                gs.append(g)
                afs.append(pop.get_allele_frequencies())
    
            # Change the fitness landscape introducing new epitopes
            if (g - last_epitope_time) >= lag_time:
                for i in xrange(int((g - last_epitope_time) / lag_time)):
                    pop.add_epitope(np.random.exponential(adaptive_effect))
                    last_epitope_time = g

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
                        fitness[pop.epitopes[ie][0]] = -0.01
                    pop.set_trait_additive(fitness)
    
            sys.stdout.write(str(100 * i / (n_time_points - 1))+'% completed\n')
            sys.stdout.flush()
    
        # Transform into ndarrays
        gs = np.array(gs)
        afs = np.array(afs)
    
        sys.stdout.write('100% completed\n')
        sys.stdout.write('Saving...')
        sys.stdout.flush()
    
        #######################################################################
        # OUTPUT
        #######################################################################
        # Save allele frequencies and generations
        np.save(results_dir+'generations.npy', gs)
        np.save(results_dir+'allele_frequencies.npy', afs)
    
        # Save initial fitness landscape
        np.savetxt(results_dir+'fitness_landscape.dat', fitness)
    
        # Save fitness changes
        n = len(pop.epitopes)
        with open(results_dir+'fitness_changes.dat', 'w') as f:
            f.write('# '+'\t'.join(['generation', 'site', 'fitness'])+'\n')
            for i in xrange(n):
                f.write('\t'.join([str(pop.epitope_times[i]),
                                   ','.join(map(str,pop.epitopes[i])),
                                   str(pop.epitope_fitness[i])])+'\n')
    
        sys.stdout.write('DONE\n')
