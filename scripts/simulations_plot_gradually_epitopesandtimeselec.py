# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/01/13
content:    Plot the fixation probability of increasingly deleterious sims and 
            allele freqs trajectories of some of its simulations.
'''
# Modules
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.cm as cm



# Globals
results_dir = '../data/simulations/gradually_epitopesandtimeselec/'
datafile = 'parameters_fixation_areas.dat'



# Functions
def filterfolders(folders):
    '''Filter and sort folders'''
    tmp = folders[0].split('_')
    if ('nescepi' in tmp) or ('timeselec' in tmp):
        tmp = tmp[:-2]
    base = '_'.join(tmp)

    ind = {}
    for i, f in enumerate(folders):
        if 'nescepi_3' in f:
            ind[1] = i
        elif 'nescepi_6' in f:
            ind[2] = i
        elif 'timeselec_0.005' in f:
            ind[3] = i
        elif 'timeselec_0.01' in f:
            ind[4] = i
        elif ('nescepi' not in f) and ('timeselec' not in f):
            ind[0] = i
    return [folders[ind[i]] for i in xrange(5)]



# Script
if __name__ == '__main__':

    # Folders
    g = os.walk(results_dir)
    folders = np.array(g.next()[1])

    # Sort them smartly
    folders = filterfolders(folders)
    labels = ['naive', '3 mut/epi', '6 mut/epi', 'tds 0.005', 'tds 0.01']
    colors = ['b', 'g', 'r', 'cyan', 'magenta']

    # Prepare figure
    fig = plt.figure()

    # Get the data
    for i, folder in enumerate(folders):

        fn = results_dir+folder
        nu0, nu1, Psyn, Pnonsyn, csyn, cnonsyn = np.genfromtxt(fn+'/Pfix.dat',
                                                               unpack=True,
                                                               usemask=True,
                                                               missing_values='--')
        # Plot fixation probabilities
        x = [0] + list(0.5 * (nu0 + nu1)) + [1]
        ysyn = [0] + list(Psyn) + [1]
        ynonsyn = [0] + list(Pnonsyn) + [1]
        plt.plot(x, ysyn,
                 color=colors[i],
                 ls='--',
                 lw=2)
        plt.plot(x, ynonsyn,
                 label=labels[i],
                 color=colors[i],
                 ls='-',
                 lw=2)

    # Plot diagonal
    plt.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100), color='k',
             ls='--', lw=1.5)
    plt.xlabel(r'$\nu$', fontsize=18)
    plt.ylabel(r'$P_{\text{fix}}$', fontsize=18)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.55)
    plt.legend(loc=2)
    plt.title('Fixation probabilities in complex scenarios', fontsize=16)

    plt.ion()
    plt.show()
