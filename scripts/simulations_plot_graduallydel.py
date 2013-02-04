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
datafile = 'parameters_fixation_areas.dat'



# Script
if __name__ == '__main__':


    # Input results dir
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please specify a results folder')
    results_dir = args[1].rstrip('/')+'/'

    # Load data using record arrays
    tb = np.genfromtxt(results_dir+datafile, unpack=False,
                       dtype=zip(('del_eff', 'ada_eff',
                                  'del_frac', 'ada_rate',
                                  'trashed', 'good',
                                  'Asyn', 'Anonsyn',
                                  'count_syn', 'counts_nonsyn',
                                  'subs_syn', 'subs_nonsyn',
                                  'poly_syn', 'poly_nonsyn',
                                  'foldername'),
                                 (float, float,
                                  float, float,
                                  int, int,
                                  float, float,
                                  float, float,
                                  float, float,
                                  float, float,
                                  'S200')))

    # Filter out NaNs
    is_nonnan = -np.isnan(tb['Asyn'])
    tb = tb[is_nonnan]
    ind = np.argsort(tb['del_eff'])

    # Filter only some parameters
    ind = [ind[1], ind[2], ind[4]]

    # Limits of deleff
    colors = ['b', 'g', 'r', 'cyan', 'magenta']
    labels = [r'$10^{-3.5}$', r'$10^{-3}$', r'$10^{-2.5}$']

    # Prepare figure
    fig = plt.figure()

    # Get the data
    for ii, i in enumerate(ind):

        folder = tb['foldername'][i]
        fn = results_dir+folder
        nu0, nu1, Psyn, Pnonsyn, csyn, cnonsyn = np.genfromtxt(fn+'/Pfix.dat',
                                                               unpack=True,
                                                               usemask=True,
                                                               missing_values='--')
        deleff = tb['del_eff'][i]

        # Plot fixation probabilities
        x = [0] + list(0.5 * (nu0 + nu1))
        ysyn = [0] + list(Psyn)
        ynonsyn = [0] + list(Pnonsyn)
        plt.plot(x, ysyn,
                 color=colors[ii],
                 lw=2)
        # Plot continuation line to one
        plt.plot([x[-1], 1], [ysyn[-1], 1], color=colors[ii], ls='--')

    # Plot diagonal
    plt.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100), color='k',
             ls='--', lw=1.5)
    plt.xlabel(r'$\nu$', fontsize=18)
    plt.ylabel(r'$P_{\text{fix}}$', fontsize=18)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.text(0.5, 0.28, labels[0], fontsize=16, color=colors[0])
    plt.text(0.7, 0.18, labels[1], fontsize=16, color=colors[1])
    plt.text(0.9, 0.08, labels[2], fontsize=16, color=colors[2])

    # Plot syn diversity
    indall = tb['del_eff'].argsort()
    ax = fig.add_axes([0.23, 0.63, 0.25, 0.22])
    ax.plot(tb['del_eff'][indall], tb['poly_syn'][indall], color='k',
            lw=2)
    ax.set_xlabel('del eff')
    ax.set_ylabel('syn div')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Add circles
    import matplotlib.transforms as transforms
    from matplotlib.patches import Ellipse
    for ii, i in enumerate(ind):
        x = np.log10(tb['del_eff'][i] / 1e-4) / 2
        y = np.log10(tb['poly_syn'][i] / 1e-4) / 3
        circle = Ellipse((x, y),
                         0.05 * 2, 0.05 * 3,
                         transform=ax.transAxes,
                         linewidth=3,
                         facecolor='none',
                         edgecolor=colors[ii])
        ax.add_artist(circle)


    plt.ion()
    plt.show()
