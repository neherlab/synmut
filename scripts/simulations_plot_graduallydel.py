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

from matplotlib import rcParams
rcParams.update({'axes.labelsize': 20, 
                 'text.fontsize': 20,
                 'legend.fontsize': 20,
                 'xtick.labelsize': 20,
                 'ytick.labelsize': 20,
                 'text.usetex': True})

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.cm as cm



# Globals
results_dir = '../data/simulations/gradually_del/'
subdirs = ['low', 'medium', 'high']
del_effs = {'low': 0.003, 'medium': 0.001, 'high': 0.0003,
            'very_high': 0.0001, 'very_low': 0.01} # Multiply by 2!
colors = ['b', 'g', 'r']
symbols = ['o', 's', '^']


# Script
if __name__ == '__main__':

    fig = plt.figure()

    # Plot fixation probability
    for i, subdir in enumerate(subdirs):
        table = np.genfromtxt(results_dir+subdir+'/Pfix.dat',
                              usecols=[0, 1, 2], unpack=True,
                              usemask=True)
        nu0 = [0] + list(0.5 * (table[0] + table[1]))

        Pfix = table[2]
        Pfix[np.isnan(Pfix)] = 0
        Pfix = [0] + list(Pfix)
        color = colors[i]

        plt.plot(nu0, Pfix, lw=2, c=color)
        plt.text(1.03 * nu0[-1], 0.9 * Pfix[-1],
                 r'$s_d = $ '+str(2 * del_effs[subdir]),
                 fontsize=20, color=color)
        plt.scatter(nu0[2], Pfix[2], s=150,
                    marker=symbols[i], edgecolor=color,
                    facecolor='none', lw=2)

    plt.xlabel(r'$\nu$', fontsize=24)
    plt.ylabel(r'$P_{fix}$', fontsize=22)
    plt.plot([0, 1], [0, 1], lw=1.5, ls='--', c='k')
    plt.arrow(0.5, 0.4, 0, -0.22,
              facecolor='k', edgecolor='k',
              width=0.007, head_width=0.05, head_length=0.12,
              overhang=0.3)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)

    # Inset: plot synonymous diversity
    from matplotlib.patches import Rectangle
    plt.gca().add_patch(Rectangle((-0.03, 0.50), 0.52, 0.53,
                                  edgecolor='k', facecolor='none', lw=1.2))
    ax = fig.add_axes([0.25, 0.67, 0.25, 0.25])
    div = []
    for i, subdir in enumerate(subdirs):
        div.append(np.loadtxt(results_dir+subdir+'/poly_25to75.dat',
                              usecols=[0]))
    div = np.array(div)
    del_eff_x = [del_effs[l] for l in subdirs]
    ax.plot(del_eff_x, div, lw=2, c='k')
    ax.plot([3e-4, 2.5e-4], [9e-3, 1.5e-3], lw=1.2, c='purple')
    ax.text(1.5e-4, 2e-4, 'observed\ndiversity', fontsize=16)

    ax.set_xlabel(r'$s_d$', fontsize=18)
    ax.set_ylabel(r'$P_{interm}$', fontsize=18)
    ax.plot([1e-4, 1e-2], [0.01] * 2, color='purple', lw=1.5, ls='--')
    for i, subdir in enumerate(subdirs):
        ax.scatter(del_effs[subdir], div[i],
                   s=100, marker=symbols[i],
                   edgecolor=colors[i], facecolor='none', lw=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-4, 1e-2)


    panel = plt.figtext(0.02, 0.92, 'A', fontsize=24)

    plt.tight_layout(rect=(0, 0, 1, 1))


    if True:
        for ext in ['svg', 'pdf']:
            plt.savefig('../figures/simulations_graduallydel.'+ext)

    plt.ion()
