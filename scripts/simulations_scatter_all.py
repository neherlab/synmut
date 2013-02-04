# vim: fdm=marker
'''
author:     Fabio Zanini
date:       14/01/13
content:    Scatter plot and further visualization of the results of random
            gradual simulations.
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

# Set limits for data-like sims
Asynmin = -0.5
Asynmax = -0.15
poly_synmin = 0.0025
poly_synmax = 0.010

# (nonsyn gates, if only nonsyn data-like results are wished)
Anonsynmax = 0.5
Anonsynmin = -0.5
tmpsnsmin = 0e-6


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

    # Select only configurations with the right output (Asyn in the right
    # window)
    is_good = (tb['Asyn'] < Asynmax) & (tb['Asyn'] > Asynmin)

    # only with enough synonymous diversity
    is_good &= (tb['poly_syn'] > poly_synmin) & (tb['poly_syn'] < poly_synmax)

    # (nonsynonymous gating)
    # only with small Anonsyn
    is_good &= (tb['Anonsyn'] < Anonsynmax) & (tb['Anonsyn'] > Anonsynmin)
    # only in the right window for nonsyn substitutions
    is_good &= tb['subs_nonsyn'] > tmpsnsmin

    # Divide in classes
    tbg = tb[is_good]
    tbng = tb[-is_good]

    # Create figure
    fig, axs = plt.subplots(1, 1, figsize=[ 8 ,   6.0875])

    # Plot the simulations
    colors = cm.jet([int(255.0 * (i + 2.5)) for i in np.log10(tbg['ada_eff'])])
    sizes = 10 + 120 * (np.log10(tbg['ada_rate']) + 3)
    ind = sizes.argsort()[::-1]
    plt.scatter(np.log10(tbg['del_eff'])[ind], tbg['del_frac'][ind],
                color=colors[ind], s=sizes[ind])
    plt.xlabel(r'$\log_{10}( \text{deleterious effect} )$', fontsize=16)
    plt.ylabel(r'deleterious fraction', fontsize=16)
    plt.gca().xaxis.set_label_coords(0.5, -0.07)
    plt.gca().yaxis.set_label_coords(-0.1, 0.5)
    plt.xlim(-4.1, -1.9)
    plt.ylim(0.72, 1.02)
    plt.title(('color is adaptive effect (1e-2.5 to 1e-1.5)\n'+
               'size is rate of new epitopes (1e-3 to 1e-2)'), fontsize=14)
    
    # Draw a rectangle around the parameter priors
    rec = Rectangle(xy=(-4, 0.75),
                    width=2,
                    height=0.25,
                    facecolor='none', edgecolor='red',
                    lw=1.5)
    plt.gca().add_patch(rec)
    txt2 = plt.text(-2.25, 0.735, 'flat prior', color='red')

    plt.ion()
    plt.show()
