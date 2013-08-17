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
import re
import numpy as np

from matplotlib import rcParams
rcParams.update({'axes.labelsize': 20, 
                 'text.fontsize': 20,
                 'legend.fontsize': 20,
                 'xtick.labelsize': 20,
                 'ytick.labelsize': 20,
                 'text.usetex': True})

from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.pyplot as plt


# Globals
datafile = 'parameters_fixation_areas.dat'
panel_name = 'B'
fig_name = 'simulations_syn_latin_contour'
figm_name = 'simulations_marginals_stratified'

# Set limits for good sims
Asynmin = -0.11
Asynmax = -0.06
poly_synmin = 0.008
poly_synmax = 0.015
i_final_offset = 2

# Condensation
recondense = False
savefig = False
plot_type = 'contour'
plot_things = True



# Functions
def condense_data(results_dir):
    '''Condense the results from many simulations into a single file'''
    results_dir = results_dir.rstrip('/')+'/'
    table = []

    # List all directories
    g = os.walk(results_dir)
    folders = g.next()[1]

    # Collect data from all simulations
    n_empty = 0
    for folder in folders:
        parfile = results_dir+folder+'/parameters.dat'
        Pfixfile = results_dir+folder+'/Pfix.dat'
        subfile = results_dir+folder+'/substitutions.dat'
        polyfile = results_dir+folder+'/poly_25to75.dat'

        if not all(map(os.path.isfile, (parfile, Pfixfile, subfile, polyfile))):
            n_empty += 1
            continue

        # Load data on parameters
        try:
            pars = np.loadtxt(parfile)
        except ValueError:
            continue

        # Load data on Pfix
        (nu0a, nu0b,
         Psyn, Pnonsyn,
         csyn, cnonsyn) = np.genfromtxt(Pfixfile, comments='#', unpack=True,
                                        missing_values=['--'], usemask=True)

        # Calculate areas
        # Note: skip the last point (via i_final < len(x) - 1) because many simulations
        # have nothing there (not enough time or very noisy)
        x = [0] + list(0.5 * (nu0a + nu0b)) + [1]
        ysyn = [0] + list(Psyn) + [1]
        ynonsyn = [0] + list(Pnonsyn) + [1]
        i_final = len(x) - 1 - i_final_offset
        Atot = 0.5 * x[i_final]**2
        Asyn = - Atot + sum([0.5 * (ysyn[i+1] + ysyn[i]) * (x[i+1] - x[i]) for i in xrange(i_final)])
        Anonsyn = - Atot + sum([0.5 * (ynonsyn[i+1] + ynonsyn[i]) * (x[i+1] - x[i]) for i in xrange(i_final)])

        # Summarize counts
        csyns = csyn.sum()
        cnonsyns = cnonsyn.sum()

        # Load substitutions
        rates = np.loadtxt(subfile, comments='#')

        # Load poly_25to75
        poly_25to75 = np.loadtxt(polyfile, comments='#')

        # Add population size, mutation and coinfection rate, and input adaptive and deleterious stuff
        fbl = os.path.basename(folder).rstrip('/').split('_')
        del_eff_in = fbl[10]
        ada_eff_in = fbl[13]
        del_frac_in = fbl[16]
        ada_rate_in = fbl[-1]

        # Latest simulations store the info
        if len(pars) == 9:
            pars = list(pars)
            (N, coi, mu) = [pars.pop(4) for ii in xrange(3)]
        # earlier, we used folder name (which is also fine)
        else:
            N = fbl[fbl.index('N') + 1]
            mu = fbl[fbl.index('mu') + 1]
            coi = fbl[fbl.index('coi') + 1]



        line = (list(pars)+
                [Asyn, Anonsyn]+
                [csyns, cnonsyns]+
                list(rates)+
                list(poly_25to75)+
                [N, mu, coi]+
                [del_eff_in, ada_eff_in, del_frac_in, ada_rate_in]+
                [folder])
        table.append(line)

    # Write results to file
    with open(results_dir+datafile, 'w') as f:
        f.write('# '+'\t'.join(['del_effect', 'ada_effect',
                                'del_frac', 'ada_frac',
                                'trashed', 'good',
                                'Asyn', 'Anonsyn',
                                'counts_syn', 'counts_nonsyn',
                                'subs_syn', 'subs_nonsyn',
                                'poly_25to75_syn', 'poly_25to75_nonsyn',
                                'N', 'mu', 'coi',
                                'del_effect_in', 'ada_eff_in',
                                'del_frac_in', 'ada_rate_in',
                                'folder'])+'\n')
        for line in table:
            f.write('\t'.join(map(str,line))+'\n')

    if n_empty:
        print str(n_empty)+' folders empty' 


def plot_contour(tb_raw, tbs, fields1, fields2, nbins=8):
    '''Plot contour plots of the two lists of observables in dataset tbs'''
    if isinstance(fields1, basestring): fields1 = [fields1]
    if isinstance(fields2, basestring): fields2 = [fields2]

    fig, axs = plt.subplots(len(fields2), len(fields1),
                            figsize=[8 * len(fields1),
                                     6.0875 * len(fields2)])
    axs = np.array(axs, ndmin=2)
    for if1, field1 in enumerate(fields1):
        for if2, field2 in enumerate(fields2):
            plt.sca(axs[if2][if1])
            if field1 == 'neu_frac':
                field1 = 'del_frac'
                inv1 = True
            else:
                inv1 = False
            if field2 == 'neu_frac':
                field2 = 'del_frac'
                inv2 = True
            else:
                inv2 = False

            x = {field1: tbs[field1].copy(), field2: tbs[field2].copy()}
            extr = {field: np.array([tb_raw[field].min(), tb_raw[field].max()]) for field in x.iterkeys()}

            if inv1:
                x[field1] = extr[field1][1] - x[field1]
                extr[field1] = extr[field1][1] - extr[field1][::-1]
            if inv2:
                x[field2] = extr[field2][1] - x[field2]
                extr[field2] = extr[field2][1] - extr[field2][::-1]

            labels = {f: f for f in x.iterkeys()}
            for field in x.iterkeys():
                if 'eff' in field:
                    x[field] *= 2
                    extr[field] *= 2
                if field != 'del_frac':
                    x[field] = np.log10(x[field])
                    extr[field] = np.log10(extr[field])
                    labels[field] = 'log10('+re.sub(r'_', r' ', field)+')'
                else:
                    labels[field] = re.sub(r'_', r' ', field)
            
            bins_x = np.linspace(extr[field1][0], extr[field1][1], nbins)
            bins_y = np.linspace(extr[field2][0], extr[field2][1], nbins)
            H, bins_x, bins_y = np.histogram2d(x[field1], x[field2],
                                               bins=(bins_x, bins_y),
                                               normed=True)
            plt.contourf(H.T, lw=2)
            #plt.matshow(H.T)
            plt.xticks(np.arange(len(bins_x) - 1),
                       map(lambda x: r'$'+str(np.round(x, 1))+r'$', 0.5 * (bins_x[:-1] + bins_x[1:])))
            plt.yticks(np.arange(len(bins_y) - 1),
                       map(lambda x: r'$'+str(np.round(x, 2))+r'$', 0.5 * (bins_y[:-1] + bins_y[1:])))
            plt.gca().xaxis.set_label_coords(0.5, -0.11)
            plt.gca().yaxis.set_label_coords(-0.12, 0.5)
            if if2 == 0:
                plt.xlabel(labels[field1])
            if if1 == 0:
                plt.ylabel(labels[field2])
            plt.colorbar()
    
    plt.tight_layout(w_pad=0.05, rect=(0.01, -0.01, 1.05, 0.99))

    plt.ion()
    plt.show()

    return H, bins_x, bins_y


def plot_marginals(tb_raw, tbss, fields, axsm=None):
    '''Plot the marginals of all the fields selected'''
    if not len(fields):
        raise ValueError('What should I plot?')

    labels = ['prior', 'both', 'Pfix', 'syndiv gtr x', 'syndiv les x']
    labels_labels = {'prior': 'prior',
                     'both': 'gated',
                     'Pfix': 'Area of Pfix',
                     'syndiv gtr x': r'syn diversity $> '+str(poly_synmin)+r'$',
                     'syndiv les x': r'syn diversity $< '+str(poly_synmax)+r'$'}

    if axsm is None:
        # Define the geometry of the plot
        n_rows = [1, 1, 1, 2, 1, 2, 2, 2]
        n_cols = [1, 2, 3, 2, 5, 3, 4, 4]
        fig, axs = plt.subplots(n_rows[len(fields) - 1], n_cols[len(fields) - 1],
                                figsize=(6 * n_cols[len(fields) - 1],
                                         5 * n_rows[len(fields) - 1]))
        axsm = np.array(axs.ravel(), ndmin=1)

    label_fields = {'neu_frac': 'neutral fraction of syn',
                    'del_frac_in': 'deleterious fraction of syn',
                    'ada_rate_in': 'rate of new epitopes',
                    'ada_eff_in': 'escape rate',
                    'coi': 'recombination rate',
                    'mu': 'mutation rate',
                    'del_eff_in': 'deleterious effect of syn',
                    'N': 'population size'}

    # Plot all fields
    h0 = []
    for i, field in enumerate(fields):
        ax = axsm[i]
        plt.sca(ax)

        label = label_fields[field]

        if field == 'neu_frac':
            field = 'del_frac'
            inv = True
        else:
            inv = False

        # Plot all partitions
        for j, label_tbs in enumerate(labels):
            x = tbss[label_tbs][field].copy()
            extr = np.array([tb_raw[field].min(), tb_raw[field].max()])
            if 'eff' in field:
                x *= 2
                extr *= 2
            if 'del_frac' not in field:
                x = np.log10(x)
                extr = np.log10(extr)

            # Use recombination instead of coinfection
            if field == 'coi':
                x -= 3
                extr -= 3

            if inv:
                x = extr[1] - x
                extr = extr[1] - extr[::-1]
    
            bins = np.linspace(extr[0], extr[1], 8)
            h = np.histogram(x, bins=bins, density=True)

            bins_c = 0.5 * (bins[1:] + bins[:-1])
            if 'del_frac' not in field:
                bins_c = 10**(bins_c)
    
            # Store the prior
            if j == 0: h0.append(h[0])

            # Plot
            plt.plot(bins_c, 1.0 * h[0],# / (1e-10 + h0[i]),
                        lw=2, c=cm.jet(int(255.0 * j / len(tbss))),
                        label=labels_labels[label_tbs])

        ax.set_xlabel(label)
        if 'del_frac' not in field:
            ax.set_xscale('log')
            ax.set_xlim(10**(extr[0]) * 0.9, 10**(extr[1]) * 1.1)

        if ax == axsm[-1]:
            ax.legend(bbox_to_anchor=(2.2, 0.7))
        
        plt.tight_layout(w_pad=0.25, rect=(0.01, 0, 0.99, 0.95))



# Script
if __name__ == '__main__':


    # Input results dir
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please specify a results folder')
    results_dir = args[1].rstrip('/')+'/'

    # Get plot type
    if len(args) > 2:
        plot_type = args[2]
        if plot_type == 'noplot':
            plot_things = False

    # Check for different generation times
    if 'gentime' in results_dir:
        tmp = os.path.basename(results_dir.rstrip('/')).split('_')
        i = tmp.index('gentime')
        gen_time = float(tmp[i+1])
    else:
        gen_time = 1.0

    # Check the existence of the summary file, if not present generate it
    if recondense or (not os.path.isfile(results_dir+datafile)):
        condense_data(results_dir)

    # Load data using record arrays
    tb_raw = np.genfromtxt(results_dir+datafile, unpack=False,
                       dtype=zip(('del_eff', 'ada_eff',
                                  'del_frac', 'ada_rate',
                                  'trashed', 'good',
                                  'Asyn', 'Anonsyn',
                                  'count_syn', 'counts_nonsyn',
                                  'subs_syn', 'subs_nonsyn',
                                  'poly_syn', 'poly_nonsyn',
                                  'N', 'mu', 'coi',
                                  'del_eff_in', 'ada_eff_in',
                                  'del_frac_in', 'ada_rate_in',
                                  'foldername'),
                                 (float, float,
                                  float, float,
                                  int, int,
                                  float, float,
                                  float, float,
                                  float, float,
                                  float, float,
                                  float, float, float,
                                  float, float,
                                  float, float,
                                  'S200')))

    # Filter out NaNs
    is_nonnan = -np.isnan(tb_raw['Asyn'])
    tb = tb_raw[is_nonnan]

    ## FIXME: Take only simulations with a high adaptive rate
    #tb_raw = tb_raw[tb_raw['ada_rate_in'] < 2e-3]

    # Make sure it's a 1D array (missing feature in genfromtxt)
    if not len(tb.shape):
        tb = np.array(tb, ndmin=1)

    ## Criterion: simulations with more than one sweep
    #is_twoplus_sweeps = (tb['ada_rate'] < 1) & (tb['ada_rate'] > 0)

    # Criterion: Asyn in the right window
    is_good_Pfix =  (tb['Asyn'] < Asynmax) & (tb['Asyn'] > Asynmin)

    # Criteria: window of synonymous diversity
    is_good_syndiv_gtr = (tb['poly_syn'] > poly_synmin)
    is_good_syndiv_les = (tb['poly_syn'] < poly_synmax)

    # Divide in classes
    is_good = is_good_Pfix & is_good_syndiv_gtr & is_good_syndiv_les
    tbg = tb[is_good]

    if plot_things:

        # Plot the gated results
        # Do we want a scatter plot?
        if plot_type == 'scatter':
            fig, axs = plt.subplots(1, 1, figsize=[8, 6.0875])
            colors = cm.jet([int(255.0 * (i + 2.5)) for i in np.log10(tbg['ada_eff'] / gen_time)])
            sizes = 10 + 120 * (np.log10(tbg['ada_rate'] / gen_time) + 3)
            ind = sizes.argsort()[::-1]
            plt.scatter(tbg['del_eff'][ind] * 2.0, tbg['del_frac'][ind],
                        color=colors[ind], s=sizes[ind])
            plt.xscale('log')
            plt.xlim(2 * 10**(-4.1) * gen_time, 2 * 10**(-1.9) * gen_time)
            plt.ylim(0.72, 1.02)
            plt.text(2 * 10**(-2.25) * gen_time, 0.735, 'flat prior', color='red')
            plt.gca().xaxis.set_label_coords(0.5, -0.07)
            plt.gca().yaxis.set_label_coords(-0.1, 0.5)

            # Draw a rectangle around the parameter priors
            rec = Rectangle(xy=(2 * 10**(-4) * gen_time, 0.75),
                            width=2 * 0.0099 * gen_time,
                            height=0.25,
                            facecolor='none', edgecolor='red',
                            lw=1.5)
            plt.gca().add_patch(rec)

            plt.xlabel(r'deleterious effect')
            plt.ylabel(r'deleterious fraction')

        # Heatmap otherwise
        elif plot_type == 'heatmap':
            bins_x = np.linspace(-4, -2, 8) + np.log10(2)
            bins_y = np.linspace(0.25, 1, 8)
            H, bins_x, bins_y = np.histogram2d(np.log10(tbg['del_eff'] * 2), tbg['del_frac'],
                                               bins=(bins_x, bins_y), normed=True)
            fig, axs = plt.subplots(1, 1, figsize=[9, 7.8])
            plt.imshow(H.T, interpolation='nearest')
            plt.xlim(-0.5, len(bins_x) - 1.5)
            plt.ylim(-0.50, len(bins_y) - 1.5)
            plt.gca().xaxis.set_label_coords(0.5, -0.12)
            plt.gca().yaxis.set_label_coords(-0.13, 0.5)
            plt.yticks(np.arange(len(bins_y)) - 0.5, map(lambda x: r'$'+str(np.round(x, 2))+r'$', bins_y))
            plt.xticks(np.arange(len(bins_x)) - 0.5, map(lambda x: r'$'+str(-np.round(x, 1))+r'$', bins_x))
            plt.colorbar()
            plt.xlabel(r'$\log_{10}$ (deleterious effect)')
            plt.ylabel(r'deleterious fraction')

        # Contour plot?
        elif plot_type == 'contour':
            plot_contour(tb_raw, tbg, 'del_eff', 'neu_frac', nbins=8)
            plt.xlabel(r'deleterious effect')
            plt.ylabel(r'neutral fraction of syn')

            # Restore log ticks on x axis
            import math
            field = 'del_eff_in'
            xm = math.ceil(10**(-3.6) * 1e4) / 1e4
            xM = math.floor(10**(-1.8) * 1e2) / 1e2
            xs = [n for n in np.arange(1, 10) * 1e-4 if n >= xm] + list(np.arange(1, 10) * 1e-3) + [xM]
            xs = np.array(xs)
            posx = (np.log10(xs) + 3.6) * plt.xlim()[1] / (3.6 - 1.8)
            ind_major = np.array([True if set(ii) == set(['0', '.', '1'])
                                   else False for ii in  map(str, xs)], bool)
            xlabels = [r'$10^{'+'{:1.0f}'.format(np.log10(xi))+r'}$' for xi in xs[ind_major]]
            ax = plt.gca()
            ax.set_xticks(posx[ind_major])
            ax.set_xticks(posx[-ind_major], minor=True)
            ax.set_xticklabels(xlabels, y=-0.01)


        #plt.title(('color is adaptive effect (2e-2.5 to 2e-1.5)\n'+
        #           'size is rate of new epitopes (1e-3 to 1e-2)'), fontsize=14)
        panel = plt.figtext(0.02, 0.92, panel_name, fontsize=24)

        if savefig:
            for ext in ['svg', 'pdf']:
                plt.savefig('/home/fabio/publications/synmut/figures/'+fig_name+'.'+ext)

        # Plot marginals
        fields = ['ada_rate_in', 'ada_eff_in', 'del_eff_in', 'neu_frac', 'N', 'mu', 'coi']
        tbss = {'prior': tb,
                'both': tbg,
                'Pfix': tb[is_good_Pfix],
                'syndiv gtr x': tb[is_good_syndiv_gtr],
                'syndiv les x': tb[is_good_syndiv_les],
               }

        # Make figure
        figm = plt.figure(figsize=(18, 11.5))
        axsm = []
        for i in xrange(2):
            for j in xrange(4):
                if i + j < 4:
                    axsm.append(plt.subplot(2, 4, 4 * i + j + 1))
        plot_marginals(tb_raw, tbss, fields, axsm=axsm)
        figm.suptitle('Marginal distributions of simulation parameters',
                      fontsize=20)
        if savefig:
            for ext in ['svg', 'pdf']:
                plt.savefig('/home/fabio/publications/synmut/figures/'+figm_name+'.'+ext)


        # Plot correlation matrices
        #import scipy.stats as stats
        #cmats = {}
        #Ps = {}
        #for j, label in enumerate(labels):
        #    tbs = tbss[label]
        #    a = np.zeros((len(tbs), len(fields)))
        #    for i, f in enumerate(fields):
        #        a[:, i] = tbs[f]
        #    cmat, P = stats.spearmanr(a)
        #    cmats[label] = cmat
        #    Ps[label] = P
        #    plt.figure(figsize=(7, 5.5))
        #    plt.imshow(cmat, interpolation='nearest', vmin=-1, vmax=1)
        #    plt.xticks(np.arange(len(fields)), [re.sub(r'_', ' ', f) for f in fields], rotation=45)
        #    plt.yticks(np.arange(len(fields)), [re.sub(r'_', ' ', f) for f in fields], rotation=45)
        #    plt.title(labels[j], fontsize=20)
        #    plt.tight_layout(rect=(0, 0, 1, 1))
        #    plt.colorbar()

        #    #plt.savefig('/home/fabio/publications/synmut/figures/simulations_rankcorr_'+labels[j]+'.pdf')


        #########################################################################
        ## FIXME
        ## Plot marginals for simulations with or without 2+ sweeps
        #########################################################################
        #fields = ['del_eff', 'del_frac', 'N', 'mu', 'coi']
        #labels = ['prior', '2+ sweeps', '1 sweep', 'no sweep']
        #tb_raws = [tb_raw,
        #           tb_raw[(tb_raw['ada_rate'] > 0) & (tb_raw['ada_rate'] < 1)],
        #           tb_raw[(tb_raw['ada_rate'] < 0) | (tb_raw['ada_rate'] > 1)],
        #           tb_raw[tb_raw['ada_rate'] == 0],
        #          ]
        #for j, tbs in enumerate(tb_raws):
        #    labels[j] = labels[j]+': '+str(len(tbs))
        #pran = np.zeros(2, dtype=[(f, float) for f in fields])
        #for f in fields:
        #    pran[f] = (tb_raw[f].min(), tb_raw[f].max())
        #figm, axsm = plt.subplots(1, 5, figsize=(18, 5.5))
        #axsm = axsm.ravel()
        #lines = []
        #h0 = []
        #for i, field in enumerate(fields):
        #    ax = axsm[i]
        #    for j, tbs in enumerate(tb_raws):
        #        a = tbs[field]
        #        extr = pran[field]
        #        if field != 'del_frac':
        #            a = np.log10(a)
        #            extr = np.log10(extr)
        #            # Fitness effects must be multiplied by 2 because of the -1/1 basis
        #            if 'eff' in field:
        #                a += np.log10(2)
        #                extr += np.log10(2)
        #            # Use recombination rate instead of coinfection
        #            if field == 'coi':
        #                label = r'$\log_{10}$ (recombination rate)'
        #                a -= 3
        #                extr -= 3
        #            else:
        #                label = r'$\log_{10}$ ('+re.sub('_', ' ', field)+r')'
        #        else:
        #            label = re.sub('_', ' ', field)
        #        if False:#field == 'ada_rate':
        #            bins = np.linspace(-3.2, -1.7, 10)
        #        else:
        #            bins = np.linspace(extr[0], extr[1], 8)
        #        h = np.histogram(a, bins=bins,
        #                         density=True)
        #        # Store the prior
        #        if j == 0:
        #            h0.append(h[0])
        #        l = ax.plot(0.5 * (h[1][1:] + h[1][:-1]), 1.0 * h[0],# / (1e-10 + h0[i]),
        #                    lw=2, c=cm.jet(int(255.0 * j / len(tb_raws))),
        #                    label=labels[j])
        #        if i == 0:
        #            lines.append(l[0])
        #    ax.set_xlabel(label)

        #axsm[-1].legend(lines, labels, bbox_to_anchor=(2.42, 0.8))
        #axsm[0].set_xticks([-3.6, -3.3, -3.0, -2.7, -2.4, -2.1, -1.8])

        #plt.tight_layout(w_pad=0.05, rect=(0, 0, 0.85, 0.95))


        #
        plt.ion()
        plt.show()
