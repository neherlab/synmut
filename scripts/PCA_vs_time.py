# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/02/2012
content:    Show the first PCs vs time.
'''
# Standard modules
from operator import *
import numpy as np
import matplotlib.pyplot as ppl
import matplotlib.cm as cm
import scipy.linalg as linalg
from mpl_toolkits.mplot3d import Axes3D

# Custom modules
import modules.parser_Shankarappa as pS
from modules.alphabet import alpha



# Functions
def plot3D(patient, ax=None, show=False):
    if not hasattr(patient,'U'):
        raise AttributeError('Please perform PCA before plotting.')

    # Time points
    visits = np.array(map(attrgetter('visit'), patient.seqs))

    # Plot different time points with different colors
    n = len(patient.visit)
    cols = cm.jet([int(255.0 * i / n) for i in xrange(n)])
    
    if ax == None:
        fig = ppl.figure(figsize=(6,5))
        ax = fig.add_subplot(1,1,1, projection='3d')

    for i,v in enumerate(patient.visit):
        ind = (visits == v)
        
        x = patient.U[ind,0]
        y = patient.U[ind,1]
        z = patient.U[ind,2]
    
        ax.scatter(x,y,z, color=cols[i], s=60)
    
    ax.set_xlabel('PC1', fontsize=18)
    ax.set_ylabel('PC2', fontsize=18)
    ax.set_zlabel('PC3', fontsize=18)
    ax.set_title('Patient '+patient.name)

    if show:
        ppl.show()



# Script
if __name__ == '__main__':
    
    patients = pS.parse_sequences()

    for p in patients[0:1]:
        print p
        p.filter_only_sequenced()
    
        # Perform PCA on all sequences
        p.PCA()

        # Plot
        plot3D(p)

    ppl.ion()
    ppl.show()
