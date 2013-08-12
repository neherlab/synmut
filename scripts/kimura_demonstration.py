import matplotlib.pyplot as plt
import numpy as np

dele_frac = np.linspace(0.05,0.95,9)
Ts = np.linspace(-10,2,24)

def fix_prob(s, af):
    return (1-np.exp(-s*af))/(1-np.exp(-s))

def area(sd, alpha, af_vec=np.linspace(0.05,0.95,19)):
    sfs_neutral = (1-alpha)/af_vec
    sfs_dele  = alpha*np.exp(sd)/af_vec

    frac_fixed = (af_vec*sfs_neutral+sfs_dele*fix_prob(sd,af_vec))\
        /(sfs_neutral+sfs_dele)

    return sum(frac_fixed/len(af_vec)), frac_fixed

A = np.zeros((len(dele_frac), len(Ts)))
for ai,alpha in enumerate(dele_frac):
    for si,sd in enumerate(Ts):
        A[ai,si], fp = area(sd,alpha)

    if ai%2:
        plt.plot(Ts, A[ai,:]-0.5)
    else:
        plt.plot(Ts, A[ai,:]-0.5, label = 'dele frac = '+str(alpha))

plt.legend(loc=2)
plt.xlabel('Ts')
plt.ylabel('Area deviation')
plt.ion()
plt.show()
