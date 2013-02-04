# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/07/12
content:    Parser for the SHAPE reactivity data, from doi:10.1038/nature08237.
'''
# Modules
import numpy as np
import re


# Globals
helixfile = '/home/fabio/university/phd/longitudinal/data/SHAPE/nature08237-s2.txt'
bpseqfile = '/home/fabio/university/phd/longitudinal/data/SHAPE/nature08237-s2_bpseq.txt'
ppfile = '/home/fabio/university/phd/longitudinal/data/SHAPE/nature08237-s3.txt'



# Classes
class secondary(object):
    '''Store the SHAPE reactivities and pairing probs.'''

    Shankarappa_indices = [6558, 7173]
    

    def __init__(self):
        '''Store sequence, reactivities, pairing probs, and map.'''
        self.r, self.p = np.genfromtxt(ppfile,
                                       usecols=(2,3),
                                       delimiter='\t',
                                       unpack=True,
                                       usemask=True)
    
        self.sequence = np.genfromtxt(ppfile,
                                      dtype='S1',
                                      usecols=[1],
                                      delimiter='\t',
                                      unpack=True,
                                      usemask=False)

        # Convert to DNA notation
        self.sequence[self.sequence == 'U'] = 'T'

        # Save map in a masked array
        bp1, bp2 = np.loadtxt(bpseqfile, unpack=True, usecols=(0,2), dtype=int)
        self.map = np.ma.masked_all(len(self.sequence), int)
        self.map[bp1 - 1] = bp2 - 1


    def __repr__(self):
        return 'secondary structure class of HIV-1 NL4-3'



# Functions
def assign_prob(seq, include_reactivity=False):
    '''Create a masked array with pairing probability and reactivity.
    
    Gap positions in the input sequence will stay masked.
    Sites where the reactivity has not been measured will stay masked as well.
    '''
    # Create secondary structure object
    sec = secondary()
    s = ''.join(sec.sequence)

    # Initialize
    assp = np.ma.masked_all(len(seq))
    assr = np.ma.masked_all(len(seq))

    # Short names
    p = sec.p
    r = sec.r

    # Check input string
    if not isinstance(seq, basestring):
        seq = ''.join(seq)

    # Take the ungapped sequence and see whether the input sequence is
    # compatible with that one (sanity check)
    start = s.find(re.sub('-','', seq))
    if start == -1:
        raise ValueError('Sequence not found!')

    # Fill results
    ii = start
    for i in xrange(len(assp)):
        if seq[i] != '-':
            assp[i] = p[ii]
            assr[i] = r[ii]
            ii += 1
    if include_reactivity:
        return (assp, assr)
    else:
        return assp


def assign_map(seq):
    '''Create a masked array with pairing map.
    
    Gap positions in the input sequence will stay masked. Non-paired positions
    will also be masked. Please look at the gaps in the sequence to distinguish
    between the two.

    Note: this algorithm is more involved than the one for pairing
    probabilities, because the second pairing partner is a index and must take
    into account gaps.

    Note: this function sets the pairing partner to -1 if it is outside the
    sequenced region.
    '''
    # Create secondary structure object
    sec = secondary()
    s = ''.join(sec.sequence)

    # Initialize
    ass = np.ma.masked_all(len(seq), int)
    attr = sec.map

    # Check input string
    if not isinstance(seq, basestring):
        seq = ''.join(seq)

    # Take the ungapped sequence and see whether the input sequence is
    # compatible with that one (sanity check)
    seqng = re.sub('-','',seq)
    start = s.find(seqng)
    if start == -1:
        raise ValueError('Sequence not found!')

    # Create index map, so that
    #     ind[nongapped] = gapped
    # and of course gapped >= nongapped
    ind = np.zeros(len(seqng), int)
    ii = start
    for i in xrange(len(seq)):
        if seq[i] != '-':
            ind[ii - start] = i
            ii += 1

    # Fill results
    ii = 0
    for i in xrange(len(seq)):
        if seq[i] != '-':
            bp2 = attr[start + ii]
            # Check whether there is any pairing
            if (bp2 is not np.ma.masked):
                # Check whether the pairing is WITHIN the sequenced region
                if (bp2 - start > 0) and (bp2 - start < len(seqng)):
                    ass[i] = ind[bp2 - start]
                # Otherwise set a special flag (-1)
                else:
                    ass[i] = -1
            ii += 1
    return ass    



def helix_to_bpseq():
    '''Convert the helix file into bpseq format'''
    # Create secondary structure object
    sec = secondary()

    # Load helix file
    h1, h2, hls = np.loadtxt(helixfile, unpack=True, dtype=int)

    # Total number of base pairs
    L = np.sum(hls)

    # Store bpseq
    bp1 = np.zeros(2 * L, int)
    bpa = np.zeros(2 * L, 'S1')
    bp2 = np.zeros(2 * L, int)
    ibp = 0
    for i in xrange(len(hls)):
        hl = hls[i]
        bp1[ibp: ibp + hl] = h1[i] + np.arange(hl)
        bp2[ibp: ibp + hl] = h2[i] - np.arange(hl)
        bpa[ibp: ibp + hl] = sec.sequence[bp1[ibp: ibp + hl] - 1]
        ibp += hl

        # Mirror image
        bp1[ibp: ibp + hl] = h2[i] - hl + 1 + np.arange(hl)
        bp2[ibp: ibp + hl] = h1[i] + hl - 1 - np.arange(hl)
        bpa[ibp: ibp + hl] = sec.sequence[bp1[ibp: ibp + hl] - 1]
        ibp += hl

    # Sort them by the first column
    ind_sort = np.argsort(bp1)
    bp1 = bp1[ind_sort]
    bp2 = bp2[ind_sort]
    bpa = bpa[ind_sort]

    # Write to output file
    with open(bpseqfile, 'w') as f:
        f.write('# RNA secondary structure pairing map from Watts et al. (2009), bpseq format')
        f.write('\n')
        f.write('# Reference sequence NL4-3')
        f.write('\n')
        for i in xrange(2 * L):
            f.write(str(bp1[i]))
            f.write('\t')
            f.write(bpa[i])
            f.write('\t')
            f.write(str(bp2[i]))
            f.write('\n')
