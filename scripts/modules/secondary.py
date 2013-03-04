# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/10/12
content:    Module for parsing and handling RNA secondary structures
'''
# Modules
import os
import numpy as np



# Classes
class entry(object):
    '''RNA secondary structure prediction entry'''
    
    def __init__(self, name, seq, contact_map=None, vienna=None):
        self.name = name
        self.seq = seq
        if contact_map is not None:
            self.contact_map = contact_map
            self.vienna = contact_map_to_vienna(contact_map, length=len(seq))
        elif vienna is not None:
            self.vienna = vienna
            self.contact_map = vienna_to_contact_map(vienna)
        else:
            raise ValueError('Please input a Vienna format or a contact map!')


    def adapt_map(self, sequence):
        '''Adapt the prediction to a gapped sequence'''
        vv = []
        ii = 0
        for s in sequence:
            if s == '-':
                vv.append('.')
            else:
                vv.append(self.vienna[ii])
                ii += 1
        return vienna_to_contact_map(''.join(vv))



# Functions
def vienna_to_contact_map(vienna):
    '''Convert a Vienna format to contact map'''
    contact_map = {}
    tmp = []
    for j, c in enumerate(vienna):
        if c == '(':
            tmp.append(j)
        elif c == ')':
            p1 = tmp.pop(-1)
            p2 = j
            contact_map[p1] = p2
            contact_map[p2] = p1
    return contact_map


def contact_map_to_vienna(contact_map, length):
    '''Convert a contact map into Vienna format'''

    if max(contact_map.keys()) >= length:
        raise ValueError('Length is smaller than maximal contact')

    vv = np.repeat('.', length)
    for (p1, p2) in contact_map.items():
        if p1 < p2:
            vv[p1] = '('
        else:
            vv[p1] = ')'

    return ''.join(vv)


def read_Vienna(filename):
    '''Parse a list of sequence-structures in Vienna format'''

    entries = []

    # Check file existance
    with open(filename, 'r') as f:
        for i, line in enumerate(f):

            # Read name
            if i % 3 == 0:
                name = line.lstrip('> ')[:-1]

            # Read sequence
            elif i % 3 == 1:
                seq = line[:-1]

            # Read map and instantiate entry class
            else:
                vienna = line.split(' ')[0]
                entries.append(entry(name, seq, vienna=vienna))

    return entries    
