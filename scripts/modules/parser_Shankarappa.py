# vim: fdm=marker:fileencoding=utf-8
'''
author:     Fabio Zanini
date:       04/01/2012
content:    Parser for various Shankarappa data files
'''
# Standard modules
import re
import os
import glob
import shutil
from operator import *
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as ppl
import matplotlib.cm as cm

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import *
from Bio.Align import MultipleSeqAlignment as MSA

import longitudinal.classes.parser_SHAPE as pSHAPE
from longitudinal.classes.helper import find_with_gaps
from longitudinal.classes.parser_reference import seq_patterns
from longitudinal.classes.parser_generic import patient_generic, alpha


# Globals
datadir = '/home/fabio/university/phd/longitudinal/data/Shankarappa'


# Classes
class clinical_parser(object):
    '''Parse the TXT file with the clinical info.'''
    def __init__(self, filename):
        cols = [['visit', 1, int],
                ['mo_to_SC', 2, float],
                ['logRNA', 4, float],
                ['CD3', 5, float],
                ['CD8', 6, float],
                ['CD4', 7, float]]

        with open(filename, 'r') as f:
            f.readline()    # Header
            lsts = [l.split('\t') for l in f.readlines()]

        names = set(map(itemgetter(0),lsts)) - set([''])
        patients = [patient(p) for p in names]

        for p in patients:
            tmpls = filter(lambda x: x[0] == p.name,lsts)
            n_points = len(tmpls)
            tmpdict = {}
            for n in cols:
                tmpdict[n[0]] = np.ma.zeros(n_points, dtype=n[2])

            for i,tmp in enumerate(tmpls):
                for n in cols:
                    try:
                        tmpdict[n[0]][i] = tmp[n[1]]
                    except ValueError:
                        tmpdict[n[0]][i] = np.ma.masked
                    
            for n,c in tmpdict.items():
                setattr(p,n,c)

            # Try to add days_from_seroconversion
            p.days_from_seroconversion = p.mo_to_SC * 30.5

        patients.sort(key=(lambda x: int(x.name[1:])))
        self.patients = patients


    def plot(self,patients=[]):
        '''Plot the clinical data as a function of time.'''

        if not patients:
            patients = self.patients
        else:
            patients = list(np.array(self.patients)[patients])

        ppl.ioff()
        for p in patients:
            p.plot()
        ppl.ion()
        ppl.show()



class patient(patient_generic):
    def __init__(self,name):
        self.name = name
        self.visit = []
        self.seqs = []
        self.file_tree = datadir+'/tree_nj_'+self.name+'.txt'


    def __call__(self):
        for a in vars(self).keys():
                if a != 'seqs':
                    print a
                    print getattr(self,a)
                elif self.seqs:
                    print a
                    print self.seqs[0]
                    print '[...]'
                print


    def filter_for_longitudinal_genetics(self):
        '''Filter only time points with sequences.'''
        self.filter_only_sequenced()


    def visit_to_time(self,visitnumber):
        return self.mo_to_SC[self.visit == visitnumber][0]


    def time_to_visit(self,mo_to_SC):
        return self.visit[self.mo_to_SC == mo_to_SC][0]


    def seqs_from_visit(self,visitnumbers):
        if np.isscalar(visitnumbers):
            visitnumbers = [visitnumbers]
        return MSA(filter(lambda x: x.visit in visitnumbers, self.seqs))


    def get_names(self):
        return get_names(self.seqs)


    def get_visits(self):
        return get_visits(self.seqs)


    def plot_clinical(self, ax=[], show=False):

        if not ax:
            fig = ppl.figure()
            ax = fig.add_subplot(1,1,1)

        ax.scatter(self.mo_to_SC,self.CD8,c='b',s=40,label='CD8')
        ax.scatter(self.mo_to_SC,self.CD4,c='g',s=40,label='CD4')
        ax.scatter(self.mo_to_SC,self.CD3,c='y',s=40,label='CD3')
        ax.scatter(self.mo_to_SC,10**self.logRNA,c='r',s=40,label='RNA')
        
        ax.plot([self.mo_to_SC.min(),self.mo_to_SC.max()],[200,200],'--k',lw=2)
        
        ax.set_yscale('log')
        ax.set_title(self.name,fontsize=20)
        ax.set_xlabel('Months to SC')

        ppl.legend(loc=2)

        if show:
            ppl.ion()
            ppl.show()


    def filter_only_sequenced(self):
        ind = -np.ma.getmaskarray(self.visit)
        seqs = self.seqs
        if seqs:
            for i,v in enumerate(self.visit):
                if ind[i]:
                    seqs = self.seqs_from_visit(v)
                    ind[i] = (len(seqs) > 0)

        for attr in ['visit','mo_to_SC','logRNA','CD3','CD4','CD8']:
            setattr(self,attr,getattr(self,attr)[ind])


    def find_V_region(self, name):
        '''Find a variable region based on the reference'''
        if name not in self.seq_patterns:
            raise ValueError(str(name)+" is not in self.seq_patterns.")
        
        patterns = self.seq_patterns[name]
        seq = str(self.reference_seq.seq)
        start = find_with_gaps(seq, patterns[0])
        end = find_with_gaps(seq[start:], patterns[1])
        if end != -1:
            end += start 
        # V5 convers the final part of the seqs
        if name == 'V5':
            end = len(seq)
        return [start, end]


    def show_V_region(self, name):
        '''Show a variable region'''
        indices = self.find_V_region(name)
        if (-1) not in indices:
            return self.reference_seq[indices[0]: indices[1]]
        else:
            return None


    def only_V_regions(self):
        '''Denote the hypervariable regions'''

        cons = self.get_consensus_initial()
        is_V = np.zeros(len(cons), bool)
        for name in self.seq_patterns:
            indices = self.find_V_region(name)
            if (-1) not in indices:
                is_V[indices[0]:indices[1]] = True
        return is_V


    def only_region(self, region_name):
        '''Start and end indices for the chosen region in the patient MSA.

        Note: this function can only be used with a reference sequence for which
              the known region patterns exist.
        
        Note: the region can be selected by calling

           (s, e) = only_region(...)
           self.seqs[:, s:e]
        '''
        if not hasattr(self, 'reference_seq'):
            raise AttributeError('No reference sequence found.')

        from longitudinal.classes.helper import find_with_gaps
        from longitudinal.classes.parser_reference import seq_patterns

        # Get the nucleotide patterns for the chosen region and reference
        patterns = seq_patterns[self.reference_name][region_name]

        # Get the start and end indices
        start = find_with_gaps(str(self.reference_seq.seq), patterns[0])
        end = find_with_gaps(str(self.reference_seq.seq), patterns[1])
        if (start == -1) and (end == -1):
            raise ValueError('Neither start nor end pattern found, take or leave')

        if start == (-1):
            start = 0
        if end == (-1):
            end = len(self.reference_seq)
        return (start, end)


    def coallele_frequencies(self, i1, a1, i2, a2):
        '''Monitor double sweeps of a1 at position i1 and a2 at i2'''
    
        visits = self.visit
        Lv = len(visits)
    
        # 4 places: WT, SM1, SM2, DM
        freqs = np.zeros((4,Lv))
    
        for j,v in enumerate(visits):
            seqs = np.array(self.seqs_from_visit(v))[:,[i1,i2]]
    
            ind1 = (seqs[:,0] == a1)
            ind2 = (seqs[:,1] == a2)
    
            freqs[0,j] = ((-ind1) * (-ind2)).mean()
            freqs[1,j] = (ind1 * (-ind2)).mean()
            freqs[2,j] = ((-ind1) * ind2).mean()
            freqs[3,j] = (ind1 * ind2).mean()
    
        return freqs


    def extend_sequences(self, seqs=None, alpha=np.array(['A', 'C', 'G', 'T', '-'])):
        '''Extend the seqs in their largest binary space.'''
        if seqs == None:
            seqs = self.seqs

        seqs = np.array(seqs)

        # Filter out gaps and stuff
        is_nongap = ((seqs == '-').mean(axis=0) == 0)
        if is_nongap.mean() < 0.5:
            raise ValueError('Most positions are recognized as gaps!?')
        seqs = seqs[:,is_nongap]
    
        # Create a tensor of rank 3
        sext = np.zeros(((alpha != '-').sum(),seqs.shape[0],seqs.shape[1]),bool)
        for j,a in enumerate(alpha[alpha != '-']):
            sext[j] = (seqs == a)
        # Swap axes
        sext = sext.swapaxes(1,0).swapaxes(1,2)
        # Flattify the last dimension
        sext = sext.reshape((sext.shape[0],sext.shape[1] * sext.shape[2]))
        return sext


    def PCA(self):
        '''Calculate PCA of the patient seqs in the extended space.'''
        # Extend (kernelize)
        X = self.extend_sequences()

        # Calculate PCA
        U,l,Vh = linalg.svd(X)
        self.U = U
        self.l = l
        self.Vh = Vh


    def plot_PCA(self, dimension='2D', **kwargs):
        if dimension in ['2','2D','2d']:
            self._plot_PCA_2D(**kwargs)
        else:
            self._plot_PCA_3D(**kwargs)


    def _plot_PCA_2D(patient, ax=None, show=False):
        try:
            self.U
        except AttributeError:
            raise AttributeError('Please perform PCA before plotting.')
    
        # Time points
        visits = np.array(map(attrgetter('visit'), self.seqs))
    
        # Plot different time points with different colors
        n = len(self.visit)
        cols = cm.jet([int(255.0 * i / n) for i in xrange(n)])
        
        if ax == None:
            fig = ppl.figure(figsize=(6,5))
            ax = fig.add_subplot(1,1,1)
    
        for i,v in enumerate(self.visit):
            ind = (visits == v)
            
            x = self.U[ind,0]
            y = self.U[ind,1]
        
            ax.scatter(x,y,color=cols[i], label='v '+str(v))
        
        ax.set_xlabel('PC1', fontsize=18)
        ax.set_ylabel('PC2', fontsize=18)
        ax.set_title('Patient '+self.name)
    
        if show:
            ppl.show()


    def _plot_PCA_3D(self, ax=None, show=False):
        try:
            self.U
        except AttributeError:
            raise AttributeError('Please perform PCA before plotting.')

        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError('3D extension to matplotlib not found.')
    
        # Time points
        visits = np.array(map(attrgetter('visit'), self.seqs))
    
        # Plot different time points with different colors
        n = len(self.visit)
        cols = cm.jet([int(255.0 * i / n) for i in xrange(n)])
        
        if ax == None:
            fig = ppl.figure(figsize=(6,5))
            ax = fig.add_subplot(1,1,1, projection='3d')
    
        for i,v in enumerate(self.visit):
            ind = (visits == v)
            
            x = self.U[ind,0]
            y = self.U[ind,1]
            z = self.U[ind,2]
        
            ax.scatter(x,y,z, color=cols[i])
        
        ax.set_xlabel('PC1', fontsize=18)
        ax.set_ylabel('PC2', fontsize=18)
        ax.set_zlabel('PC3', fontsize=18)
        ax.set_title('Patient '+self.name)
    
        if show:
            ppl.show()


    def plot_tree_NJ(self, ax=None, show=False):
        '''Plot NJ tree'''
        import matplotlib
        from Bio import Phylo
        from longitudinal.classes.phylo_draw import draw

        # Keep only sequenced time points
        self.filter_only_sequenced()

        # Get the tree
        tree = Phylo.read(self.file_tree, "newick")

        # Root tree at one of the initial seqs
        names_init = get_names(self.seqs_from_visit(self.visit[0]))
        leaves = tree.get_terminals()
        leafnames = map(attrgetter('name'), leaves)
        root = leaves[leafnames.index(names_init[0])]
        tree.root_with_outgroup(root)
        tree.ladderize()

        if ax is None:
            fig = ppl.figure()
            ax = fig.add_subplot(1,1,1)

        draw(tree, do_show=False, axes=ax)
        ax.set_title('NJ tree of '+str(self),fontsize=18)
    
        # Change the leaves' colors depending on time
        visits = list(self.visit)
        n = len(visits)
        cols = [cm.jet(int(255.0 * i / n)) for i in xrange(n)]
    
        leaves = filter(lambda x: isinstance(x, matplotlib.text.Text), ax.get_children())
    
        # Find the visit number of each leaf
        names = list(self.get_names())
        for l in leaves:
            text = l.get_text().strip(' ')
            if text in names:
                v = self.seqs[names.index(text)].visit
                try:
                    ind = visits.index(v)
                    l.set_color(cols[ind])
                except ValueError:
                    # Some seqs are weirdly labelled (!)
                    l.set_text(text+' BAD!')
        if show:
            ppl.show()


    def parse_sequences(self, reference=None, **kwargs):
        '''Parse sequences, aligned already.'''
        if reference is not None:
            if not isinstance(reference, str):
                refstr = '+'.join(reference)
                multiple_ref = True
            else:
                refstr = reference
                multiple_ref = False
            self.file_sequences = datadir+'/to_reference/'+self.name[1:]+'_to_'+refstr+'_aligned.fasta'
        else:
            self.file_sequences = datadir+'/nuc_aligned_'+self.name[1:]+'.fasta'

        alignment = AlignIO.read(self.file_sequences,
                                 format='fasta',
                                 alphabet=ambiguous_dna)
        if reference is not None:
            if multiple_ref:
                self.reference_names = reference
                self.reference_seqs = alignment[:len(reference)]
                alignment = alignment[len(reference):]

                # Identify variable loops
                for refstrtmp in reference:
                    if refstrtmp in seq_patterns:
                       self.seq_patterns = seq_patterns[refstrtmp]

            else:
                self.reference_name = refstr
                self.reference_seq = alignment[0]
                alignment = alignment[1:]

                # Identify variable loops
                if refstr in seq_patterns:
                    self.seq_patterns = seq_patterns[refstr]

            if 'SHAPE' in refstr:
                if multiple_ref:
                    for i in xrange(len(reference)):
                        if 'SHAPE' in self.reference_seqs[i].name:
                            SHAPEseq = str(self.reference_seqs[i].seq)
                else:
                    SHAPEseq = str(self.reference_seq.seq)
                (self.SHAPE_p, self.SHAPE_r) = pSHAPE.assign_prob(SHAPEseq,
                                                                  include_reactivity=True)

                # Flag for excluding the map, which takes a while
                if 'SHAPEmap' in kwargs and kwargs['SHAPEmap']:
                    self.SHAPEmap = pSHAPE.assign_map(SHAPEseq)

        for rec in alignment:
            rec.visit = int(re.findall('V\d{2}',rec.name)[0][1:])
        self.seqs = alignment



# Functions
def parse_sequences(reference=None, exclude=None, **kwargs):
    '''Parse Shankarappa sequences using clinical data in the background.
    
    Parameters:
        - reference: use one/several reference sequences?
        - exclude: ignore certain patients
    '''

    clin = clinical_parser(datadir+'/set1_clinical_data.txt')
    patients = clin.patients

    if exclude is not None:
        patients = filter(lambda x: x.name not in exclude, patients)

    for p in patients:
        p.parse_sequences(reference=reference, **kwargs)

    return patients


def align_sequences(translate=True, reference=None):
    '''Multiple Sequence Alignment using DNA or aa.
    
    Parameters:
    - translate is True if amino acids are to be aligned and the aligment to be
      back-propagated instead of a direct DNA alignment
    - reference is None if self-alignment is used, is a string if some reference
      sequence is to be used
    '''

    from Bio.Align.Applications import MuscleCommandline
    # Find the directory
    dirname = datadir

    # Find the unaligned files
    names = glob.glob(dirname+'/nuc_*.fasta')
    names = [f.rstrip('.fasta') for f in names]
    names = filter(lambda f: 'aligned' not in f, names)
    names = list(set(names))

    # Set reference stuff
    if not isinstance(reference, str):
        refstr = '+'.join(reference)
        multiple_ref = True
    else:
        refstr = reference
        multiple_ref = False

    # For each patient, do the alignment
    for name in names:
        print os.path.basename(name)

        # Check whether a reference is given
        if refstr is None:
            fileDNAin = name+'.fasta'
            fileDNAout = name+'_aligned.fasta'
            fileaain = name+'_aa.fasta'
            fileaaout = name+'_aligned_aa.fasta'

        # if there is a reference, we need to create new files which include the
        # reference in a special folder
        else:
            nname = os.path.basename(name)[4:]
            fileDNAin = dirname+'/to_reference/'+nname+'_to_'+refstr+'.fasta'
            fileDNAout = dirname+'/to_reference/'+nname+'_to_'+refstr+'_aligned.fasta'
            fileaain = dirname+'/to_reference/'+nname+'_to_'+refstr+'_aa.fasta'
            fileaaout = dirname+'/to_reference/'+nname+'_to_'+refstr+'_aligned_aa.fasta'

            # Create folder if not present
            if not os.path.isdir(dirname+'/to_reference'):
                os.mkdir(dirname+'/to_reference')

            # Copy stuff into the new folder
            if multiple_ref:
                with open(fileDNAin, 'w') as fDNAin:
                    for ref in reference:
                        reffile = dirname+'/reference/'+ref+'_trimmed.fasta'
                        with open(reffile, 'r') as fref:
                            fDNAin.write(fref.read())
            else:
                reffile = dirname+'/reference/'+str(reference)+'_trimmed.fasta'
                shutil.copy(reffile, fileDNAin)
            with open(name+'.fasta', 'r') as raw_file:
                with open(fileDNAin, 'a') as new_file:
                    for line in raw_file:
                        new_file.write(line)

        # Direct DNA alignment
        if not translate:       
            cline = MuscleCommandline(input=fileDNAin,
                                      out=fileDNAout,
                                      diags=True)
            cline()

        # Pass through the protein sequence
        else:
            from Bio.SeqRecord import SeqRecord

            # Translate
            DNArecs = list(SeqIO.parse(fileDNAin, 'fasta'))
            DNAids = map(attrgetter('id'), DNArecs)
            aarecs = []
            for rec in DNArecs:
                aarecs.append(SeqRecord(rec.seq.translate(), id=rec.id,
                                        name=rec.name,
                                        description=rec.description))
            SeqIO.write(aarecs, fileaain, 'fasta')

            # Align
            cline = MuscleCommandline(input=fileaain,
                                      out=fileaaout,
                                      diags=True)
            cline()

            # Back-translate
            aarecs = SeqIO.parse(fileaaout, 'fasta')
            bDNArecs = []
            for rec in aarecs:
                DNArec = DNArecs[DNAids.index(rec.id)]
                DNAseq = DNArec.seq
                bDNAseq = []
                i = 0
                for aa in str(rec.seq):
                    if aa != '-':
                        bDNAseq.extend(DNAseq[i:i+3])
                        i += 3
                    else:
                        bDNAseq.extend(['---'])
                bDNAseq = Seq(''.join(bDNAseq))
                bDNArec = SeqRecord(bDNAseq, id=rec.id,
                                    name=rec.name,
                                    description=rec.description)
                bDNArecs.append(bDNArec)
            SeqIO.write(bDNArecs, fileDNAout, 'fasta')


def plot_sweeps(ts,fsweeps,indices=[],labels=[],ax=[]):
    if not indices:
        indices = np.arange(len(fsweeps))
    Li = len(fsweeps)

    if not labels:
        labels = ['i = '+str(ind) for ind in indices]

    cols = [cm.jet(int(255.0 * i / Li)) for i in xrange(Li)]
    cols = np.array(cols)
    cols[:,3] = 0.5

    if not ax:
        fig = ppl.figure()
        ax = fig.add_subplot(1,1,1)
    for i,j in enumerate(indices):
            ax.plot(ts+0.1*j,fsweeps[j],lw=3, color=cols[j],label=labels[i])
    ax.set_xlabel('time [mo to SC]')
    ax.set_ylabel('f')
    ax.set_title('Sweeping loci',fontsize=18)

    ppl.legend(loc=5)
    ppl.ion()
    ppl.show()


def get_names(seqs):
    return np.array(map(attrgetter('name'), seqs))


def get_visits(seqs):
    return np.array(map(attrgetter('visit'), seqs),int)
