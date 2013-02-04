# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/02/2012
content:    Parser for the various files downloaded from the LANL HIV database.
'''
# Standard modules
import os
import re
import glob
from operator import *
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as ppl
import matplotlib.cm as cm

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import *

import parser_SHAPE as pSHAPE
from parser_generic import patient_generic, alpha



# Globals
datadir = '../data/longitudinal/LANL'

columns = [('#',int),
           ('patient_code',str),
           ('patient_id',str),
           ('accession',str),
           ('name',str),
           ('subtype',str),
           ('country',str),
           ('sampling_year',int),
           ('days_from_first_sample',int),
           ('Fiebig_stage',str),
           ('days_from_treatment_end',int),
           ('n_timepoints',int),
           ('days_from_treatment_start',int),
           ('days_from_infection',int),
           ('days_from_seroconversion',int),
           ('pubmed_ID',str),
           ('start',int),
           ('stop',int),
           ('sequence_length',int),
           ('organism',str),
           ('sequence',str)]


# Classes
class patient(patient_generic):
    def __init__(self,name):
        self.name = name
        self.visit = []
        self.seqs = []

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


    def visit_to_time(self,visitnumber):
        attr = filter(lambda x: 'days' in x, vars(self).keys())[0]
        return getattr(self, attr)[self.visit == visitnumber][0]


    def time_to_visit(self,days):
        attr = filter(lambda x: 'days' in x, vars(self).keys())[0]
        return self.visit[getattr(self, attr) == days][0]


    def seqs_from_visit(self,visitnumbers):
        if np.isscalar(visitnumbers):
            visitnumbers = [visitnumbers]
        return filter(lambda x: x.visit in visitnumbers, self.seqs)


    def get_ids(self):
        return get_ids(self.seqs)


    def get_names(self):
        return get_names(self.seqs)


    def get_visits(self):
        return get_visits(self.seqs)


    def add_features_from_seqs(self):
        seqs = self.seqs
        seq0 = seqs[0]
        for colname in map(itemgetter(0),columns)[7:-1]:
            if hasattr(seq0, colname):
                attr = set([])
                for seq in self.seqs:
                    if hasattr(seq, colname):
                        attr.add(getattr(seq,colname))

                attr = np.sort(list(attr))
                setattr(self, colname, attr)
                if 'days' in colname:
                    self.visit = np.arange(len(attr))
                    for seq in self.seqs:
                        if hasattr(seq, colname):
                            inds = (getattr(self, colname) == getattr(seq,colname))
                            seq.visit = self.visit[inds][0]


    def filter_seqs_with_attribute(self,attr):
        '''Filter only sequences that have that attribute.
        
        Note: this function is ugly because MultipleSeqAlignment does not
        support anything like the list.pop function, nor iterators for slicing
        but only pure slice objects.
        
        Result: I have to emulate the library itself...'''

        from Bio.Align import MultipleSeqAlignment
        f = lambda x: hasattr(x,attr)
        self.seqs =  MultipleSeqAlignment(filter(f,self.seqs), self.seqs._alphabet)

        # Reset the unique features
        self.add_features_from_seqs()
        

    def coallele_frequencies(self, i1, a1, i2, a2, visits=None):
        '''Monitor double sweeps of a1 at position i1 and a2 at i2'''
    
        if visits is None:
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


    def get_LD(self, visits=None, threshold=0.1):
        '''Calculate the LD (a single number) of certain visits.'''

        (dists, LDs) = self.get_LD_dist(visits, threshold)

        if not len(LDs):
            raise ValueError('No pairs found for calculating LD.')

        if np.isscalar(LDs[0]):
            return np.mean(LDs)

        LD_fin = 0.0
        n = 0
        for LD in LDs:
            LD_fin += sum(LD)
            n += len(LD)
        LD_fin /= n

        return LD_fin


    def get_LD_dist(self, visits=None, threshold=0.1):
        '''Calculate pairwise LD (Rebecca/Fabio-wise).
        
        - visits: if None (default), all timepoints are used. Otherwise only
        from those visits.
        - threshold: only use pairs of sites/alleles, whose allele freqs are
        between threshold and 1- threshold, and polarize with respect to those.
        This option avoids focusing on very rare stuff.
        '''

        if visits is None:
            visits = self.visit
        elif isinstance(visits,int):
            visits = [visits]
        
        dists = []
        LDs = []

        for i, v in enumerate(visits):
            seqs = self.seqs_from_visit(v)
            smat = np.array(seqs)

            # Get the allele freqs
            af = np.zeros((len(alpha),smat.shape[1]))
            for j, a in enumerate(alpha):
                af[j] = (smat == a).mean(axis=0)

            # Get only polymorphic sites
            ind_poly = (af.max(axis=0) < (1 - threshold))

            # Get only non-gap positions
            ind_nongaps = (af[alpha == '-'] < 0.1).flatten() 

            # Final index
            ind = (ind_poly & ind_nongaps)
            inde = ind.nonzero()[0]

            if len(inde) < 2:
                dists.append([])
                LDs.append([])
            else:

                allele_major_ind = af.argmax(axis=0)
                allele_major = alpha[allele_major_ind]

                # Polarized allele freqs of the minor variants
                af_pol = (1 - af.max(axis=0))

                # Calculate coallele freqs
                dist = []
                LD = []
                ratios = np.zeros(4)
                for j, indj in enumerate(inde):
                    i1 = (smat[:,indj] != allele_major[indj])

                    for k, indk in enumerate(inde[:j]):
                        i2 = (smat[:,indk] != allele_major[indk])

                        ratios[0] = (i1 & i2).mean() / (af_pol[indj] * af_pol[indk])
                        ratios[1] = (-i1 & i2).mean() / ((1 - af_pol[indj]) * af_pol[indk])
                        ratios[2] = (i1 & -i2).mean() / (af_pol[indj] * (1 - af_pol[indk]))
                        ratios[3] = (-i1 & -i2).mean() / ((1 - af_pol[indj]) * (1 - af_pol[indk]))

                        dist.append(np.abs(indj - indk))
                        LD.append(1 - np.min(ratios))
            dists.append(dist)
            LDs.append(LD)

        dists = map(np.array,dists)
        LDs = map(np.array,LDs)

        if len(visits) == 1:
            return (dists[0], LDs[0])
        else:
            return (dists, LDs)


    def extend_sequences(self, seqs=None, alpha=alpha):
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


    def _plot_PCA_2D(self, ax=None, show=False, colors=None):
        try:
            self.U
        except AttributeError:
            raise AttributeError('Please perform PCA before plotting.')
    
        # Time points
        visits = self.get_visits()
    
        # Plot different time points with different colors
        if colors is None:
            n = len(self.visit)
            colors = cm.jet([int(255.0 * i / n) for i in xrange(n)])
        
        if ax is None:
            fig = ppl.figure(figsize=(6,5))
            ax = fig.add_subplot(1,1,1)
    
        for i,v in enumerate(self.visit):
            ind = (visits == v)
            
            x = self.U[ind,0]
            y = self.U[ind,1]
        
            ax.scatter(x,y,color=colors[i], label='v '+str(v))
        
        ax.set_xlabel('PC1', fontsize=18)
        ax.set_ylabel('PC2', fontsize=18)
        ax.set_title('Patient '+self.name)
    
        if show:
            ppl.show()


    def _plot_PCA_3D(self, ax=None, show=False, colors=None):
        try:
            self.U
        except AttributeError:
            raise AttributeError('Please perform PCA before plotting.')

        try:
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError('3D extension to matplotlib not found.')
    
        # Time points
        visits = self.get_visits()
    
        # Plot different time points with different colors
        if colors is None:
            n = len(self.visit)
            colors = cm.jet([int(255.0 * i / n) for i in xrange(n)])
        
        if ax is None:
            fig = ppl.figure(figsize=(6,5))
            ax = fig.add_subplot(1,1,1, projection='3d')
    
        for i,v in enumerate(self.visit):
            ind = (visits == v)
            
            x = self.U[ind,0]
            y = self.U[ind,1]
            z = self.U[ind,2]
        
            ax.scatter(x,y,z, s=40, color=colors[i])
        
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
            text = l.get_text().split('.')[0].lstrip(' ')
            if text in names:
                v = self.seqs[names.index(text)].visit
                try:
                    ind = visits.index(v)
                    l.set_color(cols[ind])
                    l.set_text(text)
                except ValueError:
                    # Some seqs are weirdly labelled (!)
                    l.set_text(text+' BAD!')
        if show:
            ppl.show()

        return tree


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



# Functions
def parse_sequences_raw(filename):
    '''Parse sequences from raw LANL file and other stuff at the same time.'''
    patients = []
    names = []

    if not re.findall('\.dat$',str(filename)):
        filename = str(filename)+'.dat'

    # Track skipped columns
    skipped = np.zeros(len(columns), bool)


    # Get the sequences and add them to the patients
    with open(datadir+'/'+filename) as f:
        line = f.readline()     # Number of records
        line = f.readline()     # Header    
        for line in f:
            ls = line.split('\t')

            if len(ls) > 1:
                name = ls[1]
                if name not in names:
                    names.append(name)
                    p = patient(name)
                    p.ID = ls[2]
                    patients.append(p)
                else:
                    p = patients[names.index(name)]

                # Avoid double counting of seqs
                seqid = ls[3]
                if seqid not in p.get_ids():
                    rec = SeqRecord(Seq(ls[-1].strip('\n').upper(), ambiguous_dna),
                                    id=seqid,
                                    name=seqid,
                                    description=ls[4])
    
                    for i in xrange(7,len(ls)-1):
                        # This if is independent for each line, which allows
                        # incomplete records mixed with complete ones.
                        if ls[i]:
                            try:
                                setattr(rec,columns[i][0],columns[i][1](ls[i]))
                            except ValueError:
                                if not skipped[i]:
                                    from warnings import warn
                                    warn(columns[i][0]+' not formatted properly. Skipping.',
                                         RuntimeWarning)
                                    skipped[i] = True
        
                    p.seqs.append(rec)

    return patients


def parse_sequences(filename, reference=None, exclude=None,
                    patientclass=patient, **kwargs):
    '''Parse sequences and other stuff at the same time.
    
       Parameters:
    - filename is the number corresponding to this dataset among the ones
      downloaded
    - reference is None if self-alignment is used, is a string if some reference
      sequence is to be used 
    - patientclass is the class to use for patients. Useful only for subclassing
      parser_LANL.patient (e.g. parser_Bunnik2008.patient)
    '''
    patients = []
    names = []

    if not re.findall('\.dat$',str(filename)):
        filename = str(filename)+'.dat'

    # Create patients
    dirname = datadir+'/'+filename.strip('.dat')+'_fasta'
    names = glob.glob(dirname+'/*.fasta')
    # Get rid of dups
    names = list(set([os.path.basename(x).split('.')[0].split('_')[0] for x in names]))
    patients = [patientclass(name) for name in  names]

    # Read seqs
    for p in patients:
        # Shall we go for an alignment with reference?
        if reference is not None:
            p.file_sequences = dirname+'/to_reference/'+p.name+'_to_'+str(reference)+'_aligned.fasta'
            alignment = AlignIO.read(p.file_sequences,
                                     format='fasta',
                                     alphabet=ambiguous_dna)
            p.reference_name = str(reference)
            p.reference_seq = alignment[0]
            p.seqs = alignment[1:]
            if 'SHAPE' in reference:
                (p.SHAPE_p, p.SHAPE_r) = pSHAPE.assign_prob(str(p.reference_seq.seq),
                                                            include_reactivity=True)

            # Flag for excluding the map, which takes a while
            if 'SHAPEmap' in kwargs and kwargs['SHAPEmap']:
                p.SHAPEmap = pSHAPE.assign_map(str(p.reference_seq.seq))

        # Otherwise, do not use any reference sequence (self-aligned)
        else:
            # Try to use aligned sequences (and the AlignIO module)
            if os.path.isfile(dirname+'/'+p.name+'_aligned.fasta'):
                p.file_sequences = dirname+'/'+p.name+'_aligned.fasta'
                alignment = AlignIO.read(p.file_sequences,
                                         format='fasta',
                                         alphabet=ambiguous_dna)
                p.seqs = alignment
            else:
                p.file_sequences = dirname+'/'+p.name+'.fasta'
                iterator = SeqIO.parse(p.file_sequences,
                                       format='fasta',
                                       alphabet=ambiguous_dna)
                recstmp = []
                for seq_record in iterator:
                    recstmp.append(seq_record)
                p.seqs = recstmp

        # Jalview has a bad habit of adding lengths to the labels
        for rec in p.seqs:
            rec.id = rec.name = rec.name.split('/')[0]

    # Track skipped columns
    skipped = np.zeros(len(columns), bool)

    # Read other info from original file
    with open(datadir+'/'+filename) as f:
        line = f.readline()     # Number of records
        line = f.readline()     # Header
        for line in f:
            ls = line.split('\t')
            if len(ls) > 1:
                # Check whether or not the sequence is in any alignment, and in
                # which ones (it could be used several times after being chopped
                # into pieces)
                for p in patients:
                    tmp = (p.get_ids() == ls[3]).nonzero()[0]
                    if len(tmp):
                        s = p.seqs[tmp[0]]
                        for i in xrange(7, len(ls)-1):
                            if ls[i]:
                                try:
                                    setattr(s, columns[i][0], columns[i][1](ls[i]))
                                except ValueError:
                                    if not skipped[i]:
                                        from warnings import warn
                                        warn(columns[i][0]+' not formatted properly. Skipping.',
                                             RuntimeWarning)
                                        skipped[i] = True
    
    
    # Add file_tree if found
    for p in patients:
        if os.path.isfile(dirname+'/'+p.name+'_nj.tree'):
            p.file_tree = dirname+'/'+p.name+'_nj.tree'

    # Add path to private RNA secondary structure prediction
    for p in patients:
        RNAsecfile = dirname+'/initial_consensus/T_330.15/nuc_'+p.name+'.dbn'
        if os.path.isfile(RNAsecfile):
            p.RNAsecfile = RNAsecfile

    # Refine the patients before publishing
    for p in patients:
        p.add_features_from_seqs()

    # Filter out the excluded ones
    if exclude is not None:
        patients = filter(lambda x: x.name not in exclude, patients)

    return patients


def convert_to_fasta(filename):
    if re.findall('\.dat$',str(filename)):
        filename = str(filename)[:-4]
    else:
        filename = str(filename)
    
    dirname = datadir+'/'+filename+'_fasta'

    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    patients = parse_sequences_raw(filename)

    for p in patients:
        file_patient = dirname+'/'+p.name+'.fasta'
        SeqIO.write(p.seqs, file_patient, 'fasta')


def align_sequences(filename, translate=True, reference=None):
    '''Multiple Sequence Alignment using DNA or aa.
    
    Parameters:
    - filename is the number corresponding to this dataset among the ones
      downloaded
    - translate is True if amino acids are to be aligned and the aligment to be
      back-propagated instead of a direct DNA alignment
    - reference is None if self-alignment is used, is a string if some reference
      sequence is to be used
    '''

    from Bio.Align.Applications import MuscleCommandline
    # Find the directory
    if re.findall('\.dat$',str(filename)):
        filename = filename[:-4]
    else:
        filename = str(filename)
    dirname = datadir+'/'+filename+'_fasta'

    # Find the unaligned files
    names = glob.glob(dirname+'/*.fasta')
    names = [f.rstrip('.fasta').rstrip('_aligned') for f in names]
    names = list(set(names))

    # For each patient, do the alignment
    for name in names:
        print os.path.basename(name)

        # Check whether a reference is given
        if reference is None:
            fileDNAin = name+'.fasta'
            fileDNAout = name+'_aligned.fasta'
            fileaain = name+'_aa.fasta'
            fileaaout = name+'_aligned_aa.fasta'

        # if there is a reference, we need to create new files which include the
        # reference in a special folder
        else:
            bname = os.path.basename(name)
            sref = str(reference)
            reffile = dirname+'/reference/'+sref+'_trimmed.fasta'
            fileDNAin = dirname+'/to_reference/'+bname+'_to_'+sref+'.fasta'
            fileDNAout = dirname+'/to_reference/'+bname+'_to_'+sref+'_aligned.fasta'
            fileaain = dirname+'/to_reference/'+bname+'_to_'+sref+'_aa.fasta'
            fileaaout = dirname+'/to_reference/'+bname+'_to_'+sref+'_aligned_aa.fasta'

            # Create folder if not present
            if not os.path.isdir(dirname+'/to_reference'):
                os.mkdir(dirname+'/to_reference')

            # Copy stuff into the new folder
            import shutil
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


def get_ids(seqs):
    return np.array(map(attrgetter('id'), seqs))


def get_names(seqs):
    return np.array(map(attrgetter('name'), seqs))


def get_visits(seqs):
    return np.array(map(attrgetter('visit'), seqs),int)


def get_diversity(smat):
    '''Measure the diversity in a sample of sequences.'''
    smat= np.asarray(smat)
    N = smat.shape[0]
    div = []
    for i in xrange(N):
        for j in xrange(i):
            div.append((smat[i] != smat[j]).mean())
    return np.mean(div)


