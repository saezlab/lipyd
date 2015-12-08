#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
    # filtering features accross samples and fractions
    # to find those relevant ones 
    # belonging to LTP bound lipids
#

import os
import sys
import struct
try:
    import cPickle as pickle
except:
    import pickle

import warnings
import re
import timeit
import xlrd
from xlrd.biffh import XLRDError
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects.packages as rpackages
rvcd = rpackages.importr('vcdExtra')
rbase = rpackages.importr('base')
rutils = rpackages.importr('utils')
rstats = rpackages.importr('stats')
rococo = rpackages.importr('rococo')

import mass
import progress
import _curl
from common import *

import rlcompleter, readline
readline.parse_and_bind('tab:complete')

warnings.filterwarnings('error')
warnings.filterwarnings('default')

path_root = '/'
basedir = os.path.join(path_root, 'home', 'denes', 'Documents' , 'ltp')
ltpdirs = [os.path.join(basedir, 'share'),
    os.path.join(basedir, 'share', '2015_06_Popeye')]
samplesf = os.path.join(ltpdirs[0], 'control_sample.csv')
ppfracf = os.path.join(ltpdirs[0], 'fractions.csv')
ppsecdir = os.path.join(ltpdirs[0], 'SEC_profiles')
pptablef = os.path.join(basedir, 'proteins_by_fraction.csv')
swisslipids_url = 'http://www.swisslipids.org/php/export.php?action=get&file=lipids.csv'

metrics = [
        ('Kendall\'s tau', 'ktv', False),
        ('Spearman corr.', 'spv', False),
        ('Pearson corr.', 'pev', False),
        ('Euclidean dist.', 'euv', True),
        ('Robust corr.', 'rcv', False),
        ('Goodman-Kruskal\'s gamma', 'gkv', False),
        ('Difference', 'dfv', True)
    ]

class MolWeight():
    
    # 
    # Thanks for https://github.com/bsimas/molecular-weight/blob/master/chemweight.py
    #
    
    def __init__(self, formula = None, **kwargs):
        '''
            **kwargs: elements & counts, e.g. c = 6, h = 12, o = 6...
        '''
        if not hasattr(mass, 'massFirstIso'):
            mass.getMassFirstIso()
        self.mass = mass.massFirstIso
        self.reform = re.compile(r'([A-Za-z][a-z]*)([0-9]*)')
        if formula is None:
            formula = ''.join('%s%u'%(elem.upper(), num) for elem, num in kwargs.iteritems())
        self.formula = formula
        self.calc_weight()
    
    def __neg__(self):
        return -1 * self.weight
    
    def __add__(self, other):
        return float(other) + self.weight
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __iadd__(self, other):
        self.weight += float(other)
    
    def __sub__(self, other):
        return self.weight - float(other)
    
    def __rsub__(self, other):
        return float(other) - self.weight
    
    def __isub__(self, other):
        self.weight += float(other)
    
    def __truediv__(self, other):
        return self.weight / float(other)
    
    def __rtruediv__(self, other):
        return float(other) / self.weight
    
    def __itruediv__(self, other):
        self.weight /= float(other)
    
    def __mul__(self, other):
        return self.weight * float(other)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __imul__(self, other):
        self.weight *= float(other)
    
    def __float__(self):
        return self.weight
    
    def __eq__(self, other):
        return abs(self.weight - float(other)) <= 0.01
    
    def calc_weight(self):
        atoms = self.reform.findall(self.formula)
        w = 0.0
        for element, count in atoms:
            count = int(count or '1')
            w += self.mass[element] * count
        self.weight = w

class Mz():
    
    def __init__(self, mz, z = 1, sign = None, tolerance = 0.01):
        self.mz = mz
        self.z = z
        self.sign = sign
        self.tol = tolerance
    
    def __eq__(self, other):
        return self.z == other.z and \
            self.mz > other.mz - self.tol and \
            self.mz < other.mz + self.tol
    
    def __str__(self):
        return 'm/z = %f' % self.mz
    
    def adduct(self, m):
        return (self.mz * self.z + float(m)) / abs(self.z)
    
    def weight(self):
        return self.mz * self.z
    
    def remove_h(self):
        return self.adduct(-mass.proton)
    
    def remove_ac(self):
        m = MolWeight('H3C2O2')
        return self.adduct(-m - mass.electron)
    
    def remove_fo(self):
        m = MolWeight('HCO2')
        return self.adduct(-m - mass.electron)
    
    def remove_nh4(self):
        m = MolWeight('NH4')
        return self.adduct(-m + mass.electron)
    
    def remove_oh(self):
        m = MolWeight('OH')
        return self.adduct(-m - mass.electron)
    
    def add_h(self):
        return self.adduct(mass.proton)
    
    def add_oh(self):
        m = MolWeight('OH')
        return self.adduct(m + mass.electron)
    
    def add_fo(self):
        m = MolWeight('HCO2')
        return self.adduct(m + mass.electron)
    
    def add_ac(self):
        m = MolWeight('H3C2O2')
        return self.adduct(m + mass.electron)
    
    def add_nh4(self):
        m = MolWeight('NH4')
        return self.adduct(m - mass.electron)

#
# obtaining lipid data from SwissLipids
#

def get_id(dct, value):
    if value not in dct:
        i = 0 if len(dct) == 0 else max(dct.values()) + 1
        dct[value] = i
    return dct[value]

def get_lipidmaps(lipidmaps_url = 
    'http://www.lipidmaps.org/resources/downloads/LMSDFDownload28Jun15.tar.gz',
    fname = 'LMSDFDownload28Jun15/LMSDFDownload28Jun15FinalAll.sdf'):
    lipidmaps = []
    fields = ['lm_id', 'systematic_name', 'synonyms', 'category', 'main_class', 
        'exact_mass', 'formula', 'inchi_key', 'pubchem_sid', 'pubchem_cid']
    record = {}
    expect = False
    lmapsf = _curl.curl(lipidmaps_url, large = True, silent = False, files_needed = [fname])
    for l in lmapsf.values()[0]:
        if expect:
            record[expect] = l.strip() if expect != 'exact_mass' else float(l.strip())
            expect = False
        elif l[0] == '>':
            label = l.strip()[3:-1].lower()
            if label in fields:
                expect = label
        elif l[0] == '$':
            lipidmaps.append([record[label] if label in record else None \
                for label in fields])
            record = {}
    return lipidmaps

def lipidmaps_adducts(lipidmaps = None, pAdducts = None, nAdducts = None):
    _nAdducts = []
    _pAdducts = []
    lipidmaps = lipidmaps if lipidmaps is not None else get_lipidmaps()
    for l in lipidmaps:
        if l[5] is not None:
            _pAdducts.append([l[0], 'Species', l[1], l[6], '[M+H]+', Mz(l[5]).add_h()])
            _pAdducts.append([l[0], 'Species', l[1], l[6], '[M+NH4]+', Mz(l[5]).add_nh4()])
            _nAdducts.append([l[0], 'Species', l[1], l[6], '[M-H]-', Mz(l[5]).remove_h()])
            _nAdducts.append([l[0], 'Species', l[1], l[6], '[M+Fo]-', Mz(l[5]).add_fo()])
    _pAdducts = np.array(sorted(_pAdducts, key = lambda x: x[-1]), dtype = np.object)
    _nAdducts = np.array(sorted(_nAdducts, key = lambda x: x[-1]), dtype = np.object)
    if pAdducts is not None:
        _pAdducts = np.vstack((pAdducts, _pAdducts))
        del pAdducts
        _pAdducts = _pAdducts[_pAdducts[:,-1].argsort()]
    if nAdducts is not None:
        _nAdducts = np.vstack((nAdducts, _nAdducts))
        del nAdducts
        _nAdducts = _nAdducts[_nAdducts[:,-1].argsort()]
    return _pAdducts, _nAdducts

def lipidmaps_exact(lipidmaps = None, exacts = None):
    _exacts = []
    lipidmaps = lipidmaps if lipidmaps is not None else get_lipidmaps()
    for l in lipidmaps:
        if l[5] is not None:
            _exacts.append([l[0], 'Species', l[1], l[6], '', l[5]])
    _exactd = np.array(sorted(_exacts, key = lambda x: x[-1]), dtype = np.object)
    if exacts is not None:
        _exacts = np.vstack((exacts, _exacts))
        del exacts
        _exacts = _exacts[_exacts[:,-1].argsort()]
    return _exacts

def get_swisslipids(swisslipids_url, adducts = None, exact_mass = False):
    '''
    Downloads the total SwissLipids lipids dataset.
    Returns numpy array with all queried adducts and optionally exact masses,
    calculates the formate adduct m/z's on demand.
    Please note, for a number of lipids some adduct m/z's are missing from SL,
    while the exact masses are given for all the lipids. Also storing masses of 
    many adducts of the same lipid in an array results extensive memory usage
    (couple of hundreds MB).
    '''
    adduct_method = {
        '[M+OAc]-': 'add_ac',
        '[M+H]+': 'add_h',
        '[M-H]-': 'remove_h',
        '[M+Fo]-': 'add_fo',
        '[M+NH4]+': 'add_nh4'
    }
    if type(adducts) is list:
        adducts = set(adducts)
    readd = re.compile(r'.*(\[.*)')
    swl = _curl.curl(swisslipids_url, silent = False, compr = 'gz', large = True)
    swl.fileobj.seek(-4, 2)
    ucsize = struct.unpack('I', swl.fileobj.read(4))[0]
    swl.fileobj.seek(0)
    hdr = swl.readline().split('\t')
    positives = []
    negatives = []
    exact_masses = []
    prg = progress.Progress(ucsize, 'Processing SwissLipids', 101)
    for l in swl:
        prg.step(len(l))
        l = l.split('\t')
        if len(l) > 22:
            for i in xrange(13, 23):
                if len(l[i]) > 0:
                    add = '' if i == 13 else readd.findall(hdr[i])[0]
                    if not adducts or add in adducts or i == 13 and exact_mass or formiate and i == 13:
                        mz = to_float(l[i])
                        if i == 13:
                            if formiate:
                                mzfo = Mz(mz).add_fo()
                                foline = [l[0], l[1], l[2], l[10], '[M+Fo]-', mzfo]
                                negatives.append(foline)
                            if exact_mass:
                                exact_masses.append([l[0], l[1], l[2], l[10], add, mz])
                        elif hdr[i].endswith('+'):
                            positives.append([l[0], l[1], l[2], l[10], add, mz])
                        elif hdr[i].endswith('-'):
                            negatives.append([l[0], l[1], l[2], l[10], add, mz])
    prg.terminate()
    swl.close()
    positives = sorted(positives, key = lambda x: x[-1])
    positives = np.array(positives, dtype = np.object)
    negatives = sorted(negatives, key = lambda x: x[-1])
    negatives = np.array(negatives, dtype = np.object)
    return positives, negatives

def get_swisslipids_exact(swisslipids_url):
    '''
    Downloads the total SwissLipids lipids dataset.
    Returns numpy array with only exact mass values.
    '''
    readd = re.compile(r'.*(\[.*)')
    swl = _curl.curl(swisslipids_url, silent = False, compr = 'gz', large = True)
    swl.fileobj.seek(-4, 2)
    ucsize = struct.unpack('I', swl.fileobj.read(4))[0]
    swl.fileobj.seek(0)
    hdr = swl.readline().split('\t')
    exact_masses = []
    prg = progress.Progress(ucsize, 'Processing SwissLipids', 101)
    for l in swl:
        prg.step(len(l))
        l = l.split('\t')
        if len(l) > 22:
            mz = to_float(l[13])
            add = ''
            exact_masses.append([l[0], l[1], l[2], l[10], add, mz])
    prg.terminate()
    swl.close()
    exact_masses = sorted(exact_masses, key = lambda x: x[-1])
    exact_masses = np.array(exact_masses, dtype = np.object)
    return exact_masses

#
# reading the SEC absorption values and calculating the protein 
# profiles in the fractions
#

def protein_profiles(basedir, ppfracf):
    '''
    For each protein, for each fraction, calculated the mean of 
    absorptions of all the measurements belonging to one fraction.
    '''
    reltp = re.compile(r'.*[\s_-]([A-Za-z0-9]{3,})\.xls')
    result = {}
    fnames = os.listdir(basedir)
    with open(ppfracf, 'r') as f:
        frac = [(to_float(i[0]), to_float(i[1]), i[2]) for i in \
            [l.split(';') for l in f.read().split('\n')]]
    for fname in fnames:
        ltpname = reltp.findall(fname)[0]
        frac_abs = dict((i[2], []) for i in frac)
        gfrac = (i for i in frac)
        fr = gfrac.next()
        try:
            tbl = read_xls(os.path.join(basedir, fname))[3:]
        except xlrd.biffh.XLRDError:
            sys.stdout.write('Error reading XLS:\n\t%s\n'%os.path.join(basedir, fname))
        minabs = min(0.0, min(to_float(l[5]) for l in tbl))
        for l in tbl:
            # l[4]: volume (ml)
            ml = to_float(l[4])
            if ml < fr[0]:
                continue
            if ml >= fr[1]:
                try:
                    fr = gfrac.next()
                except StopIteration:
                    break
            # l[5] mAU UV3 215nm
            frac_abs[fr[2]].append(to_float(l[5]) - minabs)
        result[ltpname] = dict((fnum, np.mean(a)) for fnum, a in frac_abs.iteritems())
    return result

def zero_controls(samples, pprofs):
    fracs = ['a9', 'a10', 'a11', 'a12', 'b1']
    for ltpname, sample in samples.iteritems():
        for i, fr in enumerate(fracs):
            if sample[i + 1] == 0:
                pprofs[ltpname][fr] = 0.0

def write_pptable(pprofs, pptablef):
    '''
    Writes protein profiles in a table, so we don't need to read
    all the 60 XLS files every time.
    '''
    with open(pptablef, 'w') as f:
        f.write('\t%s%s' % ('\t'.join(sorted(pprofs.values()[0].keys(), 
            key = lambda x: (x[0], int(x[1:])))), '\n'))
        f.write('\n'.join('%s\t%s' % (ltp, 
                '\t'.join('%.20f'%d[fr] for fr in sorted(d.keys(), 
                    key = lambda x: (x[0], int(x[1:]))))) \
            for ltp, d in pprofs.iteritems()))

def read_pptable(pptablef):
    '''
    Reads protein profiles from table.
    '''
    with open(pptablef, 'r') as f:
        header = f.readline().strip().split('\t')
        return dict((ll[0], dict(zip(header, float_lst(ll[1:])))) \
            for ll in (l.split('\t') for l in f.read().split('\n')))

def read_xls(xls_file, sheet = '', csv_file = None, return_table = True):
    '''
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    '''
    try:
        book = xlrd.open_workbook(xls_file, on_demand = True)
        try:
            sheet = book.sheet_by_name(sheet)
        except XLRDError:
            sheet = book.sheet_by_index(0)
        table = [[str(c.value) for c in sheet.row(i)] for i in xrange(sheet.nrows)]
        if csv_file:
            with open(csv_file, 'w') as csv:
                csv.write('\n'.join(['\t'.join(r) for r in table]))
        if not return_table:
            table = None
        book.release_resources()
        return table
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
    sys.stdout.flush()

#
# reading in a small table for keeping track which fractions are samples 
# and which ones are controls
#

def read_samples(fname):
    '''
    Reads from file sample/control annotations 
    for each proteins.
    '''
    data = {}
    with open(fname, 'r') as f:
        null = f.readline()
        for l in f:
            l = l.split(',')
            data[l[0].replace('"', '')] = \
                np.array([to_int(x) if x != '' else None for x in l[1:]])
    return data

#
# scanning the directory tree and collecting the MS data csv files
#

def get_filenames(loc):
    '''
    Files are placed in 2 root directories, in specific subdirectories,
    positive and negative MS mode in separate files.
    This function collects all the file paths and returns them 
    in a dict of dicts. Keys are LTP names, and 'pos'/'neg'.
    '''
    redirname = re.compile(r'(^[0-9a-zA-Z]+)[ _]')
    fnames = {}
    ltpdd = os.listdir(loc)
    ltpdd = [i for i in ltpdd if os.path.isdir('%s/%s'%(loc,i))]
    if len(ltpdd) == 0:
        sys.stdout.write('\t:: Please mount the shared folder!\n')
        return fnames
    for ltpd in ltpdd:
        ltpname = redirname.findall(ltpd)
        pos = 'pos' if 'pos' in ltpd else 'neg' if 'neg' in ltpd else None
        if pos is not None and len(ltpname) > 0:
            ltpname = ltpname[0]
            fpath = [ltpd, 'features']
            for f in os.listdir(os.path.join(loc, *fpath)):
                if 'LABELFREE' in f:
                    fpath.append(f)
                    break
            for f in os.listdir(os.path.join(loc, *fpath)):
                if f.endswith('.csv'):
                    fpath.append(f)
                    break
            if ltpname not in fnames:
                fnames[ltpname] = {}
            if pos not in fnames[ltpname] or ltpd.endswith('update'):
                fnames[ltpname][pos] = os.path.join(loc, *fpath)
    return fnames

def read_file(fname):
    '''
    This function is deprecated.
    Returns list of lists, not numpy array.
    '''
    rehdr = re.compile(r'([_0-9a-zA-Z]+)[_\s]([/0-9a-zA-Z\s]+)')
    refra = re.compile(r'.*_[ab]{1,2}([0-9]{,2}).*')
    names = {
        'm/z': 0,
        'RT mean': 1,
        'Normalized Area': 2
    }
    data = []
    with open(fname, 'r') as f:
        hdr = f.readline().replace('Sample 5', '').split(',')[1:]
        cols = [tuple([i] + [x.strip() for x in list(rehdr.match(h).groups(0))]) \
            for i, h in enumerate(hdr[6:-7])]
        for c in cols:
            if c[0] % 3 != names[c[2]]:
                sys.stdout.write('erroneous column order: col %u with header %s'\
                    '\n\tin file %s\n' % (c[0], c[2], fname))
        for l in f:
            l = [i.strip() for i in l.split(',')][1:]
            common_cols = float_lst(l[:3]) + \
                [to_float(i.strip()) for i in l[3].split('-')] + \
                [int(l[4]), to_float(l[5])] + \
                [float_lst(l[-7].split(':'))] + \
                float_lst(l[-6:-4]) + \
                [float_lst(l[-4].split(':'))] + \
                [to_float(l[-3])]
            fractions = {}
            for c in cols:
                fnum = 0 if 'ctrl' in c[1] or 'ctlr' in c[1] \
                    else int(refra.match(c[1]).groups(0)[0])
                if fnum not in fractions:
                    fractions[fnum] = []
                fractions[fnum].append(to_float(l[c[0] + 6]))
            common_cols.append(fractions)
            data.append(common_cols)
    return data

def read_file_np(fname, read_vars = ['Normalized Area']):
    '''
    Reads one MS file, returns numpy masked array, 
    with void mask.
    Column order:
    quality, m/z, rt-min, rt-max, charge, control, a9, a10, a11, a12, b1
    '''
    # typos in the headers what need to be fixed
    typos = {
        'Sample 5': '',
        'Sample 6': '',
        '_ ': '_'
    }
    rehdr = re.compile(r'([_0-9a-zA-Z]+)[_\s]([/0-9a-zA-Z\s]+)')
    refra = re.compile(r'.*_[ab]{1,2}([0-9]{,2}).*')
    retyp = re.compile(r'(' + '|'.join(typos.keys()) + r')')
    # order of samples (fractions)
    sname = [0, 9, 10, 11, 12, 1]
    # variable names
    vname = {
        'm/z': 0,
        'RT mean': 1,
        'Normalized Area': 2
    }
    data = []
    # col nums for each variable, for each fraction
    scols = dict([(var, dict([(i, None) for i in sname])) for var in vname.keys()])
    with open(fname, 'r') as f:
        hdr = retyp.sub(lambda x: typos[x.group()], f.readline()).split(',')[1:]
        cols = [tuple([i] + [x.strip() for x in list(rehdr.match(h).groups(0))]) \
            for i, h in enumerate(hdr[6:-7])]
        # testing if order of columns is always 
        # m/z, RT, norm area; all files passed
        for c in cols:
            if c[0] % 3 != vname[c[2]]:
                sys.stdout.write('erroneous column order: col %u with header %s'\
                    '\n\tin file %s\n' % (c[0], c[2], fname))
            fnum = 0 if 'ctrl' in c[1] or 'ctlr' in c[1] \
                    else int(refra.match(c[1]).groups(0)[0])
            # assigning columns to variable->fraction slots
            # col offsets from 6
            scols[c[2]][fnum] = c
        for l in f:
            l = [i.strip() for i in l.split(',')][1:]
            vals = [to_float(l[0]), to_float(l[2])] + \
                [to_float(i.strip()) for i in l[3].split('-')] + \
                [to_float(l[4])]
            for var, scol in scols.iteritems():
                # number of columns depends on which variables we need
                if var in read_vars:
                    for s in sname:
                        if scol[s] is None:
                            vals.append(None)
                        else:
                            # col offset is 6 !
                            vals.append(to_float(l[scol[s][0] + 6]))
            data.append(vals)
    data.sort(key = lambda x: x[1])
    return np.ma.masked_array(data, dtype = 'float64')

def read_data(fnames, samples):
    '''
    Iterates through dict of dicts with file paths, reads
    each file, and makes array view and masks to make it easy
    to handle missing data, and access m/z values of samples 
    and controls, and other values.
    '''
    data = dict((ltp, {}) for ltp in fnames.keys())
    prg = progress.Progress(len(fnames) * 2, 'Reading files', 1)
    for ltp, pos_neg in fnames.iteritems():
        for p, fname in pos_neg.iteritems():
            prg.step()
            try:
                data[ltp][p] = {}
                # this returns a masked array:
                data[ltp][p]['raw'] = read_file_np(fname)
            except:
                sys.stdout.write('\nerror reading file: %s\n'%fname)
                sys.stdout.flush()
            # making a view with the intensities:
            data[ltp][p]['int'] = data[ltp][p]['raw'][:, 5:]
            # mask non measured:
            data[ltp][p]['mes'] = data[ltp][p]['int'].view()
            data[ltp][p]['mes'].mask = np.array([[x is None for x in samples[ltp]] * \
                    ((data[ltp][p]['int'].shape[1] - 5) / 6)] * \
                data[ltp][p]['int'].shape[0])
            data[ltp][p]['mes'] = np.ma.masked_invalid(data[ltp][p]['mes'])
            # controls:
            data[ltp][p]['ctr'] = data[ltp][p]['mes'][:,
                np.array([x == 0 for x in samples[ltp]] * \
                    (data[ltp][p]['mes'].shape[1] / 6))]
            # samples with lipids:
            data[ltp][p]['lip'] = data[ltp][p]['mes'][:,
                np.array([x == 1 for x in samples[ltp]] * \
                    (data[ltp][p]['mes'].shape[1] / 6))]
            # all samples except blank control:
            data[ltp][p]['smp'] = data[ltp][p]['mes'][:,
                np.array([False] + [x is not None \
                    for x in samples[ltp][1:]] * \
                    (data[ltp][p]['mes'].shape[1] / 6))]
    prg.terminate()
    return data

#
# Generic helper functions
#

def to_float(num):
    '''
    Extracts float from string, or returns None.
    '''
    renum = re.compile(r'([-]?[0-9]*[\.]?[0-9]+[eE]?[-]?[0-9]*)')
    num = renum.match(num.strip())
    if num:
        return float(num.groups(0)[0])
    else:
        return num

def to_int(num):
    '''
    Extracts int from string or returns None.
    '''
    renum = re.compile(r'([-]?[0-9]+[\.]?[0-9]*)')
    num = renum.match(num.strip())
    if num:
        return int(num.groups(0)[0])
    else:
        return num

def float_lst(l):
    '''
    Converts elements of a list to floats.
    '''
    return [to_float(x) for x in l]

#
# Filtering functions
#

def quality_filter(data, threshold = 0.2):
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            tbl['qly'] = np.array(tbl['raw'][:,0] >= threshold)

def charge_filter(data, charge = 1):
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            tbl['crg'] = np.array(tbl['raw'][:,4] == charge)

def area_filter(data, area = 10000.0):
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            tbl['are'] = np.nanmax(tbl['lip'], 1) >= area

def peaksize_filter(data, peak = 2.0):
    prg = progress.Progress(data.values()[0].values()[0]['raw'].shape[0] * 2, 
        'Peak size filter', 1)
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            prg.step()
            tbl['pks'] = np.nanmax(tbl['lip'], 1) / np.nanmax(tbl['ctr'], 1) >= peak
            # tbl['pks'] = lmax / cmax >= peak
            #tbl['pks'] = np.apply_along_axis(lambda x: , 1, tbl['lip'])
            #tbl['pks'] = np.empty([tbl['lip'].shape[0], 1])
            #for i, (c, s) in enumerate(zip(tbl['ctr'], tbl['lip'])):
                #tbl['pks'][i] = bool(sum(max(ss/cc >= peak for cc in c) for ss in s))
    prg.terminate()

def cprofile_filter(data, pprofs, samples):
    profile_filter(data, pprofs, samples, prfx = 'c')

def profile_filter(data, pprofs, samples, prfx = ''):
    frs = ['c0', 'a9', 'a10', 'a11', 'a12', 'b1']
    notf = []
    cols = 'smp' if prfx == 'c' else 'lip'
    prg = progress.Progress(len(data) * 2, 'Profile filter', 1, percent = False)
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            prg.step()
            if ltp.upper() not in pprofs:
                notf.append(ltp)
                continue
            # i != 0 : we always drop the blank control
            ppr = np.array([pprofs[ltp.upper()][frs[i]] \
                for i, fr in enumerate(samples[ltp]) if fr == 1 and i != 0])
            ppr = norm_profile(ppr).astype(np.float64)
            prr = stats.rankdata(ppr)
            pranks = sorted((i for i, x in enumerate(prr)), 
                key = lambda x: x, reverse = True)
            if tbl[cols][0,:].count() > 1:
                flatp = ppr[pranks[0]] - ppr[pranks[1]] < \
                    ppr[pranks[0]] * 0.1
            # 0 if have only one fraction in sample
            tbl['%sprf'%prfx] = np.apply_along_axis(
                lambda x: diff_profiles(ppr, norm_profile(x)), 
                    axis = 1, arr = tbl[cols])
            # all True if have only one fraction in sample
            tbl['%srpr'%prfx] = np.array([True] * tbl[cols].shape[0]) \
                if tbl[cols].shape[1] == 1 else \
                np.apply_along_axis(
                lambda x: comp_profiles(ppr, 
                    norm_profile(x), prr, flatp),
                    axis = 1, arr = tbl[cols])
    prg.terminate()
    sys.stdout.write('No protein profiles found for %s\n\n' % ', '.join(notf))
    sys.stdout.flush()

def diff_profiles(p1, p2):
    # profiles are numpy arrays
    # of equal length
    if len(p1) > 1:
        return np.nansum(np.abs(p1 - p2))
    else:
        return 1.0 * len(p2)

def _comp_profiles(p1, p2):
    p1r = stats.rankdata(p1)
    pranks = sorted((i for i, x in enumerate(p1r)), 
        key = lambda x: x, reverse = True)
    flatp = p1r[pranks[0]] - p1r[pranks[1]] < \
        p1r[pranks[0]] * 0.1
    return comp_profiles(p1, p2, p1r, flatp), 0.0

def comp_profiles(p1, p2, p1r = False, flatp = False):
    highest = False
    flat = False
    hollow = False
    p2r = stats.rankdata(p2)
    if type(p1r) is not np.ndarray:
        p1r = stats.rankdata(p1)
    highest = np.any(np.logical_and(p1r == 1, p2r == 1))
    p2ranks = sorted((i for i, x in enumerate(p2r)), 
        key = lambda x: x, reverse = True)
    if len(p1r) == 3:
        hollow = p1r[1] == 2 or p2r[1] == 2
    if not highest:
        flat = flatp and sum(p2) > 0.0 and \
            p2[p2ranks[0]] - p2[p2ranks[1]] < \
                p2[p2ranks][0] * 0.1
    return not hollow and (highest or flat)

def norm_profile(profile):
    return (profile - np.nanmin(profile.astype(np.float64))) / \
        (np.nanmax(profile) - np.nanmin(profile))

def norm_profiles(tbl):
    if type(tbl) == np.ma.core.MaskedArray:
        tbl = tbl.data
    return (tbl - tbl.min(axis = 1, keepdims = True)) / tbl.max(axis = 1, keepdims = True)

def ubiquity_filter_old(data, proximity = 0.02, only_valid = True):
    prg = progress.Progress(len(data)**2, 'Ubiquity filter', 1)
    ltps = sorted(data.keys())
    for pn in ['pos', 'neg']:
        for i, ltp1 in enumerate(ltps):
            for ltp2 in ltps[i+1:]:
                ltp1s = set([])
                ltp2s = set([])
                prg.step()
                if pn in data[ltp1] and pn in data[ltp2]:
                    ltp1t = data[ltp1][pn]
                    ltp2t = data[ltp2][pn]
                    if 'ubi' not in ltp1t:
                        ltp1t['ubi'] = np.zeros([ltp1t['raw'].shape[0], 1], dtype = np.int8)
                    if 'ubi' not in ltp2t:
                        ltp2t['ubi'] = np.zeros([ltp2t['raw'].shape[0], 1], dtype = np.int8)
                    for i1, mz1 in np.ndenumerate(ltp1t['raw'][:,1]):
                        i2u = ltp2t['raw'][:,1].searchsorted(mz1)
                        if i2u < ltp2t['raw'].shape[0] and ltp2t['raw'][i2u,1] - mz1 <= proximity:
                            if i1 not in ltp1s:
                                ltp1t['ubi'][i1] += 1
                                ltp1s.add(i1)
                            if i2u not in ltp2s:
                                ltp2t['ubi'][i2u] += 1
                                ltp2s.add(i2u)
                        if mz1 - ltp2t['raw'][i2u - 1,1] <= proximity:
                            if i1 not in ltp1s:
                                ltp1t['ubi'][i1] += 1
                                ltp1s.add(i1)
                            if i2u - 1 not in ltp2s:
                                ltp2t['ubi'][i2u - 1] += 1
                                ltp2s.add(i2u - 1)
    prg.terminate()

def ubiquity_filter(valids, proximity = 0.02, only_valid = True):
    prg = progress.Progress(len(valids)**2, 'Ubiquity filter', 1)
    ltps = sorted(valids.keys())
    for pn in ['pos', 'neg']:
        for i, ltp1 in enumerate(ltps):
            for ltp2 in ltps[i+1:]:
                ltp1s = set([])
                ltp2s = set([])
                prg.step()
                if pn in valids[ltp1] and pn in valids[ltp2]:
                    ltp1t = valids[ltp1][pn]
                    ltp2t = valids[ltp2][pn]
                    if 'ubi' not in ltp1t:
                        ltp1t['ubi'] = np.zeros(ltp1t['mz'].shape[0], dtype = np.int8)
                    if 'ubi' not in ltp2t:
                        ltp2t['ubi'] = np.zeros(ltp2t['mz'].shape[0], dtype = np.int8)
                    for i1, mz1 in np.ndenumerate(ltp1t['mz']):
                        if ltp1t['cpv'][i1] or not only_valid:
                            i2u = ltp2t['mz'].searchsorted(mz1)
                            u = 0
                            while True:
                                if i2u + u < ltp2t['mz'].shape[0] and \
                                    ltp2t['mz'][i2u + u] - mz1 <= proximity:
                                    if ltp2t['cpv'][i2u + u] or not only_valid:
                                        if i1 not in ltp1s:
                                            ltp1t['ubi'][i1] += 1
                                            ltp1s.add(i1)
                                        if i2u + u not in ltp2s:
                                            ltp2t['ubi'][i2u + u] += 1
                                            ltp2s.add(i2u + u)
                                    u += 1
                                else:
                                    break
                            l = 1
                            while True:
                                if i2u - l >= 0 and \
                                    mz1 - ltp2t['mz'][i2u - l] <= proximity:
                                    if ltp2t['cpv'][i2u - l] or not only_valid:
                                        if i1 not in ltp1s:
                                            ltp1t['ubi'][i1] += 1
                                            ltp1s.add(i1)
                                        if i2u - l not in ltp2s:
                                            ltp2t['ubi'][i2u - l] += 1
                                            ltp2s.add(i2u - l)
                                    l += 1
                                else:
                                    break
    prg.terminate()

def remove_filter(data, attr_name):
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            if attr_name in tbl:
                del tbl[attr_name]

def rprofile_filter(data):
    '''
    Dummy function for eval_filter()
    '''
    pass

class get_hits(object):
    
    def __init__(self, fun):
        self.fun = fun
    
    def __call__(self, data, **kwargs):
        result = empty_dict(data)
        for ltp in data.keys():
            print ltp
            for pn, tbl in data[ltp].iteritems():
                #try:
                result[ltp][pn] = self.fun(tbl, **kwargs)
                if result[ltp][pn] is None:
                    sys.stdout.write('No profile info for %s, %s, %s()\n' % \
                        (ltp, pn, self.fun.__name__))
                #except:
                #    pass
        return result

class combine_filters(object):
    
    def __init__(self, fun):
        self.fun = fun
    
    def __call__(self, data, **kwargs):
        for ltp in data.keys():
            for pn, tbl in data[ltp].iteritems():
                try:
                    self.fun(tbl, **kwargs)
                except:
                    pass

class count_hits(object):
    
    def __init__(self, fun):
        self.fun = fun
    
    def __call__(self, data, **kwargs):
        hits = empty_dict(data)
        phits = empty_dict(data)
        for ltp in data.keys():
            for pn, tbl in data[ltp].iteritems():
                try:
                    hits[ltp][pn] = self.fun(tbl, **kwargs)
                    phits[ltp][pn] = self.fun(tbl, **kwargs) / float(len(tbl[kwargs['name']])) * 100
                except:
                    pass
        return hits, phits

@get_hits
def val_ubi_prf_rpr_hits(tbl, ubiquity = 7, treshold = 0.15, tresholdB = 0.25, profile_best = False):
    '''
    Column order:
    quality, m/z, rt_min, rt_max, charge, control, a9, a10, a11, a12, b1, 
    profile_score, control_profile_score, rank_profile_boolean, ubiquity_score, ubiquity_score, 
    original_index
    '''
    prf = 'cprf' if tbl['lip'][1,:].count() == 1 else 'prf'
    _treshold = (treshold, tresholdB)
    if prf not in tbl:
        return None
    # dict of profile matching score values and their indices
    selected = np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
            tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['rpr']),
            tbl['ubi'] <= ubiquity), 
            np.logical_or(np.logical_not(np.isnan(tbl['prf'])), np.logical_not(np.isnan(tbl['cprf']))))
    prf_values = dict(zip(np.where(selected)[0], zip(tbl['prf'][selected], tbl['cprf'][selected])))
    if profile_best:
        # select the best scoring records
        prf_values = sorted(prf_values.items(), key = lambda x: x[1])
        try:
            score_treshold = prf_values[min(profile_best - 1, len(prf_values) - 1)][1]
        except IndexError:
            print prf_values
            return None
        # comments?
        indices = np.array(sorted(map(lambda x: x[0], \
            filter(lambda x: x[1] <= score_treshold, prf_values))))
    else:
        # or select the records with profile match less then or equal to treshold
        indices = np.array(sorted((i for i, v in prf_values.iteritems() if v <= _treshold)))
    return (tbl['raw'][indices,:], tbl['prf'][indices], tbl['cprf'][indices], 
        tbl['rpr'][indices], tbl['ubi'][indices], tbl['uby'][indices] if 'uby' in tbl else np.zeros(len(indices)), indices)

@get_hits
def pass_through(tbl, ubiquity = 7, treshold = 0.15, tresholdB = 0.25, profile_best = False):
    '''
    Column order:
    quality, m/z, rt_min, rt_max, charge, control, a9, a10, a11, a12, b1, 
    profile_score, control_profile_score, rank_profile_boolean, ubiquity_score, ubiquity_score
    '''
    prf = 'cprf' if tbl['lip'][1,:].count() == 1 else 'prf'
    _treshold = (treshold, tresholdB)
    if prf not in tbl:
        return None
    # dict of profile matching score values and their indices
    selected = np.logical_or(np.logical_not(np.isnan(tbl['prf'])), np.logical_not(np.isnan(tbl['cprf'])))
    prf_values = dict(zip(np.where(selected)[0], zip(tbl['prf'][selected], tbl['cprf'][selected])))
    if profile_best:
        # select the best scoring records
        prf_values = sorted(prf_values.items(), key = lambda x: x[1])
        try:
            score_treshold = prf_values[min(profile_best - 1, len(prf_values) - 1)][1]
        except IndexError:
            print prf_values
            return None
        # comments?
        indices = np.array(sorted(map(lambda x: x[0], \
            filter(lambda x: x[1] <= score_treshold, prf_values))))
    else:
        # or select the records with profile match less then or equal to treshold
        indices = np.array(sorted((i for i, v in prf_values.iteritems() if v <= _treshold)))
    print 'Selected: ', np.nansum(tbl['cprf'][indices])
    print 'Total: ', np.nansum(tbl['cprf'])
    return (tbl['raw'][indices,:], tbl['prf'][indices], tbl['cprf'][indices], 
        tbl['rpr'][indices], tbl['ubi'][indices], tbl['uby'][indices] if 'uby' in tbl else np.zeros(len(indices)))

@combine_filters
def validity_filter(tbl):
    tbl['vld'] = np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks'])

@combine_filters
def val_ubi_filter(tbl, ubiquity = 7):
    tbl['vub'] = np.logical_and(np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['ubi'][:,0] <= ubiquity)

@combine_filters
def val_prf_filter(tbl, treshold = 0.15):
    tbl['vpr'] = np.logical_and(np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['prf'] <= treshold)

@combine_filters
def val_rpr_filter(tbl):
    tbl['vrp'] = np.logical_and(np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['rpr'])

@combine_filters
def val_ubi_prf_filter(tbl, ubiquity = 7, treshold = 0.15):
    tbl['vup'] = np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['prf'] <= treshold), tbl['ubi'][:,0] <= ubiquity)

@combine_filters
def val_ubi_prf_rprf_filter(tbl, ubiquity = 7, treshold = 0.15):
    tbl['vur'] = np.logical_and(np.logical_and(np.logical_and(
        np.logical_and(np.logical_and(np.logical_and(
        tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['rpr']),
        tbl['prf'] <= treshold), tbl['ubi'][:,0] <= ubiquity)


def apply_filters(data, filters = ['quality', 'charge', 'area', 'peaksize'], param = {}):
    for f in filters:
        p = param[f] if f in param else {}
        fun = globals()['%s_filter' % f]
        fun(data, **p)

def empty_dict(data):
    return dict((ltp, {'pos': None, 'neg': None}) for ltp in data.keys())

def eval_filter(data, filtr, param = {}, runtime = True, repeat = 3, number = 10,
    hit = lambda x: x):
    t = None
    fun = globals()['%s_filter' % filtr]
    filters = {
            'quality': 'qly',
            'charge': 'crg',
            'peaksize': 'pks',
            'ubiquity': 'ubi',
            'area': 'are',
            'profile': 'prf',
            'rprofile': 'rpr',
            'cprofile': 'cprf',
            'validity': 'vld',
            'val_ubi': 'vub',
            'val_prf': 'vpr',
            'val_ubi_prf': 'vup',
            'val_ubi_prf_rprf': 'vur',
            'val_rpr': 'vrp'
        }
    name = filters[filtr]
    if runtime:
        t = min(timeit.Timer(lambda: fun(data, **param)).\
            repeat(repeat = repeat, number = number)) / float(number)
    @count_hits
    def counter(tbl, name = None):
        return np.nansum(hit(tbl[name]))
    hits, phits = counter(data, name = name)
    nums = np.array(list(j for i in  hits.values() for j in i.values()), dtype = np.float)
    pcts = np.array(list(j for i in  phits.values() for j in i.values()), dtype = np.float)
    return t, hits, phits, (np.nanmin(nums), np.nanmax(nums)), (np.nanmin(pcts), np.nanmax(pcts))

def combined_hits(data, profile = 0.15, ubiquity = 20, verbose = True):
    @count_hits
    def counter(tbl):
        return np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
                tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['prf'] <= profile), tbl['ubi'][:,0] <= ubiquity).sum()
    hits = counter(data)
    if verbose:
        sys.stdout.write('\n\tLTP\t\t+\t-\n\t%s\n'%('='*30))
        for ltp in sorted(hits.keys()):
            sys.stdout.write('\t%s\t'%ltp)
            for pn in ['pos', 'neg']:
                num = hits[ltp][pn]
                sys.stdout.write('\t%s' % (str(num) if num is not None else 'n/a'))
            sys.stdout.write('\n')
    return hits

def save_data(data, basedir, fname = 'features.pickle'):
    pickle.dump(data, open(os.path.join(basedir, fname), 'wb'))

def load_data(basedir, fname = 'features.pickle'):
    return pickle.load(open(os.path.join(basedir, fname), 'rb'))

def save(fnames, samples, pprofs, basedir, fname = 'save.pickle'):
    pickle.dump((fnames, samples, pprofs), open(os.path.join(basedir, fname), 'wb'))

def load(basedir, fname = 'save.pickle'):
    return pickle.load(open(os.path.join(basedir, fname), 'rb'))

def find_lipids(hits, pAdducts, nAdducts, levels = ['Species'], tolerance = 0.02):
    '''
    Column order:
    
    in:
    [0] quality, m/z, rt_min, rt_max, charge, control, a9, a10, a11, a12, b1
    [1] profile_score
    [2] control_profile_score
    [3] rank_profile_boolean
    [4] ubiquity_score
    [5] ubiquity_score
    
    out:
    ltp_name, m/z, 
    profile_score, control_profile_score, rank_profile_boolean, ubiquity_score, ubiquity_score,
    swisslipids_ac, level, lipid_name, lipid_formula, adduct, adduct_m/z
    '''
    # levels: 'Structural subspecies', 'Isomeric subspecies', 
    # 'Species', 'Molecular subspecies'
    levels = levels if type(levels) is set \
        else set(levels) if type(levels) is list \
        else set([levels])
    lipids = dict((ltp.upper(), {}) for ltp in hits.keys())
    for ltp, d in hits.iteritems():
        for pn, tbl in d.iteritems():
            adducts = pAdducts if pn == 'pos' else nAdducts
            result = []
            if tbl[0] is not None:
                for i in xrange(tbl[0].shape[0]):
                    swlipids = adduct_lookup(tbl[0][i,1], adducts, levels, tolerance)
                    if swlipids is not None:
                        for lip in swlipids:
                            result.append(np.concatenate(
                                # tbl[1] and tbl[2] are the profile and cprofile scores
                                (np.array([ltp.upper(), tbl[0][i,1], tbl[1][i], tbl[2][i], 
                                    tbl[3][i], tbl[4][i], tbl[5][i]], dtype = np.object), lip), axis = 0))
                lipids[ltp.upper()][pn] = np.vstack(sorted(result, key = lambda x: x[1]))
    return lipids

def find_lipids_exact(valids, exacts, levels = ['Species'], tolerance = 0.02):
    '''
    Looks up lipids by m/z among database entries in
    `exacts`, and stores the result in dict under key
    `lip`, where keys are the original indices (`i`).
    '''
    adducts = {
        'pos': {
            '[M+H]+': 'remove_h',
            '[M+NH4]+': 'remove_nh4'
        },
        'neg': {
            '[M-H]-': 'add_h',
            '[M+Fo]-': 'remove_fo'
        }
    }
    levels = levels if type(levels) is set \
        else set(levels) if type(levels) is list \
        else set([levels])
    for ltp, d in valids.iteritems():
        for pn, tbl in d.iteritems():
            tbl['lip'] = {}
            for i in xrange(tbl['mz'].shape[0]):
                tbl['lip'][tbl['i'][i]] = adduct_lookup_exact(tbl['mz'][i], exacts, levels, adducts[pn], tolerance)

def adduct_lookup(mz, adducts, levels, tolerance = 0.02):
    '''
    Looks up m/z values in the table containing the reference database
    already converted to a pool of all possible adducts.
    (Does not convert the m/z to other adducts.)
    '''
    result = []
    iu = adducts[:,-1].searchsorted(mz)
    if adducts.shape[0] > iu:
        u = 0
        while True:
            if adducts[iu + u,-1] - mz <= tolerance:
                if adducts[iu + u,1] in levels:
                    result.append(adducts[iu + u,:])
                u += 1
            else:
                break
    if iu > 0:
        l = 1
        while True:
            if iu - l >= 0 and mz - adducts[iu - l,-1] <= tolerance:
                if adducts[iu - l,1] in levels:
                    result.append(adducts[iu - l,:])
                l += 1
            else:
                break
    return None if len(result) == 0 else np.vstack(result)

def adduct_lookup_exact(mz, exacts, levels, adducts, tolerance = 0.02):
    '''
    Looks up m/z values in the table containing the reference database
    casting the m/z to specific adducts.
    '''
    result = []
    for addName, addFun in adducts.iteritems():
        addMz = getattr(Mz(mz), addFun)()
        iu = exacts[:,-1].searchsorted(addMz)
        if exacts.shape[0] > iu:
            u = 0
            while True:
                if exacts[iu + u,-1] - addMz <= tolerance:
                    if exacts[iu + u,1] in levels:
                        match = np.append(exacts[iu + u,:], addMz)
                        match[4] = addName
                        result.append(match)
                    u += 1
                else:
                    break
        if iu > 0:
            l = 1
            while True:
                if iu - l >= 0 and addMz - exacts[iu - l,-1] <= tolerance:
                    if exacts[iu - l,1] in levels:
                        match = np.append(exacts[iu - l,:], addMz)
                        match[4] = addName
                        result.append(match)
                    l += 1
                else:
                    break
    return None if len(result) == 0 else np.vstack(result)

def negative_positive(lipids, tolerance = 0.02, add_col = 12, mz_col = 1, swl_col = 8):
    '''
    Column order:
    
    in:
    ltp_name, m/z,
    profile_score, control_profile_score, rank_profile_boolean, ubiquity_score, ubiquity_score,
    original_index, 
    swisslipids_ac, level, lipid_name, lipid_formula, adduct, adduct_m/z
    
    out:
    [0] pos_m/z, pos_profile_score, pos_control_profile_score, pos_rank_profile_boolean, 
    [4] pos_ubiquity_score, pos_ubiquity_score, pos_original_index, pos_swisslipids_ac, pos_level, 
    [9] pos_lipid_name, pos_lipid_formula, pos_adduct, pos_adduct_m/z
    [13] neg_m/z, neg_profile_score, neg_control_profile_score, neg_rank_profile_boolean, 
    [17] neg_ubiquity_score, neg_ubiquity_score, neg_original_index, neg_swisslipids_ac, neg_level, 
    [22] neg_lipid_name, neg_lipid_formula, neg_adduct, neg_adduct_m/z
    '''
    result = dict((ltp.upper(), []) for ltp in lipids.keys())
    prg = progress.Progress(len(result), 'Matching positive & negative', 1, percent = False)
    for ltp, tbl in lipids.iteritems():
        prg.step()
        if 'neg' in tbl and 'pos' in tbl:
            for neg in tbl['neg']:
                adds = [neg[add_col]] if neg[add_col] != 'Unknown' else ['[M-H]-', '[M+Fo]-']
                for add in adds:
                    if add == '[M-H]-':
                        poshmz = Mz(Mz(neg[mz_col]).add_h()).add_h()
                        posnh3mz = Mz(Mz(neg[mz_col]).add_h()).add_nh4()
                    elif add == '[M+Fo]-':
                        poshmz = Mz(Mz(neg[mz_col]).remove_fo()).add_h()
                        posnh3mz = Mz(Mz(neg[mz_col]).remove_fo()).add_nh4()
                    else:
                        continue
                    iu = tbl['pos'][:,mz_col].searchsorted(poshmz)
                    u = 0
                    if iu < len(tbl['pos']):
                        while iu + u < len(tbl['pos']):
                            if tbl['pos'][iu + u,mz_col] - poshmz <= tolerance:
                                if tbl['pos'][iu + u,swl_col] == neg[swl_col] and \
                                    tbl['pos'][iu + u,add_col] == '[M+H]+' or \
                                    tbl['pos'][iu + u,add_col] == 'Unknown':
                                    result[ltp].append(np.concatenate(
                                        (tbl['pos'][iu + u,1:], neg[1:]), axis = 0))
                                u += 1
                            else:
                                break
                    if iu > 0:
                        l = 1
                        while iu >= l:
                            if poshmz - tbl['pos'][iu - l,mz_col] <= tolerance:
                                if tbl['pos'][iu - l,swl_col] == neg[swl_col] and \
                                    tbl['pos'][iu - l,add_col] == '[M+H]+' or \
                                    tbl['pos'][iu - l,add_col] == 'Unknown':
                                    result[ltp].append(np.concatenate(
                                        (tbl['pos'][iu - l,1:], neg[1:]), axis = 0))
                                l += 1
                            else:
                                break
                    iu = tbl['pos'][:,1].searchsorted(posnh3mz)
                    u = 0
                    if iu < len(tbl['pos']):
                        while iu + u < len(tbl['pos']):
                            if tbl['pos'][iu + u,mz_col] - posnh3mz <= tolerance:
                                if tbl['pos'][iu + u,swl_col] == neg[swl_col] and \
                                    tbl['pos'][iu + u,add_col] == '[M+NH4]+' or \
                                    tbl['pos'][iu + u,add_col] == 'Unknown':
                                    result[ltp].append(np.concatenate(
                                        (tbl['pos'][iu + u,1:], neg[1:]), axis = 0))
                                u += 1
                            else:
                                break
                    if iu > 0:
                        l = 1
                        while iu >= l:
                            if posnh3mz - tbl['pos'][iu - l,mz_col] <= tolerance:
                                if tbl['pos'][iu - l,swl_col] == neg[swl_col] and \
                                    tbl['pos'][iu - l,add_col] == '[M+NH4]+' or \
                                    tbl['pos'][iu - l,add_col] == 'Unknown':
                                    result[ltp].append(np.concatenate(
                                        (tbl['pos'][iu - l,1:], neg[1:]), axis = 0))
                                l += 1
                            else:
                                break
            if len(result[ltp]) > 0:
                result[ltp] = np.vstack(result[ltp])
    prg.terminate()
    return result

def best_matches(lipids, matches, minimum = 2, unknowns = None, unknown_matches = None):
    result = dict((ltp, {'pos': None, 'neg': None, 'both': None}) for ltp in lipids.keys())
    sort_lipids(lipids, by_mz = True)
    # returns a new array, or only lipids sorted
    alll = sort_lipids(lipids, by_mz = True, unknowns = unknowns)
    for ltp, tbl in matches.iteritems():
        utbl = unknown_matches[ltp] if unknown_matches is not None else None
        result[ltp]['both'] = _best_matches(tbl, minimum, utbl)
    if unknown is None:
        sort_lipids(lipids)
    for ltp, d in lipids.iteritems():
        for pn, tbl in d.iteritems():
            result[ltp][pn] = _best_lipids(tbl, minimum)
    return result

def _best_lipids(tbl, minimum = 2):
    key = lambda x: (x[4], x[3], x[5])
    return _best(tbl, key, minimum)

def _best(tbl, key, minimum):
    if minimum is None:
        return tbl
    k = None
    for i in xrange(tbl.shape[0]):
        if i > minimum - 1:
            thisKey = key(tbl[i,:])
            if k is not None and thisKey > k:
                break
            k = thisKey
    return tbl[:i,:]

def _best_matches(tbl, minimum = 2, utbl = None):
    if len(tbl) == 0 and utbl is None:
        return tbl
    if utbl is not None:
        if len(tbl) == 0:
            tbl = utbl
        else:
            tbl = np.vstack((tbl, utbl))
    # sorting: rprf, prf, cprf, ubi
    # key = lambda x: ((x[3] + x[16]), (x[1] + x[14]), (x[2] + x[15]), (x[4] + x[17]))
    key = lambda x: ((x[3] + x[16]), (x[2] + x[15]), (x[4] + x[17]))
    tbl = sorted(tbl, key = key)
    # print type(tbl), len(tbl)
    tbl = np.vstack(tbl)
    return _best(tbl, key, minimum, utbl)

def sort_lipids(lipids, by_mz = False, unknowns = None):
    in_place = unknowns is None
    result = dict((ltp, {'pos': None, 'neg': None}) for ltp in lipids.keys())
    if by_mz:
        key = lambda x: x[1]
    else:
        key = lambda x: (x[4], x[3], x[5])
    for ltp, d in lipids.iteritems():
        for pn, tbl in d.iteritems():
            _tbl = tbl if in_place else np.vstack((tbl, unknowns[ltp][pn]))
            if in_place:
                lipids[ltp][pn] = np.vstack(sorted(_tbl, key = key))
            else:
                result[ltp][pn] = np.vstack(sorted(_tbl, key = key))
    return lipids if in_place else result

#
# functions for MS2
#

def read_metabolite_lines(fname):
    '''
    From Toby Hodges.
    Reads metabolite fragments data.
    '''
    Metabolites = []
    with open(fname, 'r') as Handle:
        for Line in Handle.readlines():
            Line = Line.strip()
            MetabMass, MetabType, MetabCharge = Line.split('\t')
            Metabolites.append([to_float(MetabMass), MetabType, MetabCharge])
            if '+' not in MetabCharge and '-' not in MetabCharge and 'NL' not in MetabCharge:
                sys.stdout.write('WARNING: fragment %s has no valid charge information!\n' % \
                    metabolites[metabolite][1])
    return np.array(Metabolites, dtype = np.object)

def ms2_filenames(ltpdirs):
    '''
    Files are placed in 2 root directories, in specific subdirectories,
    positive and negative MS mode in separate files.
    This function collects all the file paths and returns them 
    in a dict of dicts. Keys are LTP names, and 'pos'/'neg'.
    '''
    redirname = re.compile(r'(^[0-9a-zA-Z]+)[ _](pos|neg)')
    refractio = re.compile(r'.*_([AB][0-9]{1,2}).*')
    fnames = {}
    for d in ltpdirs:
        ltpdd = os.listdir(d)
        ltpdd = [i for i in ltpdd if os.path.isdir(os.path.join(d, i))]
        if len(ltpdd) == 0:
            sys.stdout.write('\t:: Please mount the shared folder!\n')
            return fnames
        for ltpd in ltpdd:
            try:
                ltpname, pos = redirname.findall(ltpd)[0]
                ltpname = ltpname.upper()
                if ltpname not in fnames:
                    fnames[ltpname] = {}
                if pos not in fnames[ltpname] or ltpd.endswith('update'):
                    fnames[ltpname][pos] = {}
            except:
                # other dirs are irrelevant:
                continue
            fpath = [ltpd, 'Results']
            if os.path.isdir(os.path.join(d, *fpath)):
                for f in os.listdir(os.path.join(d, *fpath)):
                    if f.endswith('mgf') and 'buffer' not in f:
                        try:
                            fr = refractio.match(f).groups()[0]
                        except AttributeError:
                            sys.stdout.write('could not determine fraction number: %s'%f)
                        fnames[ltpname][pos][fr] = os.path.join(*([d] + fpath + [f]))
    return fnames

def ms2_map(ms2files):
    stRrtinseconds = 'RTINSECONDS'
    stRtitle = 'TITLE'
    stRbe = 'BE'
    stRen = 'EN'
    stRch = 'CH'
    stRpepmass = 'PEPMASS'
    stRpos = 'pos'
    stRneg = 'neg'
    stRcharge = 'CHARGE'
    stRempty = ''
    reln = re.compile(r'^([A-Z]+).*=([\d\.]+)[\s]?([\d\.]*)["]?$')
    redgt = re.compile(r'[AB]([\d]+)$')
    result = dict((ltp.upper(), {'pos': None, 'neg': None, 
            'ms2files': {'pos': {}, 'neg': {}}}) \
        for ltp, d in ms2files.iteritems())
    prg = progress.Progress(len(ms2files) * 2, 'Indexing MS2 data', 1, percent = False)
    for ltp, d in ms2files.iteritems():
        pFeatures = []
        nFeatures = []
        ultp = ltp.upper()
        for pn, files in d.iteritems():
            prg.step()
            features = pFeatures if pn == stRpos else nFeatures
            #fractions = set([])
            for fr, fl in files.iteritems():
                fr = int(redgt.match(fr).groups()[0])
                #fractions.add(fr)
                offset = 0
                with open(fl, 'r', 8192) as f:
                    for l in f:
                        if not l[0].isdigit() and not l[:2] == stRch and \
                            not l[:2] == stRbe and not l[:2] == stRen:
                            try:
                                m = reln.match(l).groups()
                            except:
                                print fl, l
                                continue
                            if m[0] == stRtitle:
                                scan = float(m[1])
                            if m[0] == stRrtinseconds:
                                rtime = float(m[1])
                            if m[0] == stRpepmass:
                                pepmass = float(m[1])
                                intensity = 0.0 if m[2] == stRempty else float(m[2])
                                features.append([pepmass, intensity,
                                    rtime, scan, offset + len(l), fr])
                                scan = None
                                rtime = None
                                intensity = None
                                pepmass = None
                        offset += len(l)
                result[ultp]['ms2files'][pn][int(fr)] = fl
        pFeatures = np.array(sorted(pFeatures, key = lambda x: x[0]), dtype = np.float64)
        nFeatures = np.array(sorted(nFeatures, key = lambda x: x[0]), dtype = np.float64)
        result[ultp]['pos'] = pFeatures
        result[ultp]['neg'] = pFeatures
    prg.terminate()
    return result

def ms2_main(valids, samples, ms2map, pFragments, nFragments, 
    tolerance = 0.01, verbose = False):
    '''
    For all LTPs and modes obtains the MS2 data from original files.
    '''
    # ms2map columns: pepmass, intensity, rtime, scan, offset, fraction
    prg = progress.Progress(len(valids) * 2, 'Looking up MS2 fragments', 1, 
        percent = False)
    for ltp, d in valids.iteritems():
        for pn, tbl in d.iteritems():
            prg.step()
            fragments = pFragments if pn == 'pos' else nFragments
            ms2matches = ms2_match(tbl['mz'], tbl['rt'], tbl['i'], ms2map, 
                ltp, pn, tolerance)
            # this already returns final result, from protein containing
            # fractions and with the relevant retention times
            tbl['ms2'] = ms2_lookup(ms2map[ltp][pn], ms2matches, samples[ltp], 
                ms2map[ltp]['ms2files'][pn], fragments)
            tbl['ms2r'] = ms2_result(tbl['ms2'])
            if verbose:
                print '\n'
                print 'number of positive mzs:', len(posMzs)
                print 'negative matching:', len(pos_matches)
                print 'number of negative mzs:', len(negMzs)
                print 'negative matching:', len(neg_matches)
    prg.terminate()

def ms2_result(ms2matches):
    '''
    Extracts the most relevant information from the MS2
    result arrays, throwing away the rest.
    '''
    result = dict((oi, []) for oi in ms2matches.keys())
    for oi, ms2s in ms2matches.iteritems():
        # columns: MS2 fragment name, MS2 adduct name, MS2 m/z, MS2 intensity
        result[oi].append([ms2s[7], ms2s[8], ms2s[1], ms2s[2]])
    for oi, ms2s in result.iteritems():
        result[oi] = np.vstack(ms2s)
    return result

def ms2_match(ms1Mzs, ms1Rts, ms1is, ms2map, ltp, pos, tolerance = 0.01):
    '''
    Looks up matching pepmasses for a list of MS1 m/z's.
    '''
    matches = []
    ms2tbl = ms2map[ltp][pos]
    # iterating over MS1 m/z, MS1 original index, and retention time
    for ms1Mz, ms1i, rt in zip(ms1Mzs, ms1is, ms1Rts):
        # if error comes here, probably MS2 files are missing
        try:
            iu = ms2tbl[:,0].searchsorted(ms1Mz)
        except IndexError:
            sys.stdout.write('\nMissing MS2 files for %s-%s?\n' % (ltp, pos))
        u = 0
        if iu < ms2tbl.shape[0]:
            while iu + u < ms2tbl.shape[0]:
                if ms2tbl[iu + u,0] - ms1Mz <= tolerance:
                    # checking retention time
                    if ms2tbl[iu + u, 2] >= rt[0] and ms2tbl[iu + u, 2] <= rt[1]:
                        matches.append((ms1Mz, iu + u, ms1i))
                        u += 1
                else:
                    break
        l = 1
        if iu > 0:
            while iu >= l:
                if ms1Mz - ms2tbl[iu - l,0] <= tolerance:
                    # checking retention time
                    if ms2tbl[iu - l, 2] >= rt[0] and ms2tbl[iu - l, 2] <= rt[1]:
                        matches.append((ms1Mz, iu - l, ms1i))
                        l += 1
                else:
                    break
    return sorted(uniqList(matches), key = lambda x: x[0])

def ms2_lookup(ms2map, ms1matches, samples, ms2files, fragments):
    '''
    For the matching MS2 m/z's given, reads and identifies
    the list of fragments.
    '''
    # indices of fraction numbers
    sample_i = {
        9: 1, 10: 2, 11: 3, 12: 4, 1: 5
    }
    stRcharge = 'CHARGE'
    # opening all files to have file pointers ready
    files = dict((fr, open(fname, 'r')) for fr, fname in ms2files.iteritems())
    # Initializing dict with original indices
    ms2matches = dict((m[2], []) for m in ms1matches)
    # iterating over MS1 m/z, MS2 table index, MS1 original index
    for ms1mz, ms2i, ms1oi in ms1matches:
        # ms2map columns: pepmass, intensity, rtime, scan, offset, fraction
        ms2item = ms2map[ms2i,:]
        fr = int(ms2item[5])
        # only samples with the LTP
        if samples[sample_i[fr]] == 1 and fr in files:
            f = files[fr]
            # jumping to offset
            f.seek(ms2item[4])
            # zero means no clue about charge
            charge = 0
            for l in f:
                if l[:6] == stRcharge:
                    # one chance to obtain the charge
                    charge = int(l[-1])
                    continue
                if not l[0].isdigit():
                    # finish at next section
                    break
                else:
                    # reading fragment masses
                    mass, intensity = l.strip().split()
                    mass = float(mass)
                    intensity = float(intensity)
                    # matching fragment --- direct
                    ms2hit1 = ms2_identify(mass, fragments, compl = False)
                    # matching fragment --- inverted
                    ms2hit2 = ms2_identify(ms2item[0] - mass, fragments, compl = True)
                    # matched fragment --- direct
                    # columns (14): 
                    # MS1 m/z, MS2 fragment m/z, MS2 fragment intensity, 
                    # MS2 table index, direct/inverted match, fraction number, 
                    # MS2 fragment m/z in annotation, MS2 fragment name,
                    # MS2 adduct name, 
                    # MS1 pepmass, MS1 intensity, rtime, MS2 scan, 
                    # MS2 file offset, fraction number
                    if ms2hit1 is not None:
                        ms2matches[oi].append(np.concatenate((np.array([ms1mz, mass, 
                            intensity, ms2i, 0, fr]), ms2hit1, ms2item)))
                    # matched fragment --- inverted
                    if ms2hit2 is not None:
                        ms2matches[oi].append(np.concatenate((np.array([ms1mz, mass, 
                            intensity, ms2i, 1, fr]), ms2hit2, ms2item)))
                    # no fragment matched --- unknown fragment
                    if ms2hit1 is None and ms2hit2 is None:
                        ms2matches[oi].append(np.concatenate((np.array([ms1mz, mass, 
                            intensity, ms2i, 0, fr]), 
                            np.array([None, 'unknown', 'unknown']), 
                            ms2item)))
    # removing file pointers
    for f in files.values():
        f.close()
    for oi, ms2match in ms2matches.iteritems():
        ms2matches[oi] = np.vstack(ms2match)
    return ms2matches

def ms2_identify(mass, fragments, compl, tolerance = 0.01):
    #print 'ms2_indentify() starts'
    iu = fragments[:,0].searchsorted(mass)
    if iu < len(fragments):
        if compl and 'NL' in fragments[iu,2] or \
            not compl and ('+' in fragments[iu,2] or '-' in fragments[iu,2]):
            if fragments[iu,0] - mass <= tolerance:
                #print 'ms2_identify() returns result'
                return fragments[iu,:]
            elif iu > 0 and mass - fragments[iu - 1,0] <= tolerance:
                #print 'ms2_identify() returns result'
                return fragments[iu - 1,:]
    #print 'ms2_identify() returns None'
    return None

def ms2_collect(ms2matches, ms1mz, unknown = False):
    #print 'ms2_collect() starts'
    result = []
    fragments = ms2matches[ms2matches[:,0] == ms1mz,:]
    for frag in uniqList(fragments[:,7]):
        if frag != 'Unknown':
            thisFrag = fragments[fragments[:,7] == frag,:]
            maxInt = thisFrag[:,2].max()
            thisFragMass = thisFrag[0,1]
            result.append((frag, maxInt, thisFragMass))
    if unknown:
        unknowns = fragments[fragments[:,7] == 'Unknown',:]
        result += [('Unknown', l[2], l[1]) for l in unknowns]
    #print 'ms2_collect() ends'
    return result

#
# Pipeline elements
#

def write_out(matches, fname):
    '''
    In:
    [0] pos_m/z, pos_profile_score, pos_control_profile_score, pos_rank_profile_boolean, 
    [4] pos_ubiquity_score, pos_ubiquity_score, pos_original_index, pos_swisslipids_ac, pos_level, 
    [9] pos_lipid_name, pos_lipid_formula, pos_adduct, pos_adduct_m/z
    [13] neg_m/z, neg_profile_score, neg_control_profile_score, neg_rank_profile_boolean, 
    [17] neg_ubiquity_score, neg_ubiquity_score, neg_original_index, neg_swisslipids_ac, neg_level, 
    [22] neg_lipid_name, neg_lipid_formula, neg_adduct, neg_adduct_m/z
    '''
    with open(fname, 'w') as f:
        hdr = ['LTP',
            'Positive_m/z', 'Positive_m/z_in_SwissLipids', 'Positive_adduct',
            'Negative_m/z', 'Negative_m/z_in_SwissLipids', 'Negative_adduct',
            'SwissLipids_AC', 'SwissLipids_formula', 'SwissLipids_name']
        f.write('\t'.join(hdr) + '\n')
        for ltp, tbl in matches.iteritems():
            for l in tbl:
                f.write('\t'.join([ltp, '%08f'%l[0], '%08f'%l[12], l[11], '%08f'%l[13], 
                    '%08f'%l[25], l[24], l[7], l[10], str(l[9])]) + '\n')

def counts_redundancy_table(lipids, unknowns):
    with open('unique_features_counts.csv', 'w') as f:
        f.write('\t'.join(['LTP-mode', 'unknown_features', 'lipid_matching_features', 'lipids']) + '\n')
        f.write('\n'.join('%s\t%u\t%u\t%u'%i for i in zip(['%s-%s'%(b,a) for b in unknowns.keys() for a in unknowns[b].keys()],
            [len(uniqList(list(a[:,7]))) for b in lipids.values() for a in b.values()], 
            [len(uniqList(list(a[:,7]))) for b in unknowns.values() for a in b.values()],
            [len(list(a[:,7])) for b in lipids.values() for a in b.values()])))

def samples_with_controls(samples):
    return dict((k, np.array([1 if x is not None else None for x in v])) \
        for k, v in samples.iteritems())

def upper_samples(samples):
    '''
    Returns the dict of samples with
    uppercase LTP names.
    '''
    return dict((l.upper(), s) for l, s in samples.iteritems())

def init_from_scratch(basedir, ltpdirs, pptablef, samplesf):
    '''
    Does all the initial preprocessing.
    Saves intermediate data, so it can be loaded faster 
    for next sessions.
    '''
    fnames = dict(get_filenames(ltpdirs[0]).items() + get_filenames(ltpdirs[1]).items())
    samples = read_samples(samplesf)
    # at first run, after reading from saved textfile
    pprofs = protein_profiles(ppsecdir, ppfracf)
    write_pptable(pprofs, pptablef)
    del fnames['ctrl']
    save(fnames, samples, pprofs, basedir)
    csamples = samples_with_controls(samples)
    samples_upper = upper_samples(samples)
    data = read_data(fnames, samples)
    save_data(data, basedir)
    return data, fnames, samples, csamples, samples_upper, pprofs

def init_reinit(basedir):
    '''
    Initializing from preprocessed and dumped data.
    Pickle file has a 2.0GB size.
    '''
    fnames, samples, pprofs = load(basedir)
    data = load_data(basedir)
    csamples = samples_with_controls(samples)
    samples_upper = upper_samples(samples)
    return data, fnames, samples, csamples, samples_upper, pprofs

def stage0_filters(data, pprofs, samples, csamples):
    '''
    Deprecated with new data structure.
    '''
    apply_filters(data)
    validity_filter(data)
    profile_filter(data, pprofs, samples)
    profile_filter(data, pprofs, csamples, prfx = 'c')
    ubiquity_filter(data)
    val_ubi_filter(data, ubiquity = ubiquity_treshold)
    rprofile_filter(data)
    val_prf_filter(data, treshold = profile_treshold)
    val_rpr_filter(data)
    val_ubi_prf_filter(data, treshold = profile_treshold, ubiquity = ubiquity_treshold)
    val_ubi_prf_rprf_filter(data, treshold = profile_treshold, ubiquity = ubiquity_treshold)

def basic_filters(data, pprofs, samples, csamples, 
    profile_treshold = 0.25, ubiquity_treshold = 7):
    '''
    Deprecated with new data structure.
    '''
    apply_filters(data)
    validity_filter(data)
    profile_filter(data, pprofs, samples)
    profile_filter(data, pprofs, csamples, prfx = 'c')
    ubiquity_filter(data)
    val_ubi_filter(data, ubiquity = ubiquity_treshold)
    rprofile_filter(data)
    val_prf_filter(data, treshold = profile_treshold)
    val_rpr_filter(data)
    val_ubi_prf_filter(data, treshold = profile_treshold, ubiquity = ubiquity_treshold)
    val_ubi_prf_rprf_filter(data, treshold = profile_treshold, ubiquity = ubiquity_treshold)

def basic_filters_with_evaluation(data, pprofs, samples, 
    profile_treshold = 0.25, ubiquity_treshold = 7):
    '''
    Deprecated with new data structure.
    '''
    filtr_results = {}
    for f in ['quality', 'charge', 'area', 'peaksize', 'validity']:
        filtr_results[f] = eval_filter(data, f)
    filtr_results['ubiquity'] = eval_filter(data, 'ubiquity', 
    param = {'only_valid': False},
    runtime = True, repeat = 1, number = 1, 
    hit = lambda x: x < ubiquity_treshold)
    filtr_results['val_ubi'] = eval_filter(data, 'val_ubi',
        param = {'ubiquity': ubiquity_treshold})
    filtr_results['profile025'] = eval_filter(data, 'profile', 
        param = {'pprofs': pprofs, 'samples': samples, 'prfx': ''},
        runtime = True, repeat = 1, number = 1, 
        hit = lambda x: x <= profile_treshold)
    filtr_results['cprofile025'] = eval_filter(data, 'cprofile', 
        param = {'pprofs': pprofs, 'samples': csamples},
        runtime = True, repeat = 1, number = 1, 
        hit = lambda x: x <= profile_treshold)
    filtr_results['rprofile'] = eval_filter(data, 'rprofile', runtime = True)
    filtr_results['val_prf'] = eval_filter(data, 'val_prf', 
        param = {'treshold': profile_treshold})
    filtr_results['val_rpr'] = eval_filter(data, 'val_rpr')
    filtr_results['val_ubi_prf'] = eval_filter(data, 'val_ubi_prf',
        param = {'treshold': profile_treshold, 'ubiquity': ubiquity_treshold})
    filtr_results['val_ubi_prf_rprf'] = eval_filter(data, 'val_ubi_prf_rprf',
        param = {'treshold': profile_treshold, 'ubiquity': ubiquity_treshold})
    return filtr_results

def positive_negative_runtime():
    return timeit.timeit('negative_positive(lipids)', 
        setup = 'from __main__ import negative_positive, lipids', number = 1)

def valid_features(data):
    '''
    Creates new dict of arrays with only valid features.
    Keys:
        'fe': features
        'mz': m/z values
        'i': original index
    '''
    apply_filters(data)
    validity_filter(data)
    valids = dict((ltp.upper(), {'pos': {}, 'neg': {}}) for ltp in data.keys())
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            valids[ltp.upper()][pn]['fe'] = np.array(tbl['smp'][tbl['vld']])
            valids[ltp.upper()][pn]['mz'] = np.array(tbl['raw'][tbl['vld'], 1])
            valids[ltp.upper()][pn]['rt'] = np.array(tbl['raw'][tbl['vld'], 2:4])
            valids[ltp.upper()][pn]['i'] = np.where(tbl['vld'])[0]
    return valids

def norm_all(valids):
    '''
    Creates table with all the profiles normalized
    Keys:
        'no': normalized profiles
    '''
    for ltp, d in valids.iteritems():
        for pn, tbl in d.iteritems():
            tbl['no'] = norm_profiles(tbl['fe'])

def profiles_corr(valids, pprofs, samples, metric, prfx):
    '''
    Calculates spearman's rank correlation
    between each feature and the protein profile.
    '''
    frs = ['c0', 'a9', 'a10', 'a11', 'a12', 'b1']
    for ltp, d in valids.iteritems():
        ppr = np.array([pprofs[ltp.upper()][frs[i]] \
            for i, fr in enumerate(samples[ltp]) if fr is not None and i != 0])
        ppr = norm_profile(ppr).astype(np.float64)
        for pn, tbl in d.iteritems():
            tbl['%sv'%prfx] = np.zeros((tbl['no'].shape[0],), dtype = np.float64)
            tbl['%sp'%prfx] = np.zeros((tbl['no'].shape[0],), dtype = np.float64)
            for i, fe in enumerate(tbl['no']):
                vp = metric(fe, ppr)
                tbl['%sv'%prfx][i] = vp[0]
                tbl['%sp'%prfx][i] = vp[1]

def gkgamma(x, y):
    '''
    Calls Goodman-Kruskal's gamma from vcdExtra
    R package.
    '''
    gkg = rvcd.GKgamma(rbase.matrix(rbase.c(*(list(x) + list(y))), nrow = 2))
    # gamma, 0.0, C, D, sigma, CIlevel, CI
    return tuple([gkg[0][0], 0.0] + [i[0] for i in gkg[1:]])

def roco(x, y):
    '''
    Calls R function robust correlation coefficient
    with test from rococo R package.
    '''
    x = rbase.c(*x)
    y = rbase.c(*y)
    rt = rococo.rococo_test( \
        rbase.as_vector(rstats.na_omit(x.rx( \
            rbase.c(*list(np.where(map(lambda i: not i, list(rbase.is_na(y))))[0] + 1))))), \
        rbase.as_vector(rstats.na_omit(y.rx( \
            rbase.c(*list(np.where(map(lambda i: not i, list(rbase.is_na(x))))[0] + 1))))))
    return rt.slots['sample.gamma'][0], rt.slots['p.value'][0]

def _diff_profiles(x, y):
    '''
    Wrapper for diff_profiles() to return a tuple.
    '''
    return diff_profiles(x, y), 0.0

def euclidean_dist(x, y):
    _x = x[np.where(~np.logical_or(np.isnan(x), np.isnan(y)))]
    _y = y[np.where(~np.logical_or(np.isnan(x), np.isnan(y)))]
    return (sp.spatial.distance.euclidean(_x, _y), 0.0) \
        if len(_x) > 0 else (np.inf, 0.0)

def profiles_corrs(valids, pprofs, samples):
    '''
    Calculates an array of similarity/correlation
    metrics between each MS intensity profile and 
    the corresponding protein concentration profile.
    '''
    metrics = [
        (stats.spearmanr, 'sp'),
        (stats.kendalltau, 'kt'),
        (_diff_profiles, 'df'),
        (gkgamma, 'gk'),
        (roco, 'rc'),
        (euclidean_dist, 'eu'),
        (stats.pearsonr, 'pe'),
        (_comp_profiles, 'cp')
    ]
    for metric, prfx in metrics:
        sys.stdout.write('Calculating %s\n' % metric.__name__)
        profiles_corr(valids, pprofs, samples, metric, prfx)

def distance_matrix(valids, metrics = ['eu']):
    _metrics = {
        'eu': ('Euclidean distance', euclidean_dist)
    }
    for m in metrics:
        sys.stdout.write('Calculating %s\n'%_metrics[m][0])
        for ltp, d in valids.iteritems():
            for pn, tbl in d.iteritems():
                fnum = tbl['no'].shape[0]
                tbl['%sd'%m] = np.empty((fnum, fnum))
                # to keep track the order accross sorting
                tbl['%so'%m] = np.copy(tbl['i'])
                for i in xrange(fnum - 1):
                    for j in xrange(i, fnum):
                        tbl['%sd'%m][i,j] = \
                            _metrics[m][1](tbl['no'][i,:], tbl['no'][j,:])

def kmeans(valids, pprofs, samples):
    frs = ['c0', 'a9', 'a10', 'a11', 'a12', 'b1']
    prg = progress.Progress(len(valids) * 2, 'Calculating k-means', 1, percent = False)
    for ltp, d in valids.iteritems():
        ppr = np.array([pprofs[ltp.upper()][frs[i]] \
            for i, fr in enumerate(samples[ltp]) if fr == 1 and i != 0])
        ppr = norm_profile(ppr).astype(np.float64)
        for pn, tbl in d.iteritems():
            prg.step()
            with_protein = np.vstack(tbl['no'], ppr)
            whitened = sp.cluster.vq.whiten(with_protein)
            code_book = sp.cluster.vq.kmeans(whitened, 2)
            # tbl['km'] = 
    prg.terminate()

def read_positives(basedir, fname = 'manual_positive.csv'):
    '''
    Reads manually annotated positive hits
    from file.
    '''
    with open(os.path.join(basedir, fname), 'r') as f:
        _result = [[c.strip() for c in l.split('\t')] for l in f.read().split('\n')]
    result = dict((ltp, []) for ltp in uniqList([x[0] for x in _result if x[0] != '']))
    for l in _result:
        if l[0] != '':
            result[l[0]].append(l)
    return result

def spec_sens(valids, ltp, pos, metric, asc):
    '''
    Calculates specificity, sensitivity, precision,
    false discovery rate as a function of critical
    values of a score, for one LTP for one mode.
    Returns dict of lists.
    '''
    result = {'spec': [], 'sens': [], 'prec': [], 'fdr': [], 'cutoff': [], 'n': 0}
    tbl = valids[ltp][pos]
    ioffset = 0 if pos == 'pos' else 6
    if len(tbl['std']) > 0:
        sort_all(tbl, metric, asc)
        p = len(tbl['std'])
        for i, cutoff in enumerate(tbl[metric]):
            tp = sum([1 for oi in tbl['i'][:i + 1] if oi in tbl['std']])
            fp = sum([1 for oi in tbl['i'][:i + 1] if oi not in tbl['std']])
            fn = sum([1 for oi in tbl['i'][i + 1:] if oi in tbl['std']])
            tn = sum([1 for oi in tbl['i'][i + 1:] if oi not in tbl['std']])
            result['cutoff'].append(cutoff)
            result['spec'].append(tn / float(tn + fp))
            result['sens'].append(tp / float(tp + fn))
            result['prec'].append(tp / float(tp + fp))
            result['fdr'].append(fp / float(tp + fp))
        result['n'] = p
    return result

def evaluate_scores(valids, stdltps):
    '''
    Calculates specificity, sensitivity, precision,
    false discovery rate as a function of critical
    values of each scores, for all LTPs with manual
    positive annotations, for all modes.
    Returns 4x embedded dicts of
    LTPs/modes/metrics/performance metrics.
    E.g. result['STARD10']['pos']['ktv']['sens']
    '''
    result = dict((ltp, {'pos': {}, 'neg': {}}) for ltp in stdltps)
    metrics = [
        ('Kendall\'s tau', 'ktv', False),
        ('Spearman corr.', 'spv', False),
        ('Pearson corr.', 'pev', False),
        ('Euclidean dist.', 'euv', True),
        ('Robust corr.', 'rcv', False),
        ('Goodman-Kruskal\'s gamma', 'gkv', False),
        ('Difference', 'dfv', True)
    ]
    for ltp in stdltps:
        for pos in ['pos', 'neg']:
            tbl = valids[ltp][pos]
            for m in metrics:
                result[ltp][pos][m[1]] = spec_sens(valids, ltp, pos, m[1], m[2])
    return result

def fractions_barplot(samples, pprofs, pprofs_original):
    fracs = ['a9', 'a10', 'a11', 'a12', 'b1']
    font_family = 'Helvetica Neue LT Std'
    sns.set(font = font_family)
    fig, axs = plt.subplots(8, 8, figsize = (20, 20))
    ltps = sorted(samples.keys())
    for i in xrange(len(ltps)):
        ax = axs[i / 8, i % 8]
        ltpname = ltps[i]
        ppr = [pprofs[ltpname][fr] for fr in fracs]
        ppr_o = [pprofs_original[ltpname][fr] for fr in fracs]
        #B6B7B9 gray
        #6EA945 mantis
        #007B7F teal
        col = ['#6EA945' if s == 1 else '#B6B7B9' if s is None else '#007B7F'\
            for s in samples[ltpname][1:]]
        ax.bar(np.arange(len(ppr)), ppr_o, color = col, alpha = 0.1, edgecolor = 'none')
        ax.bar(np.arange(len(ppr)), ppr, color = col, edgecolor = 'none')
        ax.set_xticks(np.arange(len(ppr)) + 0.4)
        ax.set_xticklabels(fracs)
        ax.set_title('%s protein conc.'%ltpname)
    fig.tight_layout()
    fig.savefig('protein_profiles.pdf')

# ltp.fractions_barplot(samples_upper, pprofs)


def plot_score_performance(perf):
    metrics = [
        ('Kendall\'s tau', 'ktv', False),
        ('Spearman corr.', 'spv', False),
        ('Pearson corr.', 'pev', False),
        ('Euclidean dist.', 'euv', True),
        ('Robust corr.', 'rcv', False),
        ('Goodman-Kruskal\'s gamma', 'gkv', False),
        ('Difference', 'dfv', True)
    ]
    perfmetrics = [
        ('sens', '#6EA945', 'Sensitivity'),
        ('spec', '#FCCC06', 'Specificity'),
        ('prec', '#DA0025', 'Precision'),
        ('fdr', '#007B7F', 'FDR')
    ]
    font_family = 'Helvetica Neue LT Std'
    sns.set(font = font_family)
    # plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])
    for ltp, d1 in perf.iteritems():
        for pos, d2 in d1.iteritems():
            fig, axs = plt.subplots(len(metrics), figsize = (10, 20), sharex = False)
            plt.suptitle('Performance of scores :: %s, %s'%(ltp, pos), size = 16)
            for i, m in enumerate(metrics):
                for pm in perfmetrics:
                    axs[i].plot(d2[m[1]]['cutoff'], d2[m[1]][pm[0]], 
                        '-', linewidth = 0.33, color = pm[1], label = pm[2])
                leg = axs[i].legend()
                axs[i].set_ylabel('Performance\nmetrics')
                axs[i].set_xlabel('%s value'%m[0])
                axs[i].set_title(m[0])
                if m[2]:
                    axs[i].set_xlim(axs[i].get_xlim()[::-1])
            fig.tight_layout()
            plt.subplots_adjust(top = 0.95)
            fig.savefig('score_performance_%s-%s.pdf'%(ltp, pos))
            plt.close(fig)

def plot_roc(perf):
    metrics = [
        ('Kendall\'s tau', 'ktv', False),
        ('Spearman corr.', 'spv', False),
        ('Pearson corr.', 'pev', False),
        ('Euclidean dist.', 'euv', True),
        ('Robust corr.', 'rcv', False),
        ('Goodman-Kruskal\'s gamma', 'gkv', False),
        ('Difference', 'dfv', True)
    ]
    def colors():
        _colors = ['#6EA945', '#FCCC06', '#DA0025', '#007B7F', '#454447', '#996A44']
        for c in _colors:
            yield c
    font_family = 'Helvetica Neue LT Std'
    sns.set(font = font_family)
    # plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])
    fig, axs = plt.subplots((len(metrics) + len(metrics) % 2) / 2, 2,
        figsize = (10, 20), sharex = False, sharey = False)
    plt.suptitle('Reciever operating characteristics', size = 16)
    for i, m in enumerate(metrics):
        c = colors()
        ax = axs[i/2][i%2]
        for ltp, d1 in perf.iteritems():
            for pos, d2 in d1.iteritems():
                ax.plot(1 - np.array(d2[m[1]]['spec']), np.array(d2[m[1]]['sens']), 
                    '-', linewidth = 0.33, color = c.next(), 
                    label = '%s-%s n = %u' % (ltp, pos, d2[m[1]]['n']))
        ax.plot([0,1], [0,1], '--', linewidth = 0.33, color = '#777777')
        leg = ax.legend(loc = 4)
        ax.set_ylabel('Sensitivity')
        ax.set_xlabel('1 - Specificity')
        ax.set_title('ROC :: %s'%m[0])
    fig.tight_layout()
    plt.subplots_adjust(top = 0.95)
    fig.savefig('roc.pdf')
    plt.close(fig)

def best_combined(valids, scores, best = 10, ubiquity_treshold = 5):
    for ltp, d in valids.iteritems():
        for pn, tbl in d.iteritems():
            sorted_scores = []
            for sc in scores:
                oi_sorted = np.copy(tbl['i'][tbl[sc[0]].argsort()])
                if sc[1]:
                    sorted_scores.append(oi_sorted)
                else:
                    sorted_scores.append(oi_sorted[::-1])
            all_ranks = np.empty((len(scores), len(tbl['i'])))
            for j, oi in enumerate(tbl['i']):
                for i, scs in enumerate(sorted_scores):
                    for si, sc in enumerate(scs):
                        if oi == sc:
                            if tbl['cpv'][j] < 1.0 or \
                                np.isnan(np.sum([tbl[s[0]][j] for s in scores])) or \
                                ('ubi' in tbl and tbl['ubi'][j] > ubiquity_treshold):
                                all_ranks[i,j] = np.inf
                            else:
                                all_ranks[i,j] = si
            combined_ranks = np.sum(all_ranks, axis = 0)
            best_oi = tbl['i'][np.argsort(combined_ranks)][:best]
            tbl['best%u'%best] = best_oi

def best_gk_eu(valids, best = 10):
    best_combined(valids, [('gkv', False), ('euv', True)], best = best)

def best_table(valids, fname, pos, best = 10):
    hdr = ['LTP', 'm/z', 'Database_ID', 'Level', 'Full_name', 'Formula', 'Adduct', 'Database_m/z', 'Adduct_m/z']
    sort_alll(valids, 'mz')
    with open(fname, 'w') as f:
        f.write('\t'.join(hdr) + '\n')
        for ltp, d in valids.iteritems():
            tbl = d[pos]
            for oi in tbl['best%u'%best]:
                i = np.where(tbl['i'] == oi)[0]
                mz = tbl['mz'][i]
                if tbl['lip'][oi] is not None:
                    for lip in tbl['lip'][oi]:
                        f.write('\t'.join([ltp, '%.08f'%mz] + [str(x) for x in list(lip)]) + '\n')
                else:
                    f.write('\t'.join([ltp, '%.08f'%mz] + ['unknown']*6) + '\n')

def true_positives(valids, stdpos, ltp, pos = 'pos', tolerance = 0.02):
    '''
    For one LTP having known binders looks up these
    among the valid features either in positive OR
    negative mode. (So you need to run twice, once
    for +, once for -.)
    '''
    ioffset = 0 if pos == 'pos' else 6
    valids[ltp][pos]['std'] = {}
    for fe in stdpos[ltp]:
        if fe[4 + ioffset] != '':
            mz = float(fe[ioffset + 4])
            iu = valids[ltp][pos]['mz'].searchsorted(mz)
            idx = []
            if iu < len(valids[ltp][pos]['mz']):
                u = 0
                while valids[ltp][pos]['mz'][iu + u] - mz <= tolerance:
                    idx.append(iu + u)
                    u += 1
            if iu > 0:
                l = 1
                while mz - valids[ltp][pos]['mz'][iu -l] <= tolerance:
                    idx.append(iu - l)
                    l += 1
            for i in idx:
                identify_feature(valids[ltp][pos], fe, i, ioffset)

def identify_feature(tbl, fe, i, ioffset):
    '''
    For one feature looks up if the lipid database ID
    and the adduct do match.
    Adds the matching record from the gold standard
    to lists in the `std` dict. There can be accessed
    by the original index.
    '''
    result = []
    oi = tbl['i'][i]
    if len(np.logical_and(tbl['lip'][oi][:,0] == fe[1 + ioffset], \
        tbl['lip'][oi][:,4] == fe[5 + ioffset])) > 0:
        if oi not in tbl['std']:
            tbl['std'][oi] = []
        tbl['std'][oi].append(fe)
    else:
        sys.stdout.write('Feature %.05f matched, lipid ID hasn\'t been found: %s' % \
            (fe[4 + ioffset], fe[1 + ioffset]))
        sys.stdout.flush()

def sortby_score(tbl, attr, asc = True):
    '''
    Returns sorted views of normalized features,
    m/z's and the attribute the sorting carried
    out by.
    '''
    sorted_no = tbl['no'][np.argsort(tbl[attr]),:]
    sorted_mz = tbl['mz'][np.argsort(tbl[attr])]
    sorted_attr = np.sort(tbl[attr])
    if asc:
        return sorted_no, sorted_mz, sorted_attr
    else:
        return sorted_no[::-1], sorted_mz[::-1], sorted_attr[::-1]

def sort_all(tbl, attr, asc = True):
    '''
    Sorts all arrays in one table by one specified 
    attribute.
    Either ascending or descending.
    '''
    ind = tbl[attr].argsort()
    dim = tbl[attr].shape[0]
    for k, a in tbl.iteritems():
        if k != attr and type(tbl[k]) == np.ndarray and tbl[k].shape[0] == dim:
            if len(a.shape) == 1:
                if asc:
                    tbl[k] = tbl[k][ind]
                else:
                    tbl[k] = tbl[k][ind][::-1]
            else:
                if asc:
                    tbl[k] = tbl[k][ind,:]
                else:
                    tbl[k] = tbl[k][ind,:][::-1]
    tbl[attr].sort()
    if not asc:
        tbl[attr] = tbl[attr][::-1]

def sort_alll(valids ,attr, asc = True):
    '''
    Sorts all arrays in all tables by one specified 
    attribute.
    Either ascending or descending.
    '''
    for d in valids.values():
        for tbl in d.values():
            sort_all(tbl, attr, asc)

def get_scored_hits(data):
    hits = val_ubi_prf_rpr_hits(data, ubiquity = 70, profile_best = 50000)
    [v[0].shape[0] if v is not None else None for vv in hits.values() for v in vv.values()]
    hits_upper = dict((l.upper(), {'pos': None, 'neg': None}) for l in hits.keys())
    for l, d in hits.iteritems():
        for pn, tbl in d.iteritems():
            if hits_upper[l.upper()][pn] is None:
                hits_upper[l.upper()][pn] = tbl
    return hits_upper

def lipid_lookup(stage0, runtime = None):
    '''
    Obtains full SwissLipids data
    '''
    pAdducts, nAdducts = get_swisslipids(swisslipids_url, 
        adducts = ['[M+H]+', '[M+NH4]+', '[M-H]-'], formiate = True)
    lipids = find_lipids(stage0, pAdducts, nAdducts)
    if runtime:
        runtime = timeit.timeit('find_lipids(final_hits, pAdducts, nAdducts)',
            setup = 'from __main__ import find_lipids, final_hits,'\
            'pAdducts, nAdducts', number = 1)
    return pAdducts, nAdducts, lipids, runtime

def lipid_lookup_exact(valids, swisslipids_url, exacts = None, runtime = None):
    '''
    Fetches data from SwissLipids and LipidMaps
    if not given.
    Looks up lipids based on exact masses (but
    taking into account the adducts).
    Writes hits into arrays in `lip` dict, having
    original indices as keys.
    Returns an array of lipid databases with exact
    masses.
    '''
    if exacts is None:
        la5Exacts = get_swisslipids_exact(swisslipids_url)
        la5Exacts = lipidmaps_exact(exacts = la5Exacts)
    else:
        la5Exacts = exacts
    find_lipids_exact(valids, la5Exacts)
    if runtime:
        runtime = timeit.timeit('find_lipids(valids, la5Exacts)',
            setup = 'from __main__ import find_lipids, valids, la5Exacts', number = 1)
    return la5Exacts, runtime

def evaluate_results(stage0, stage2, lipids, samples_upper, letter = 'e'):
    '''
    Tracks the number of features/lipids along stages.
    Does some plotting.
    '''
    stage3 = [(
            l.upper(), 
            [len(i['neg']), len(i['pos'])], 
            len(stage2[l.upper()]), 
            [len(stage0[l]['neg'][0]), len(stage0[l]['pos'][0])], 
            np.nansum([x for x in samples_upper[l] if x is not None])
        ) if len(i.values()) > 0 \
        else (
            l.upper(), 0, 0, 0, np.nansum([x for x in samples_upper[l] if x is not None])
        ) \
        for l, i in lipids.iteritems()]
    #
    [sys.stdout.write(str(i) + '\n') for i in stage3]
    sum([i[2] > 0 for i in stage3])
    #
    font_family = 'Helvetica Neue LT Std'
    sns.set(font = font_family)
    fig, ax = plt.subplots()
    color = ['#FCCC06', '#6EA945', '#007B7F']
    for f in [1, 2, 3]:
        x1 = [i[1][0] for i in stage3 if i[-1] == f]
        x2 = [i[1][1] for i in stage3 if i[-1] == f]
        y = [i[2] for i in stage3 if i[-1] == f]
        p = plt.scatter(x1, y, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
        p = plt.scatter(x2, y, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
        for xa, xb, yy in zip(x1, x2, y):
            p = plt.plot([xa, xb], [yy, yy], lw = .3, c = '#CCCCCC')
    plt.xlabel('Lipids (square: --; dot: +)')
    plt.ylabel('Lipids matching (square: -; dot: +)')
    fig.tight_layout()
    fig.savefig('hits-1-1000-%s.pdf'%letter)
    plt.close()
    #
    fig, ax = plt.subplots()
    color = ['#FCCC06', '#6EA945', '#007B7F']
    for f in [1, 2, 3]:
        x1 = [i[3][0] for i in stage3 if i[-1] == f]
        y1 = [i[1][0] for i in stage3 if i[-1] == f]
        x2 = [i[3][1] for i in stage3 if i[-1] == f]
        y2 = [i[1][1] for i in stage3 if i[-1] == f]
        p = plt.scatter(x1, y1, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
        p = plt.scatter(x2, y2, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
        for xa, xb, ya, yb in zip(x1, x2, y1, y2):
            p = plt.plot([xa, xb], [ya, yb], lw = .3, c = '#CCCCCC')
    #
    plt.xlabel('Num of m/z`s (green: --; yellow: +)')
    plt.ylabel('Lipids (green: --; yellow: +)')
    fig.tight_layout()
    fig.savefig('hits-2-1000-%s.pdf'%letter)
    plt.close()
    # 
    fig, ax = plt.subplots()
    color = ['#FCCC06', '#6EA945', '#007B7F']
    for f in [1, 2, 3]:
        x1 = [i[3][0] for i in stage3 if i[-1] == f]
        x2 = [i[3][1] for i in stage3 if i[-1] == f]
        y = [i[2] for i in stage3 if i[-1] == f]
        p = plt.scatter(x1, y, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
        p = plt.scatter(x2, y, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
        for xa, xb, yy in zip(x1, x2, y):
            p = plt.plot([xa, xb], [yy, yy], lw = .3, c = '#CCCCCC')
    plt.xlabel('Num of m/z`s (square: --; dot: +)')
    plt.ylabel('Lipids matching (square: -; dot: +)')
    fig.tight_layout()
    fig.savefig('hits-5-1000-%s.pdf'%letter)
    plt.close()
    #
    fig, ax = plt.subplots()
    y = flatList([i[3] for i in stage3 if type(i[3]) is not int])
    x = flatList([[i[4]]*2 for i in stage3 if type(i[3]) is not int])
    p = plt.scatter(x, y, c = '#6ea945', edgecolors = 'none')
    plt.xlabel('Fractions in sample')
    plt.ylabel('M/z min(pos, neg)')
    fig.tight_layout()
    fig.savefig('hits-3-150-%s.pdf'%letter)
    plt.close()
    #
    fig, ax = plt.subplots()
    y = [i[2] for i in stage3]
    x = [i[4] for i in stage3]
    p = plt.scatter(x, y, c = '#6ea945', edgecolors = 'none', alpha = .2)
    plt.xlabel('Fractions in sample')
    plt.ylabel('Lipids matched (pos, neg)')
    fig.tight_layout()
    fig.savefig('hits-4-1000-%s.pdf'%letter)
    plt.close()
    #
    fig, ax = plt.subplots()
    p = plt.bar(xrange(len(stage3)), 
        [i[2] for i in sorted(stage3, key = lambda x: x[2], reverse = True)], 
        0.5, color = '#6EA945', edgecolor = 'none')
    plt.xlabel('Lipid transfer protein')
    plt.ylabel('Number of valid lipid matches\n(with no filtering based on profiles)')
    ax.set_xticks(np.arange(len(stage3)) + .25)
    ax.set_xticklabels([i[0] for i in sorted(stage3, key = lambda x: x[2], reverse = True)], 
        rotation = 90, fontsize = 5)
    fig.tight_layout()
    fig.savefig('lipid-matches-3000-%s.pdf'%letter)
    plt.close()

# ## ## ## ## ## ## ## ## #
# The pipeline   ## ## ## #
# ## ## ## ## ## ## ## ## #

if __name__ == '__main__':
    #data, fnames, samples, csamples, pprofs = \
    #    init_from_scratch(basedir, ltpdirs, pptablef, samplesf)
    data, fnames, samples, csamples, samples_upper, pprofs = init_reinit(basedir)
    basic_filters(data, pprofs, samples, csamples)
    # stage0 :: feature filtering
    stage0 = get_scored_hits(data)
    # stage1 :: lipids
    #pAdducts, nAdducts, lipids, runtime = lipid_lookup(stage0, pAdducts, nAdducts)
    exacts, lipids, unknowns, runtime = lipid_lookup_exact(stage0, swisslipids_url)
    # stage2 :: positive-negative
    stage2 = negative_positive(lipids)
    stage2_unknown = negative_positive(unknowns)
    stage2_best = best_matches(lipids, stage2, minimum = 100000)
    stage2_best_unknown = best_matches(unknowns, stage2_unknown, minimum = 100000)
    stage2_best_all = best_matches((lipids, unknowns), (stage2, stage2_unknown), minimum = 100000)
    # output
    write_out(stage2_best, 'lipid_matches_nov25.csv')
    write_out(np.vstack((stage2_best, stage2_unknown)), 'all_sorted_nov25.csv')
    # evaluation
    evaluate_results(stage0, stage2, lipids, samples_upper, 'f')
    # LTP, num of lipid hits (min), num of positive-negative matching, 
    # num of selected features [neg, pos], num of fractions in sample
    [(l.upper(), p, len(j), len(stage2[l.upper()])) if j is not None else (l.upper(), p, 0, 0) \
        for l, i in lipids.iteritems() for p, j in i.iteritems()]
    # processing MS2
    pFragments = read_metabolite_lines('lipid_fragments_positive_mode.txt')
    nFragments = read_metabolite_lines('lipid_fragments_negative_mode.txt')
    ms2files = ms2_filenames(ltpdirs)
    ms2map = ms2_map(ms2files)
    ms2_main(ms2map, stage2_best, pFragments, nFragments)
    ms2_main(ms2map, stage2_best_unknown, pFragments, nFragments)

### ## ## ##
### E N D ##
### ## ## ##

    #swnames = set(i[2].split('(')[0].strip() for i in pAdducts if i[1] == 'Species')

    #sorted(list(swnames), key = lambda x: x.lower())

    #lipidnames = {
        #'PC': 'Phosphatidylcholine',
        #'PE': 'Phosphatidylethanolamine',
        #'Ceramide': 'Ceramide',
        #'PI': 'Phosphatidylinositol',
        #'PS': 'Phosphatidylserine',
        #'Palmitate': 'fatty acid',
        #'Oleate': 'fatty acid',
        #'Laurate': 'fatty acid',
        #'Stearate': 'fatty acid',
        #'Arachidonate': 'fatty acid',
        #'Sphingomyelin': ''
    #}

#'SLM:000335981', 'Structural subspecies',
#       'Monoacylglycerol (20:0/0:0/0:0)', 'C23H46O4', '[M+H]+', 387.346893]

'''
All lipid species names in SwissLipids:

set(['Lysophosphatidylcholine', 'Ganglioside GT2', 'Sulfohexosyl ceramide', 
'hexacosatetraenoate', 'tetracosahexaenoate', 'fatty acid', 'docosatrienoate', 
'decanoate', 'triacontapentaenoate', 'Phosphatidylinositol monophosphate', 
'Ceramide phosphate', 'hexadecanol', 'Ganglioside GA1', 'Diacylglycerol', 
'Triacylglycerol', 'eicosatrienoate', 'Lysophosphatidate', 'hexadecanoate', 
'Ganglioside GM4', 'octadecatetraenoate', 'octanoate', 
'dotriacontatetraenoate', 'Ganglioside GM2', 'tridecanol', 'Dihexosyl ceramide', 
'hexacosahexaenoate', 'tetracosanoate', 'eicosapentaenoate', 'octadecenol', 
'Monoalkyldiacylglycerol', 'octatriacontatetraenoate', 'Ganglioside GD1', 
'hexadecenoate', 'Ganglioside GD2', 'Ceramide', 'Phosphatidylglycerophosphate',
'Ganglioside GM3', 'octacosapentaenoate', 'eicosanol', 'Ganglioside GQ1', 
'octadecatrienoate', 'Ganglioside GM1', 'hexacosenoate', 'docosanoate', 
'octadecadienol', 'docosatetraenoate', 'heptadecanoate', 'octadecanoate', 
'Phosphatidylinositol bisphosphate', 'docosenol', 'tetracosatetraenoate', 
'hexadecadienoate', 'hexacosapentaenoate', 'Phosphatidylinositol', 'Ganglioside 
GD3', 'dotriacontahexaenoate', 'docosapentaenoate', 'triacontanoate', 
'tetradecanoate', 'pentadecanol', 'octacosanoate', 'Monoalkylglycerol', 
'eicosatetraenoate', 'Phosphatidylethanolamine', 'eicosenol', 
'Phosphatidylinositol trisphosphate', 'tetracosanol', 'docosadienoate', 
'tetracosenoate', 'tetratriacontapentaenoate', 'octacosatetraenoate', 
'Lysophosphatidylethanolamine', 'tetradecenoate', 'Lysophosphatidylglycerol', 
'docosenoate', 'octadecanol', 'Sulfodihexosyl ceramide', 'Ganglioside Gb3', 
'Phosphatidylglycerol', 'octacosanol', 'eicosenoate', 'Ganglioside GT1', 
'octadecenoate', 'hexacosanol', 'triacontatetraenoate', 
'octatriacontapentaenoate', 'tridecanoate', 'dotriacontapentaenoate', 
'Monoalkylmonoacylglycerol', 'Hexosyl ceramide', 'Lysophosphatidylserine', 
'heneicosanoate', 'octadecadienoate', 'Ganglioside GA2', 
'tetratriacontahexaenoate', 'Phosphatidylcholine', 'tetracosapentaenoate', 
'hexacosanoate', 'Bismonoacylglycerolphosphate', 'Sphingomyelin', 
'octacosahexaenoate', 'Ganglioside GP1', 'heptadecanol', 'Phosphatidylserine', 
'tetratriacontatetraenoate', 'docosanol', 'nonadecanoate', 'Ganglioside GT3', 
'dodecanoate', 'hexatriacontapentaenoate', 'Lysophosphatidylinositol', 
'hexatriacontahexaenoate', 'eicosanoate', 'Ceramide phosphoethanolamine', 
'triacontanol', 'pentadecanoate', 'tetradecanol', 'triacontahexaenoate', 
'Phosphatidate', 'hexanoate', 'eicosadienoate', 'Monoacylglycerol', 
'docosahexaenoate', 'hexatriacontatetraenoate'])
'''