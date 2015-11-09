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

import re
import timeit
import xlrd
from xlrd.biffh import XLRDError
import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import mass
import progress
import _curl

path_root = '/'
basedir = os.path.join(path_root, 'home', 'denes', 'Dokumentumok' , 'ltp')
ltpdirs = [os.path.join(basedir, 'share'),
    os.path.join(basedir, 'share', '2015_06_Popeye')]
samplesf = os.path.join(ltpdirs[0], 'control_sample.csv')
ppfracf = os.path.join(ltpdirs[0], 'fractions.csv')
ppsecdir = os.path.join(ltpdirs[0], 'SEC_profiles')
pptablef = os.path.join(basedir, 'proteins_by_fraction.csv')
swisslipids_url = 'http://www.swisslipids.org/php/export.php?action=get&file=lipids.csv'

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

def get_swisslipids(swisslipids_url, adducts = None, 
    formiate = True, exact_mass = False):
    if type(adducts) is list:
        adducts = set(adducts)
        if formiate:
            adducts.add('[M+OAc]+')
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

def adduct_lookup(mz, adducts, levels, tolerance = 0.02):
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
        tbl = read_xls(os.path.join(basedir, fname))[3:]
        minabs = min(0.0, min(to_float(l[5]) for l in tbl))
        for l in tbl:
            ml = to_float(l[4])
            if ml < fr[0]:
                continue
            if ml >= fr[1]:
                try:
                    fr = gfrac.next()
                except StopIteration:
                    break
            frac_abs[fr[2]].append(to_float(l[5]) - minabs)
        result[ltpname] = dict((fnum, np.mean(a)) for fnum, a in frac_abs.iteritems())
    return result

def write_pptable(pprops, pptablef):
    '''
    Writes protein profiles in a table, so we don't need to read
    all the 60 XLS files every time.
    '''
    with open(pptablef, 'w') as f:
        f.write('\t%s%s' % ('\t'.join(sorted(pprops.values()[0].keys(), 
            key = lambda x: (x[0], int(x[1:])))), '\n'))
        f.write('\n'.join('%s\t%s' % (ltp, 
                '\t'.join('%.20f'%d[fr] for fr in sorted(d.keys(), 
                    key = lambda x: (x[0], int(x[1:]))))) \
            for ltp, d in pprops.iteritems()))

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

def cprofile_filter(data, pprops, samples):
    profile_filter(data, pprops, samples, prfx = 'c')

def profile_filter(data, pprops, samples, prfx = ''):
    frs = ['c0', 'a9', 'a10', 'a11', 'a12', 'b1']
    notf = []
    cols = 'smp' if prfx == 'c' else 'lip'
    prg = progress.Progress(len(data) * 2, 'Profile filter', 1, percent = False)
    for ltp, d in data.iteritems():
        for pn, tbl in d.iteritems():
            prg.step()
            if ltp.upper() not in pprops:
                notf.append(ltp)
                continue
            # i != 0 : we always drop the blank control
            ppr = np.array([pprops[ltp.upper()][frs[i]] \
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
        return 1.0

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
    return (profile - np.nanmin(profile.astype(np.float64))) / np.nanmax(profile)

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

def ubiquity_filter(data, proximity = 0.02, only_valid = True):
    attr = 'ubi' if only_valid else 'uby'
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
                    if 'uby' not in ltp1t:
                        ltp1t['uby'] = np.zeros([ltp1t['raw'].shape[0], 1], dtype = np.int8)
                    if 'uby' not in ltp2t:
                        ltp2t['uby'] = np.zeros([ltp2t['raw'].shape[0], 1], dtype = np.int8)
                    for i1, mz1 in np.ndenumerate(ltp1t['raw'][:,1]):
                        if ltp1t['vld'][i1] or not only_valid:
                            i2u = ltp2t['raw'][:,1].searchsorted(mz1)
                            u = 0
                            while True:
                                if i2u + u < ltp2t['raw'].shape[0] and \
                                    ltp2t['raw'][i2u + u,1] - mz1 <= proximity:
                                    if ltp2t['vld'][i2u + u] or not only_valid:
                                        if i1 not in ltp1s:
                                            ltp1t['uby'][i1] += 1
                                            ltp1s.add(i1)
                                        if i2u + u not in ltp2s:
                                            ltp2t['uby'][i2u + u] += 1
                                            ltp2s.add(i2u + u)
                                    u += 1
                                else:
                                    break
                            l = 1
                            while True:
                                if i2u - l >= 0 and \
                                    mz1 - ltp2t['raw'][i2u - l,1] <= proximity:
                                    if ltp2t['vld'][i2u - l] or not only_valid:
                                        if i1 not in ltp1s:
                                            ltp1t['uby'][i1] += 1
                                            ltp1s.add(i1)
                                        if i2u - l not in ltp2s:
                                            ltp2t['uby'][i2u - l] += 1
                                            ltp2s.add(i2u - l)
                                    l += 1
                                else:
                                    break
    prg.terminate()

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
    prf = 'cprf' if tbl['lip'][1,:].count() == 1 else 'prf'
    _treshold = (treshold, tresholdB)
    if prf not in tbl:
        return None
    # dict of profile matching score values and their indices
    selected = np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
            tbl['qly'], tbl['crg']), tbl['are']), tbl['pks']), tbl['rpr']),
            tbl['ubi'][:,0] <= ubiquity), 
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
    return (tbl['raw'][indices,:], tbl['prf'][indices], tbl['cprf'][indices], tbl['qly'][indices], tbl['are'][indices], tbl['pks'][indices])

@get_hits
def pass_through(tbl, ubiquity = 7, treshold = 0.15, tresholdB = 0.25, profile_best = False):
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
    return (tbl['raw'][indices,:], tbl['prf'][indices], tbl['cprf'][indices], tbl['rpr'][indices], tbl['ubi'][indices])

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

def save(fnames, samples, pprops, basedir, fname = 'save.pickle'):
    pickle.dump((fnames, samples, pprops), open(os.path.join(basedir, fname), 'wb'))

def load(basedir, fname = 'save.pickle'):
    return pickle.load(open(os.path.join(basedir, fname), 'rb'))

def find_lipids(hits, pAdducts, nAdducts, levels = ['Species'], tolerance = 0.02):
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
                                (np.array([ltp.upper(), tbl[0][i,1], tbl[1][i], tbl[2][i]], 
                                    dtype = np.object), lip), axis = 0))
                lipids[ltp.upper()][pn] = np.vstack(sorted(result, key = lambda x: x[1]))
    return lipids

def negative_positive(lipids, tolerance = 0.02):
    result = dict((ltp.upper(), []) for ltp in lipids.keys())
    prg = progress.Progress(len(result), 'Matching positive & negative', 1, percent = False)
    for ltp, tbl in lipids.iteritems():
        prg.step()
        if 'neg' in tbl and 'pos' in tbl:
            for neg in tbl['neg']:
                add = neg[8]
                if add == '[M-H]-':
                    poshmz = Mz(Mz(neg[1]).add_h()).add_h()
                    posnh3mz = Mz(Mz(neg[1]).add_h()).add_nh4()
                elif add == '[M+Fo]-':
                    poshmz = Mz(Mz(neg[1]).remove_fo()).add_h()
                    posnh3mz = Mz(Mz(neg[1]).remove_fo()).add_nh4()
                else:
                    continue
                iu = tbl['pos'][:,1].searchsorted(poshmz)
                u = 0
                if iu < len(tbl['pos']):
                    while True:
                        if tbl['pos'][iu + u,1] - poshmz <= tolerance:
                            if tbl['pos'][iu + u,4] == neg[4] and \
                                tbl['pos'][iu + u,8] == '[M+H]+':
                                result[ltp].append(np.concatenate(
                                    (tbl['pos'][iu + u,1:], neg[1:]), axis = 0))
                            u += 1
                        else:
                            break
                if iu > 0:
                    l = 1
                    while iu >= l:
                        if poshmz - tbl['pos'][iu - l,1] <= tolerance:
                            if tbl['pos'][iu - l,4] == neg[4] and \
                                tbl['pos'][iu - l,8] == '[M+H]+':
                                result[ltp].append(np.concatenate(
                                    (tbl['pos'][iu - l,1:], neg[1:], ), axis = 0))
                            l += 1
                        else:
                            break
                iu = tbl['pos'][:,1].searchsorted(posnh3mz)
                u = 0
                if iu < len(tbl['pos']):
                    while True:
                        if tbl['pos'][iu + u,1] - posnh3mz <= tolerance:
                            if tbl['pos'][iu + u,4] == neg[4] and \
                                tbl['pos'][iu + u,8] == '[M+NH4]+':
                                result[ltp].append(np.concatenate(
                                    (tbl['pos'][iu + u,1:], neg[1:]), axis = 0))
                            u += 1
                        else:
                            break
                if iu > 0:
                    l = 1
                    while iu >= l:
                        if posnh3mz - tbl['pos'][iu - l,1] <= tolerance:
                            if tbl['pos'][iu - l,4] == neg[4] and \
                                tbl['pos'][iu - l,8] == '[M+NH4]+':
                                result[ltp].append(np.concatenate(
                                    (tbl['pos'][iu - l,1:], neg[1:]), axis = 0))
                            l += 1
                        else:
                            break
            if len(result[ltp]) > 0:
                result[ltp] = np.vstack(result[ltp])
    prg.terminate()
    return result

## functions for MS2

# from Toby Hodges:
def read_metabolite_lines(fname):
    Metabolites = []
    with open(fname, 'r') as Handle:
        for Line in Handle.readlines():
            Line = Line.strip()
            MetabMass, MetabType, MetabCharge = Line.split('\t')
            Metabolites.append([to_float(MetabMass), MetabType, MetabCharge])
            if '+' not in MetabCharge and '-' not in MetabCharge and 'NL' not in MetabCharge:
                sys.stdout.write('WARNING: fragment %s has no valid charge information!\n' % metabolites[metabolite][1])
    return np.array(Metabolites, dtype = np.object)

pFragments = read_metabolite_lines('lipid_fragments_positive_mode.txt')
nFragments = read_metabolite_lines('lipid_fragments_negative_mode.txt')

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
                        if ltpname not in fnames:
                            fnames[ltpname] = {}
                        if pos not in fnames[ltpname] or ltpd.endswith('update'):
                            fnames[ltpname][pos] = {}
                        fnames[ltpname][pos][fr] = os.path.join(*([d] + fpath + [f]))
        return fnames

ms2files = ms2_filenames(ltpdirs)

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
    for ltp, d in ms2files.iteritems():
        pFeatures = []
        nFeatures = []
        ultp = ltp.upper()
        for pn, files in d.iteritems():
            features = pFeatures if pn == stRpos else nFeatures
            #fractions = set([])
            for fr, fl in files.iteritems():
                fr = int(redgt.match(fr).groups()[0])
                #fractions.add(fr)
                offset = 0
                with open(fl, 'r') as f:
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
    return result

def ms2_main(ms2map, ms1matches, pFragments, nFragments):
    result = dict((ltp, {}) for ltp in ms1matches.keys())
    return result

def ms2_lookup(ms2map, ms1mz, fragments):
    pass

def ms2_collect(ms1mz, fragments):
    pass

ms2map = ms2_map(ms2files)

# ## ## ## ## ## ## ## ## #
# The pipeline   ## ## ## #
# ## ## ## ## ## ## ## ## #

#fnames = dict(get_filenames(ltpdirs[0]).items() + get_filenames(ltpdirs[1]).items())
#samples = read_samples(samplesf)
# at first run, after reading from saved textfile
#pprops = protein_profiles(ppsecdir, ppfracf)
#write_pptable(pprops, pptablef)
#pprops = read_pptable(pptablef)
#del fnames['ctrl']
#save(fnames, samples, pprops, basedir)
fnames, samples, pprops = load(basedir)
# at first run, after loading from pickle
#data = read_data(fnames, samples)
#save_data(data, basedir)
data = load_data(basedir)

csamples = dict((k, np.array([1 if x is not None else None for x in v])) for k, v in samples.iteritems())

# running and analysing filters

apply_filters(data)
validity_filter(data)
profile_filter(data, pprops, samples)
profile_filter(data, pprops, csamples, prfx = 'c')


filtr_results = {}
for f in ['quality', 'charge', 'area', 'peaksize', 'validity']:
    filtr_results[f] = eval_filter(data, f)

filtr_results['ubiquity'] = eval_filter(data, 'ubiquity', 
    param = {'only_valid': False},
    runtime = True, repeat = 1, number = 1, 
    hit = lambda x: x < 7)

ubiquity_filter(data, only_valid = False)

filtr_results['val_ubi'] = eval_filter(data, 'val_ubi', param = {'ubiquity': 7})

filtr_results['profile025'] = eval_filter(data, 'profile', 
    param = {'pprops': pprops, 'samples': samples, 'prfx': ''},
    runtime = True, repeat = 1, number = 1, 
    hit = lambda x: x <= 0.25)

filtr_results['cprofile025']

t, hits, phits, hitn, hitp = eval_filter(data, 'cprofile', 
    param = {'pprops': pprops, 'samples': csamples},
    runtime = True, repeat = 1, number = 1, 
    hit = lambda x: x <= 0.25)

t1, hits1, phits1, hitn1, hitp1 = eval_filter(data, 'profile', 
    param = {'pprops': pprops, 'samples': samples},
    runtime = True, repeat = 1, number = 1, 
    hit = lambda x: x <= 0.25)

t, hits, phits, hitn, hitp = eval_filter(data, 'profile', 
    param = {'pprops': pprops, 'samples': samples},
    runtime = False, repeat = 1, number = 1, 
    hit = lambda x: x <= 0.25)

filtr_results['rprofile'] = eval_filter(data, 'rprofile', runtime = True)
filtr_results['val_prf'] = eval_filter(data, 'val_prf', param = {'treshold': 0.25})
filtr_results['val_rpr'] = eval_filter(data, 'val_rpr')
filtr_results['val_ubi_prf'] = eval_filter(data, 'val_ubi_prf', param = {'treshold': 0.25, 'ubiquity': 7})
filtr_results['val_ubi_prf_rprf'] = eval_filter(data, 'val_ubi_prf_rprf', param = {'treshold': 0.25, 'ubiquity': 7})

# selecting (predicted) positive features

final_hits2 = val_ubi_prf_rpr_hits(data, ubiquity = 70, profile_best = 1000)
[v[0].shape[0] if v is not None else None for vv in final_hits2.values() for v in vv.values()]
final_hits_upper = dict((l.upper(), {'pos': None, 'neg': None}) for l in final_hits.keys())
for l, d in final_hits.iteritems():
    for pn, tbl in d.iteritems():
        if final_hits_upper[l.upper()][pn] is None:
            final_hits_upper[l.upper()][pn] = tbl

final_hits = final_hits_upper

final_hits2_upper = dict((l.upper(), {'pos': None, 'neg': None}) for l in final_hits2.keys())
for l, d in final_hits2.iteritems():
    for pn, tbl in d.iteritems():
        if final_hits2_upper[l.upper()][pn] is None:
            final_hits2_upper[l.upper()][pn] = tbl

final_hits2 = final_hits2_upper

samples_upper = dict((l.upper(), s) for l, s in samples.iteritems())

# number of hits:
# [len(j) if j is not None else (l, p) for l, i in final_hits.iteritems() for p, j in i.iteritems()]

# obtaining full SwissLipids data
pAdducts, nAdducts = get_swisslipids(swisslipids_url, 
    adducts = ['[M+H]+', '[M+NH4]+', '[M-H]-'], formiate = True)

fliptime = timeit.timeit('find_lipids(final_hits, pAdducts, nAdducts)', 
    setup = 'from __main__ import find_lipids, final_hits, pAdducts, nAdducts', number = 1)
lipids = find_lipids(final_hits2, pAdducts, nAdducts)

# number of lipids:
# [len(j) if j is not None else (l, p) for l, i in lipids.iteritems() for p, j in i.iteritems()]

pntime = timeit.timeit('negative_positive(lipids)', setup = 'from __main__ import negative_positive, lipids', number = 1)
pnmatched = negative_positive(lipids)

# number of matching lipids:
# [(p, len(l)) for p, l in pnmatched.iteritems()]

[(l.upper(), p, len(j), len(pnmatched[l.upper()])) if j is not None else (l.upper(), p, 0, 0) for l, i in lipids.iteritems() for p, j in i.iteritems()]
# LTP, num of lipid hits (min), num of positive-negative matching, 
# num of selected features [neg, pos], num of fractions in sample
filtered_matched = [(
        l.upper(), 
        [len(i['neg']), len(i['pos'])], 
        len(pnmatched[l.upper()]), 
        [len(final_hits2[l]['neg'][0]), len(final_hits2[l]['pos'][0])], 
        np.nansum([x for x in samples_upper[l] if x is not None])
    ) if len(i.values()) > 0 \
    else (
        l.upper(), 0, 0, 0, np.nansum([x for x in samples_upper[l] if x is not None])
    ) \
    for l, i in lipids.iteritems()]

[sys.stdout.write(str(i) + '\n') for i in filtered_matched]
sum([i[2] > 0 for i in filtered_matched])

pnum = 'e'

font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
fig, ax = plt.subplots()
color = ['#FCCC06', '#6EA945', '#007B7F']
for f in [1, 2, 3]:
    x1 = [i[1][0] for i in filtered_matched if i[-1] == f]
    x2 = [i[1][1] for i in filtered_matched if i[-1] == f]
    y = [i[2] for i in filtered_matched if i[-1] == f]
    p = plt.scatter(x1, y, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
    p = plt.scatter(x2, y, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
    for xa, xb, yy in zip(x1, x2, y):
        p = plt.plot([xa, xb], [yy, yy], lw = .3, c = '#CCCCCC')

plt.xlabel('Lipids (square: --; dot: +)')
plt.ylabel('Lipids matching (square: -; dot: +)')
fig.tight_layout()
fig.savefig('hits-1-1000-%s.pdf'%pnum)
plt.close()

fig, ax = plt.subplots()
color = ['#FCCC06', '#6EA945', '#007B7F']
for f in [1, 2, 3]:
    x1 = [i[3][0] for i in filtered_matched if i[-1] == f]
    y1 = [i[1][0] for i in filtered_matched if i[-1] == f]
    x2 = [i[3][1] for i in filtered_matched if i[-1] == f]
    y2 = [i[1][1] for i in filtered_matched if i[-1] == f]
    p = plt.scatter(x1, y1, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
    p = plt.scatter(x2, y2, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
    for xa, xb, ya, yb in zip(x1, x2, y1, y2):
        p = plt.plot([xa, xb], [ya, yb], lw = .3, c = '#CCCCCC')

plt.xlabel('Num of m/z`s (green: --; yellow: +)')
plt.ylabel('Lipids (green: --; yellow: +)')
fig.tight_layout()
fig.savefig('hits-2-1000-%s.pdf'%pnum)
plt.close()


fig, ax = plt.subplots()
color = ['#FCCC06', '#6EA945', '#007B7F']
for f in [1, 2, 3]:
    x1 = [i[3][0] for i in filtered_matched if i[-1] == f]
    x2 = [i[3][1] for i in filtered_matched if i[-1] == f]
    y = [i[2] for i in filtered_matched if i[-1] == f]
    p = plt.scatter(x1, y, marker = 's', c = color[f - 1], edgecolors = 'none', alpha = .5)
    p = plt.scatter(x2, y, marker = 'o', c = color[f - 1], edgecolors = 'none', alpha = .5)
    for xa, xb, yy in zip(x1, x2, y):
        p = plt.plot([xa, xb], [yy, yy], lw = .3, c = '#CCCCCC')

plt.xlabel('Num of m/z`s (square: --; dot: +)')
plt.ylabel('Lipids matching (square: -; dot: +)')
fig.tight_layout()
fig.savefig('hits-5-1000-%s.pdf'%pnum)
plt.close()


fig, ax = plt.subplots()
y = flatList([i[3] for i in filtered_matched if type(i[3]) is not int])
x = flatList([[i[4]]*2 for i in filtered_matched if type(i[3]) is not int])
p = plt.scatter(x, y, c = '#6ea945', edgecolors = 'none')
plt.xlabel('Fractions in sample')
plt.ylabel('M/z min(pos, neg)')
fig.tight_layout()
fig.savefig('hits-3-150-%s.pdf'%pnum)
plt.close()

fig, ax = plt.subplots()
y = [i[2] for i in filtered_matched]
x = [i[4] for i in filtered_matched]
p = plt.scatter(x, y, c = '#6ea945', edgecolors = 'none', alpha = .2)
plt.xlabel('Fractions in sample')
plt.ylabel('Lipids matched (pos, neg)')
fig.tight_layout()
fig.savefig('hits-4-1000-%s.pdf'%pnum)
plt.close()

with open('lipids_matched.csv', 'w') as f:
    hdr = ['LTP',
        'Positive_m/z', 'Positive_m/z_in_SwissLipids', 'Positive_adduct',
        'Negative_m/z', 'Negative_m/z_in_SwissLipids', 'Negative_adduct',
        'SwissLipids_AC', 'SwissLipids_formula', 'SwissLipids_name']
    f.write('\t'.join(hdr) + '\n')
    for ltp, tbl in pnmatched.iteritems():
        for l in tbl:
            f.write('\t'.join([ltp, str(l[0]), str(l[6]), l[5], str(l[7]), 
                str(l[13]), l[12], l[1], l[4], l[3]]) + '\n')