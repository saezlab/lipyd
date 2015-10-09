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
    # Some functions to use for processing the ANOVA results
    # from GDSC drug sensitivity screening.
#

import os
import sys
import re
import numpy as np

import mass
import progress

ltpdirs = ['/home/denes/Dokumentumok/ltp/share',
    '/home/denes/Dokumentumok/ltp/share/2015_06_Popeye']
samplesf = os.path.join(ltpdirs[0], 'control_sample.csv')

class MolWeight():
    
    # 
    # Thanks for https://github.com/bsimas/molecular-weight/blob/master/chemweight.py
    #
    
    def __init__(self, formula = None, h = 0, c = 0, o = 0, n = 0):
        self.reform = re.compile(r'([A-Za-z]+)([0-9]*)')
        if formula is None:
            formula = '%s%u%s%u%s%u%s%u'%('H', h, 'C', c, 'O', o, 'N', n)
        self.formula = formula
        self.calc_weight()
    
    def calc_weight(self):
        atoms = self.reform.findall(self.formula)
        w = 0.0
        for element, count in atoms:
            count = int(count or '1')
            w += mass.mass[element] * count
        self.weight = w

class Mz():
    
    def __init__(self, mz, z = 1, tolerance = 0.01):
        self.mz = mz
        self.z = z
        self.tol = tolerance
    
    def __eq__(self, other):
        return self.z == other.z and \
            self.mz > other.mz - self.tol and \
            self.mz < other.mz + self.tol
    
    def __str__(self):
        return 'm/z = %f' % self.mz
    
    def adduct(self, m):
        return self.mz + m / z
    
    def weight(self):
        return self.mz * self.z
    
    def remove_h(self):
        return self.adduct(-mass.mass['proton'])
    
    def add_h(self):
        return self.adduct(mass.mass['proton'])
    
    def add_oh(self):
        m = MolWeight('OH')
        return self.adduct(m)
    
    def add_fo(self):
        m = MolWeight('HCO2')
        return self.adduct(m)
    
    def add_ac(self):
        m = MolWeight('H3C2O2')
        return self.adduct(m)
    
    def add_nh4(self):
        m = MolWeight('NH3')
        return self.adduct(m)

def read_samples(fname):
    data = {}
    with open(fname, 'r') as f:
        null = f.readline()
        for l in f:
            l = l.split(',')
            data[l[0].replace('"', '')] = \
                np.array([to_int(x) if x != '' else None for x in l[1:]])
    return data

def get_filenames(loc):
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
    renum = re.compile(r'([-]?[0-9]*[\.]?[0-9]+[eE]?[-]?[0-9]*)')
    num = renum.match(num.strip())
    if num:
        return float(num.groups(0)[0])
    else:
        return num

def to_int(num):
    renum = re.compile(r'([-]?[0-9]+[\.]?[0-9]*)')
    num = renum.match(num.strip())
    if num:
        return int(num.groups(0)[0])
    else:
        return num

def float_lst(l):
    return [to_float(x) for x in l]

def read_file(fname):
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
    data = dict((ltp, {}) for ltp in fnames.keys())
    prg = progress.Progress(len(fnames) * 2, 'Reading files', 1)
    for ltp, pos_neg in fnames.iteritems():
        for p, fname in pos_neg.iteritems():
            prg.step()
            try:
                data[ltp][p] = {}
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
    prg.terminate()
    return data

def diff_profiles(p1, p2):
    # profiles are numpy arrays
    # of equal length
    return np.abs(p1 / p1.max() - p2 / p2.max()).sum()

def across_samples(data, threshold = 0.1):
    for ltp, p in data.iteritems():
        closests = 

def enr_profile(p, c):
    pass

# ## ## ## ## ## ## ## ## #
fnames = dict(get_filenames(ltpdirs[0]).items() + get_filenames(ltpdirs[1]).items())
samples = read_samples(samplesf)
del fnames['ctrl']
# to test:
# d = read_file_np(fnames.values()[0]['pos'], samples)
# array header: 
data = read_data(fnames, samples)



#with open('share/control_sample.csv', 'w') as f:
    #f.write(','.join(['', '9', '10', '11', '12', '0']) + '\n')
    #for ltp in sorted(fnames.keys(), key = lambda x: x.lower()):
        #f.write(','.join([ltp] + [''] * 5) + '\n')
