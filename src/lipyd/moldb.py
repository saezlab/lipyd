#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function
from past.builtins import xrange, range, reduce

from future.utils import iteritems

import os
import sys
import imp
import re
import copy
import struct
import itertools
import collections
from argparse import Namespace

import numpy as np
import pandas as pd

try:
    import pybel
except:
    sys.stdout.write('\t:: No module `pybel` available.\n')

import lipyd._curl as _curl
import lipyd.common as common
import lipyd.settings as settings
import lipyd.mz as mzmod
import lipyd.progress as progress
import lipyd.sdf as sdf
import lipyd.lipid as lipid
import lipyd.lookup as lookup
import lipyd.name as lipidname
import lipyd.formula as formula
import lipyd.lipproc as lipproc


class Reader(object):
    
    def __init__(self):
        
        self.load()
        self.process()
        return self.__iter__()
    
    def reload(self, children = False):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def iterrows(self):
        
        for mol in self:
            
            yield [
                mol[2].data['LM_ID'],
                'Species',
                '|'.join(self.names(mol[2])),
                mol[2].data['INCHI_KEY'],
                '',
                mol[0],
                mol[2].data['PUBCHEM_CID']
            ]


class LipidMaps(sdf.SdfReader):
    
    def __init__(self, extract_file = True):
        
        self.url   = settings.get('lipidmaps_url')
        self.fname = settings.get('lipidmaps_fname')
        self.curl  = _curl.Curl(self.url, large = True, silent = False)
        
        if extract_file:
            
            self.efname = os.path.join('cache', self.fname.split('/')[-1])
            
            with open(self.efname, 'wb') as efp:
                
                for l in self.curl.result[self.fname]:
                    
                    efp.write(l)
            
            efp = open(os.path.join('cache', self.fname.split('/')[-1]), 'rb')
            sdf.SdfReader.__init__(self, efp)
        
        else:
            sdf.SdfReader.__init__(self, self.curl.result[self.fname])
        
        self.nameproc = lipidname.LipidNameProcessor(
            database = 'lipidmaps',
            iso = True
        )
    
    def __iter__(self):
        
        for rec in sdf.SdfReader.__iter__(self):
            
            if (
                'EXACT_MASS' not in rec['annot'] or
                float(rec['annot']['EXACT_MASS']) == 0
            ):
                
                try:
                    exmass = formula.Formula(rec['annot']['FORMULA']).mass
                except KeyError:
                    # if no exact mass it means
                    # this is a higher level category
                    continue
            else:
                exmass = float(rec['annot']['EXACT_MASS'])
            
            names = [
                rec['name'][nametype]
                for nametype in ('COMMON_NAME', 'SYSTEMATIC_NAME')
                if nametype in rec['name']
            ]
            if 'SYNONYMS' in rec['name']:
                names.extend(
                    n.strip() for n in rec['name']['SYNONYMS'].split(';')
                )
            
            hg, chainsum, chains = self.nameproc.process(names)
            
            liprec = lipproc.LipidRecord(
                lab = lipproc.LipidLabel(
                    db_id = rec['id'],
                    db    = 'LipidMaps',
                    names = tuple(names),
                ),
                hg  = hg,
                chainsum = chainsum,
                chains = chains,
            )
            
            yield exmass, liprec


class SwissLipids(Reader):
    
    def __init__(self, levels = set(['Species']), silent = False,
                 nameproc_args = None, branched = False,
                 exact_mass_formula_fallback = True):
        """
        Downloads and serves the SwissLipids database.
        
        Automatically downloads the data at the first time and stores it in a
        cache file to be read from there at next usage. Scans the entire file
        and builds multiple indices in order to quickly access records upon
        request. Provides a number of methods to retrieve records either as
        lines or openbabel OBMol instances.
        
        Args
        ----
        :param set levels:
            Levels in SwissLipids hierarchy. By default only "species".
        :param bool branched:
            Include lipids with branched alkyl chain (iso).
        :param dict nameproc_args:
            Arguments passed to the name processor.
        :param bool exact_mass_formula_fallback:
            If exact mass not available form SwissLipids calculate it from
            the formula. This is dangerous because the formula is sometimes
            dehydrogenated and charged state while exact mass should be
            uncharged with all hydrogenes
        """
        
        self.silent = silent
        self.exact_mass_formula_fallback = exact_mass_formula_fallback
        self.nameproc_args = nameproc_args or {}
        self.set_levels(levels)
        self.url = settings.get('swisslipids_url')
        self.load()
        self.make_index()
    
    def set_levels(self, levels):
        """
        Sets the levels to be processed. Levels in SwissLipids are `Species`,
        `Molecular subspecies`, `Structural subspecies` and
        `Isomeric subspecies`.
        
        Args
        ----
        :param set levels:
            A set of one or more of the levels above.
        """
        
        self.levels = levels
        self.init_name_processor()
    
    def init_name_processor(self):
        """
        Creates a `LipidNameProcessor` instance to process lipid names
        at indexing.
        """
        
        self.nameproc = lipidname.LipidNameProcessor(
            iso = 'Isomeric subspecies' in self.levels,
            **self.nameproc_args
        )
    
    def load(self):
        
        self.close_gzfile()
        
        self._curl = _curl.Curl(self.url, silent = False,
                                compr = 'gz', large = True)
        self._gzfile = self._curl.result
    
    def iterfile(self):
        
        self._plainfile.seek(0)
        _ = self._plainfile.readline()
        
        for line in self._plainfile:
            
            yield line
    
    @staticmethod
    def names(line):
        
        return '|'.join(line[2:5])
    
    def make_index(self):
        
        def cc2str(cc):
            
            return (
                '%s%s%u:%u' % (
                    cc[0],
                    '-' if cc[0] in {'O', 'P'} else '',
                    cc[1],
                    cc[2]
                )
            )
        
        self.close_plainfile()
        
        self.load()
        self.index = collections.defaultdict(lambda: set([]))
        self.hg_index      = collections.defaultdict(lambda: set([]))
        self.species_index = collections.defaultdict(lambda: set([]))
        self.subspec_index = collections.defaultdict(lambda: set([]))
        self.isomer_index  = collections.defaultdict(lambda: set([]))
        
        if not self.silent:
            
            self._gzfile.fileobj.seek(-4, 2)
            ucsize = struct.unpack('I', self._gzfile.fileobj.read(4))[0]
            self.prg = progress.Progress(ucsize, 'Indexing SwissLipids', 101)
        
        self._gzfile.fileobj.seek(0)
        
        self._plainfilename = '%s.extracted' % self._gzfile.name
        
        with open(self._plainfilename, 'wb') as fpp:
            
            offset = self._gzfile.tell()
            
            for l in self._gzfile:
                
                if not self.silent:
                    self.prg.step(len(l))
                
                ll = l.decode('ascii').split('\t')
                
                if ll[1] in self.levels:
                    
                    names = self.names(ll)
                    self.index[ll[0]].add(offset) # SwissLipids ID
                    self.index[ll[8]].add(offset) # SMILES
                    self.index[ll[10]].add(offset) # InChI key
                    
                    for n in names.split('|'):
                        
                        self.index[n].add(offset)
                    
                    hg, chainsum, chains = self.nameproc.process(names)
                    
                    if hg:
                        
                        self.hg_index[hg].add(offset)
                    
                    if chainsum:
                        
                        self.species_index[
                            lipproc.summary_str(hg, chainsum)
                        ].add(offset)
                    
                    if chains:
                        
                        self.subspec_index[
                            lipproc.full_str(hg, chains)
                        ].add(offset)
                    
                    if chains:
                        
                        self.isomer_index[
                            lipproc.full_str(hg, chains, iso = True)
                        ].add(offset)
                
                offset = self._gzfile.tell()
                fpp.write(l)
        
        if not self.silent:
            self.prg.terminate()
        
        self.index = dict(self.index)
        self._plainfile = open(self._plainfilename, 'r')
    
    def get_hg(self, hg):
        
        return self.get_record(hg, index = 'hg')
    
    def get_hg_obmol(self, hg):
        
        return self.get_obmol(hg, index = 'hg')
    
    def get_species(self, name):
        
        return self.get_record(name, index = 'species')
    
    def get_subspec(self, name):
        
        return self.get_record(name, index = 'subspec')
    
    def get_isomer(self, name):
        
        return self.get_record(name, index = 'isomer')
    
    def get_hg_obmol(self, hg):
        
        return self.get_obmol(hg, index = 'hg')
    
    def get_species_obmol(self, name):
        
        return self.get_obmol(name, index = 'species')
    
    def get_subspec_obmol(self, name):
        
        return self.get_obmol(name, index = 'subspec')
    
    def get_isomer_obmol(self, name):
        
        return self.get_obmol(name, index = 'isomer')
    
    def get_record(self, name, index = ''):
        
        indexname = '%s%sindex' % (index, '_' if index else '')
        index = getattr(self, indexname)
        
        if name in index:
            
            for offset in index[name]:
                
                self._plainfile.seek(offset)
                yield self._plainfile.readline().strip().split('\t')
    
    def get_obmol(self, name, index = ''):
        
        return [
            self.to_obmol(rec)
            for rec in self.get_record(name, index = index)
        ]
    
    @staticmethod
    def to_obmol(record):
        
        if not record[8]:
            
            return None
        
        return pybel.readstring('inchi', record[9])
    
    @staticmethod
    def add_annotations(mol, record, exact_mass_formula_fallback = True):
        
        mol.db_id = record[0]
        mol.name  = record[3]
        mol.title = tuple(record[2:5])
        mol.chebi = record[24] if len(record) > 24 else ''
        mol.lipidmaps = record[25] if len(record) > 25 else ''
        mol.hmdb = record[26] if len(record) > 26 else ''
        mol.smiles = record[8]
        mol.swl_exact_mass = float(record[14]) if record[14] else None
        if not mol.swl_exact_mass and exact_mass_formula_fallback:
            
            try:
                # note: this is dangerous because the formula is sometimes
                # dehydrogenated and charged state while exact mass
                # should be uncharged with all hydrogenes
                mol.swl_exact_mass = formula.Formula(record[11]).mass
            except KeyError:
                pass
        
        mol.swl_formula = record[11]
        mol.inchi = record[9]
        mol.inchikey = record[10]
        mol.level = record[1]
        
        return mol
    
    def itermol(self, obmol = False):
        """
        Iterates the database either by yielding `pybel.Molecule` objects
        with extra attributes or dummy objects with the same attributes
        containing the various IDs and structure representations
        (see code of `add_annotation`) for details.
        
        Args
        ----
        :param bool obmol:
            Yield `pybel.Molecule` objects. By default simple dummy objects
            produced.
        """
        
        nosmiles = 0
        
        for line in self.iterfile():
            
            line = line.strip().split('\t')
            
            if len(line) > 22 and line[1] in self.levels:
                
                if not line[8]:
                    nosmiles += 1
                    continue
                
                if obmol:
                    
                    mol = self.to_obmol(line)
                    
                else:
                    
                    mol = Namespace()
                
                mol = self.add_annotations(
                    mol, line, self.exact_mass_formula_fallback
                )
                
                yield mol
    
    def __iter__(self):
        
        for mol in self.itermol(obmol = False):
            
            if not mol.swl_exact_mass:
                continue
            
            hg, chainsum, chains = self.nameproc.process(mol.title)
            
            rec = lipproc.LipidRecord(
                lab = lipproc.LipidLabel(
                    db_id = mol.db_id,
                    db    = 'SwissLipids',
                    names = mol.title,
                ),
                hg = hg,
                chainsum = chainsum,
                chains = chains
            )
            
            # yielding mass and record
            yield mol.swl_exact_mass or np.nan, rec
    
    def __del__(self):
        
        self.close_gzfile()
        self.close_plainfile()
    
    def close_gzfile(self):
        
        if hasattr(self, '_gzfile'):
            
            self._gzfile.close()
    
    def close_plainfile(self):
        
        if hasattr(self, '_plainfile'):
            
            self._plainfile.close()
    
    def export_names(self, proc):
        
        with open('names.tmp', 'w') as fp:
            
            for i in self.__iter__():
                
                n = proc.process(i.title)
                
                fp.write('%s\t%s\n' % (i.title, str(n)))


class MoleculeDatabaseAggregator(object):
    
    def __init__(
            self,
            resources = None,
            tolerance = 20,
            fa_args = None,
            sph_args = None,
            build = True,
            verbose = False
        ):
        """
        Builds a database of molecules and provides methods for look up by
        masses and names. Metabolites are processed from databases like
        SwissLipids and LipidMaps and also autogenerated using classes
        defined in the `lipid` module.
        
        Args
        ----
        :param dict resources:
            Databases to use with arguments. Keys are database names, values
            are tuples of classes and arguments. Default SwissLipids and
            LipidMaps.
        :param int tolerance:
            Mass lookup tolerance in ppm.
        :param dict fa_args:
            Fatty acyl arguments for autogenerated metabolites.
        :param dict sph_args:
            Sphingosine base arguments for autogenerated metabolites.
        """
        
        self.verbose = verbose
        self.resources = resources or {
            'SwissLipids': (SwissLipids, {}),
            'LipidMaps': (LipidMaps, {})
        }
        self.tolerance = tolerance
        
        self.fa_args  = fa_args or {'c': (4, 36), 'u': (0, 10)}
        self.sph_args = sph_args or {'c': (16, 22), 'u': (0, 1)}
        
        if build:
            
            self.build()
    
    def reload(self, children = False):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def init_rebuild(self):
        """
        Creates an empty list where this object collects the masses and
        molecule annotations. This needs to be done before (re)building the
        database in order to start from an empty array.
        """
        
        self._mass_data = []
    
    def build(self):
        """
        Executes the workflow of the entire database building process.
        
        First loads the third party databases (SwissLipids and LipidMaps),
        then autogenerates glycerophospholipid, glycerolipid and sphingolipid
        series. At the end all data merged into common `masses` and `data`
        arrays and sorted by increasing mass. At this point the instance
        is able to do lookups.
        """
        
        self.init_rebuild()
        self.load_databases()
        self.auto_glycerophospholipids()
        self.auto_glycerolipids()
        self.auto_sphingolipids()
        self.auto_fattyacids()
        self.mass_data_arrays()
        self.sort()
    
    def load_databases(self):
        """
        Loads all databases and generates main array.
        """
        
        for cls, resargs in self.resources.values():
            
            resource = cls(**resargs)
            
            self._mass_data.extend(resource)
    
    def mass_data_arrays(self):
        
        if hasattr(self, '_mass_data'):
            
            self.masses = np.array(
                [i[0] for i in self._mass_data],
                dtype = np.float
            )
            self.data = np.zeros(len(self._mass_data), dtype = np.object)
            self.data[:] = [i[1] for i in self._mass_data]
            
            self.sort()
        
        delattr(self, '_mass_data')
    
    def auto_metabolites(
            self,
            fa_args = None,
            sph_args = None,
            sum_only = True,
            classes = None,
            **kwargs
        ):
        
        fa_args  = fa_args  or self.fa_args
        sph_args = sph_args or self.sph_args
        
        masses = []
        data   = []
        
        prg = progress.Progress(len(classes), 'Generating metabolites', 1)
        
        for clsname in classes:
            
            prg.step()
            
            if self.verbose:
                
                sys.stdout.write('\t:: Generating `%s`\n' % clsname)
            
            cls = getattr(lipid, clsname)
            
            self.add_metabolite_series(
                cls,
                fa_args = fa_args,
                sph_args = sph_args,
                sum_only = sum_only,
                **kwargs
            )
        
        prg.terminate()
    
    def add_metabolite_series(
            self,
            cls,
            fa_args = None,
            sph_args = None,
            sum_only = True,
            **kwargs
        ):
        """
        Adds metabolites of a single class generated along a defined
        range of homolog series.
        """
        
        if not hasattr(cls, '__call__'):
            
            if hasattr(cls, 'lower') and hasattr(lipid, cls):
                
                cls = getattr(lipid, cls)
                
            else:
                
                raise ValueError(
                    'Don\'t know how to '
                    'generate lipids from this: %s' % str(cls)
                )
        
        fa_args  = fa_args or self.fa_args
        sph_args = sph_args or self.sph_args
        
        gen = cls(
            fa_args  = copy.copy(fa_args),
            sph_args = copy.copy(sph_args),
            sum_only = sum_only,
            **kwargs
        )
        
        self._mass_data.extend(gen.iterlines())
    
    def auto_fattyacids(self, **kwargs):
        """
        Autogenerates all fatty acids from classes listed in
        `lipid.fattyacids`.
        
        Args
        ----
        :param **kwargs:
            Arguments for fatty acid classes: `c`, `u`, `fa_counts`, etc.
        """
        
        self.auto_metabolites(classes = lipid.fattyacids, **kwargs)
    
    def auto_sphingolipids(self, **kwargs):
        """
        Autogenerates all sphingolipids from classes listed in
        `lipid.sphingolipids`.
        
        Args
        ----
        :param **kwargs:
            Arguments for sphingolipid classes (`fa_args`, `sph_args`, etc).
        """
        
        self.auto_metabolites(classes = lipid.sphingolipids, **kwargs)
    
    def auto_glycerolipids(self, **kwargs):
        """
        Autogenerates all glycerolipids from classes listed in
        `lipid.glycerolipids`.
        
        Args
        ----
        :param **kwargs:
            Arguments for glycerolipid classes
            (`fa_args`, `sn2_fa_args`, etc).
        """
        
        self.auto_metabolites(classes = lipid.glycerolipids, **kwargs)
    
    def auto_glycerophospholipids(self, **kwargs):
        """
        Autogenerates all glycerophospholipids from classes listed in
        `lipid.glycerophospholipids`.
        
        Args
        ----
        :param **kwargs:
            Arguments for glycerophospholipid classes
            (`fa_args`, `sn2_fa_args`, etc).
        """
        
        self.auto_metabolites(classes = lipid.glycerophospholipids, **kwargs)
    
    def sort(self):
        """
        Sorts the `masses` and `data` arrays by increasing mass in order to
        make fast lookups possible.
        """
        
        self.data = self.data[self.masses.argsort()]
        self.masses.sort()
    
    def ilookup(self, m):
        
        return lookup.findall(self.masses, m, t = self.tolerance)
    
    def lookup(self, m):
        
        i = self.ilookup(m)
        
        return (
            self.masses[i],
            self.data[i]
        )
    
    def lookup_accuracy(self, m):
        """
        Performs a lookup and adds accuracy information to the result.
        """
        
        r = self.lookup(m)
        
        a = np.array([
            (m - rm) / m * 10**6 for rm in r[0]
        ])
        
        return r[0], r[1], a
    
    def adduct_lookup(self, mz, adducts = None, ionm = None, charge = None):
        """
        Does a series of lookups in the database assuming various adducts.
        Calculates the exact mass for the m/z for each possible adduct
        and searches these exact masses in the database.
        
        Returns a dict of tuples with 3-3 numpy arrays.
        Keys of the dict are adduct types. The arrays are exact masses,
        database record details and accuracies (ppm).
        """
        
        result = {}
        
        mz = mzmod.Mz(mz)
        charge = charge if charge is not None else 1 if ionm == 'pos' else -1
        
        if ionm in {'pos', 'neg'}:
            
            adducts = list(common.ad2ex[abs(charge)][ionm].keys())
        
        methods = dict((ad, common.exact_method[ad]) for ad in adducts)
        
        for ad, method in iteritems(methods):
            
            exmz = getattr(mz, method)()
            
            result[ad] = self.lookup_accuracy(exmz)
        
        return result
    
    def export_db(self, fname = 'molecule_database.tsv'):
        
        hdr = [
            'exact_mass', 'category', 'std_name', 'database_names',
            'database_id', 'database', 'prefix', 'sum_cc', 'sum_unsat',
            'fa1_cc', 'fa1_unsat', 'fa2_cc', 'fa2_unsat',
            'fa3_cc', 'fa3_unsat'
        ]
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('%s\n' % '\t'.join(hdr))
            
            for mass, data in zip(self.masses, self.data):
                
                _ = fp.write('%.12f\t%s\n' % (
                    mass,
                    '\t'.join(str(f) for f in data)
                ))
