#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

from __future__ import print_function
from past.builtins import xrange, range, reduce

from future.utils import iteritems

import os
import sys
import imp
import re
import copy
import itertools
import collections
import functools
from argparse import Namespace

import numpy as np
import pandas as pd

try:
    import pybel
    if 'ipykernel' not in sys.modules and pybel.tk is None:
        try:
            import tkinter
            import PIL
            import PIL.ImageTk
            pybel.tk = tkinter
            pybel.PIL = PIL.Image
            pybel.piltk = PIL.ImageTk
        except:
            sys.stdout.write(
                '\t:: `PIL` or `tkinter` not available.\n'
                '     `pybel` won\'t be able to draw molecules.\n'
            )
except:
    sys.stdout.write(':: Module `pybel` not available.\n')

import lipyd._curl as _curl
import lipyd.common as common
import lipyd.settings as settings
import lipyd.mz as mzmod
import lipyd.progress as progress
import lipyd.sdf as sdf
import lipyd.lipid as lipid
import lipyd.lookup as _lookup
import lipyd.name as lipidname
import lipyd.formula as formula
import lipyd.lipproc as lipproc


class Reader(object):
    """ """
    
    def __init__(self):
        
        self.load()
        self.process()
        return self.__iter__()
    
    def reload(self, children = False):
        """

        Parameters
        ----------
        children :
             (Default value = False)

        Returns
        -------

        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def iterrows(self):
        """ """
        
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
    """ """
    
    def __init__(self, extract_file = True):
        
        self.url   = settings.get('lipidmaps_url')
        self.fname = settings.get('lipidmaps_fname')
        self.curl  = _curl.Curl(self.url, large = True, silent = False)
        
        if extract_file:
            
            self.efname = os.path.join(
                settings.get('cachedir'),
                self.fname.split('/')[-1]
            )
            
            with open(self.efname, 'wb') as efp:
                
                for l in self.curl.result[self.fname]:
                    
                    efp.write(l)
            
            efp = open(
                os.path.join(
                    settings.get('cachedir'),
                    self.fname.split('/')[-1]
                ),
                'rb'
            )
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
                rec['name'][nametype].strip()
                for nametype in ('COMMON_NAME', 'SYSTEMATIC_NAME')
                if nametype in rec['name']
            ]
            if 'SYNONYMS' in rec['name']:
                names.extend(
                    n.strip() for n in rec['name']['SYNONYMS'].split(';')
                )
            
            names = [n.strip() for n in names if n.strip()]
            
            hg, chainsum, chains = self.nameproc.process(names)
            
            liprec = lipproc.LipidRecord(
                lab = lipproc.LipidLabel(
                    db_id   = rec['id'],
                    db      = 'LipidMaps',
                    names   = tuple(names),
                    formula = rec['annot']['FORMULA'],
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
        Sets the levels to be processed. Levels in SwissLipids are
        `Species`, `Molecular subspecies`, `Structural subspecies` and
        `Isomeric subspecies`.
        
        Parameters
        ----------
        levels : set
            A set of one or more of the levels above.
        """
        
        if isinstance(levels, common.basestring):
            
            levels = [levels]
        
        self.levels = set(l.capitalize() for l in levels)
        self.init_name_processor()
    
    def init_name_processor(self):
        """
        Creates a ``LipidNameProcessor`` instance to process lipid names
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
        self._gzfile = self._curl.gzfile
    
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
            
            self.prg = progress.Progress(self._curl.size, 'Indexing SwissLipids', 101)
        
        self._plainfilename = '%s.extracted' % self._gzfile.name
        
        with open(self._plainfilename, 'wb') as fpp:
            
            offset = self._gzfile.tell()

            for l in self._gzfile:
                
                if not self.silent:
                    self.prg.step(len(l))
                
                ll = l.decode('utf-8').split('\t')
                
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
                    
                    if hg and chainsum:
                        
                        self.species_index[
                            lipproc.summary_str(hg, chainsum)
                        ].add(offset)
                    
                    if hg and chains:
                        
                        self.subspec_index[
                            lipproc.full_str(hg, chains)
                        ].add(offset)
                    
                    if hg and chains:
                        
                        self.isomer_index[
                            lipproc.full_str(hg, chains, iso = True)
                        ].add(offset)
                
                offset = self._gzfile.tell()
                fpp.write(l)
        
        if not self.silent:
            self.prg.terminate()
        
        self.index = dict(self.index)
        self._plainfile = open(self._plainfilename, 'r')
    
    def get_hg(self, hg, sub = ()):
        """

        Parameters
        ----------
        hg :
            
        sub :
             (Default value = ())

        Returns
        -------

        """
        
        if isinstance(hg, common.basestring):
            
            hg = lipproc.Headgroup(hg, sub = sub)
        
        return self.get_record(hg, index = 'hg')
    
    def get_species(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_record(name, index = 'species')
    
    def get_subspec(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_record(name, index = 'subspec')
    
    def get_isomer(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_record(name, index = 'isomer')
    
    def get_hg_obmol(self, hg, sub = ()):
        """

        Parameters
        ----------
        hg :
            
        sub :
             (Default value = ())

        Returns
        -------

        """
        
        if isinstance(hg, common.basestring):
            
            hg = lipproc.Headgroup(hg, sub = sub)
        
        return self.get_obmol(hg, index = 'hg')
    
    def get_species_obmol(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_obmol(name, index = 'species')
    
    def get_subspec_obmol(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_obmol(name, index = 'subspec')
    
    def get_isomer_obmol(self, name):
        """

        Parameters
        ----------
        name :
            

        Returns
        -------

        """
        
        return self.get_obmol(name, index = 'isomer')
    
    def get_record(self, name, index = ''):
        """

        Parameters
        ----------
        name :
            
        index :
             (Default value = '')

        Returns
        -------

        """
        
        indexname = '%s%sindex' % (index, '_' if index else '')
        index = getattr(self, indexname)
        
        if name in index:
            
            for offset in index[name]:
                
                self._plainfile.seek(offset)
                yield self._plainfile.readline().strip().split('\t')
    
    def get_obmol(self, name, index = ''):
        """

        Parameters
        ----------
        name :
            
        index :
             (Default value = '')

        Returns
        -------

        """
        
        for rec in self.get_record(name, index = index):
            
            obmol = self.to_obmol(rec)
            self.add_annotations(obmol, rec)
            yield obmol
    
    @staticmethod
    def to_obmol(record):
        """

        Parameters
        ----------
        record :
            

        Returns
        -------

        """
        
        if not record[9] or record[9] == 'InChI=none':
            
            if record[8]:
                
                # processing from SMILES
                return pybel.readstring('smi', record[8])
            
            return None
        
        # processing from InChI
        return pybel.readstring('inchi', record[9])
    
    @staticmethod
    def add_annotations(mol, record, exact_mass_formula_fallback = True):
        """

        Parameters
        ----------
        mol :
            
        record :
            
        exact_mass_formula_fallback :
             (Default value = True)

        Returns
        -------

        """
        
        mol.db_id = record[0]
        mol.name  = record[3]
        if hasattr(mol, 'OBMol'):
            mol.OBMol.SetTitle(record[2])
        else:
            mol.title = tuple(n.strip() for n in record[2:5] if n.strip())
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
        """Iterates the database either by yielding `pybel.Molecule` objects
        with extra attributes or dummy objects with the same attributes
        containing the various IDs and structure representations
        (see code of `add_annotation`) for details.
        
        Args
        ----

        Parameters
        ----------
        bool :
            obmol:
            Yield `pybel.Molecule` objects. By default simple dummy objects
            produced.
        obmol :
             (Default value = False)

        Returns
        -------

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
                    db_id   = mol.db_id,
                    db      = 'SwissLipids',
                    names   = mol.title,
                    formula = mol.swl_formula,
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
        """ """
        
        if hasattr(self, '_gzfile'):
            
            self._gzfile.close()
    
    def close_plainfile(self):
        """ """
        
        if hasattr(self, '_plainfile'):
            
            self._plainfile.close()
    
    def export_names(self, proc):
        """

        Parameters
        ----------
        proc :
            

        Returns
        -------

        """
        
        with open('names.tmp', 'w') as fp:
            
            for i in self.__iter__():
                
                n = proc.process(i.title)
                
                fp.write('%s\t%s\n' % (i.title, str(n)))


class MoleculeDatabaseAggregator(object):
    
    
    def __init__(
            self,
            resources = None,
            tolerance = None,
            fa_args = None,
            sph_args = None,
            build = True,
            verbose = False,
            database_preference = None,
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
        self.tolerance = tolerance or settings.get('ms1_tolerance')
        
        self.fa_args  = fa_args or {'c': (4, 36), 'u': (0, 10)}
        self.sph_args = sph_args or {'c': (16, 22), 'u': (0, 1)}
        
        self.database_preference = (
            database_preference or settings.get('database_preference')
        )
        
        if build:
            
            self.build()
    
    def reload(self, children = False):
        """

        Parameters
        ----------
        children :
             (Default value = False)

        Returns
        -------

        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def init_rebuild(self):
        """Creates an empty list where this object collects the masses and
        molecule annotations. This needs to be done before (re)building the
        database in order to start from an empty array.

        Parameters
        ----------

        Returns
        -------

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
        self.auto_misc()
        self.mass_data_arrays()
        self.sort()
        self.build_names()
    
    def load_databases(self):
        """Loads all databases and generates main array."""
        
        for cls, resargs in self.resources.values():
            
            resource = cls(**resargs)
            
            self._mass_data.extend(resource)
    
    def mass_data_arrays(self):
        """ """
        
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
        """

        Parameters
        ----------
        fa_args :
             (Default value = None)
        sph_args :
             (Default value = None)
        sum_only :
             (Default value = True)
        classes :
             (Default value = None)
        **kwargs :
            

        Returns
        -------

        """
        
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
        """Adds metabolites of a single class generated along a defined
        range of homolog series.

        Parameters
        ----------
        fa_args :
             (Default value = None)
        sph_args :
             (Default value = None)
        sum_only :
             (Default value = True)
        **kwargs :
            

        Returns
        -------

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
        """Autogenerates all fatty acids from classes listed in
        `lipid.fattyacids`.
        
        Args
        ----

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        self.auto_metabolites(classes = lipid.fattyacids, **kwargs)
    
    def auto_misc(self, **kwargs):
        """Autogenerates all miscellanous classes listed in `lipid.misc`.
        
        Args
        ----

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        self.auto_metabolites(classes = lipid.misc, **kwargs)
    
    def auto_sphingolipids(self, **kwargs):
        """Autogenerates all sphingolipids from classes listed in
        `lipid.sphingolipids`.
        
        Args
        ----

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        self.auto_metabolites(classes = lipid.sphingolipids, **kwargs)
    
    def auto_glycerolipids(self, **kwargs):
        """Autogenerates all glycerolipids from classes listed in
        `lipid.glycerolipids`.
        
        Args
        ----

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        self.auto_metabolites(classes = lipid.glycerolipids, **kwargs)
    
    def auto_glycerophospholipids(self, **kwargs):
        """Autogenerates all glycerophospholipids from classes listed in
        `lipid.glycerophospholipids`.
        
        Args
        ----

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        self.auto_metabolites(classes = lipid.glycerophospholipids, **kwargs)
    
    def sort(self):
        """Sorts the `masses` and `data` arrays by increasing mass in order to
        make fast lookups possible.

        Parameters
        ----------

        Returns
        -------

        """
        
        self.data = self.data[self.masses.argsort()]
        self.masses.sort()
    
    def ilookup(self, m, tolerance = None):
        """

        Parameters
        ----------
        m :
            
        tolerance :
             (Default value = None)

        Returns
        -------

        """
        
        return _lookup.findall(
            self.masses,
            m,
            t = tolerance or self.tolerance,
        )
    
    def lookup(self, m, tolerance = None):
        """

        Parameters
        ----------
        m :
            
        tolerance :
             (Default value = None)

        Returns
        -------

        """
        
        i = self.ilookup(m, tolerance = tolerance)
        
        return (
            self.masses[i],
            self.data[i],
        )
    
    def lookup_accuracy(self, m, tolerance = None):
        """Performs a lookup and adds accuracy information to the result.

        Parameters
        ----------
        m :
            
        tolerance :
             (Default value = None)

        Returns
        -------

        """
        
        r = self.lookup(m, tolerance = tolerance)
        
        a = np.array([
            (m - rm) / m * 10**6 for rm in r[0]
        ])
        
        return r[0], r[1], a
    
    def adduct_lookup(
            self,
            mz,
            adducts = None,
            ionmode = None,
            charge = None,
            adduct_constraints = True,
            tolerance = None,
        ):
        """Does a series of lookups in the database assuming various adducts.
        Calculates the exact mass for the m/z for each possible adduct
        and searches these exact masses in the database.
        
        Returns a dict of tuples with 3-3 numpy arrays.
        Keys of the dict are adduct types. The arrays are exact masses,
        database record details and accuracies (ppm).

        Parameters
        ----------
        mz :
            
        adducts :
             (Default value = None)
        ionmode :
             (Default value = None)
        charge :
             (Default value = None)
        adduct_constraints :
             (Default value = True)
        tolerance :
             (Default value = None)

        Returns
        -------

        """
        
        result = {}
        
        mz = mzmod.Mz(mz)
        charge = (
            charge
                if charge is not None else
            1
                if ionmode == 'pos' else
            -1
        )
        
        if not adducts and ionmode in {'pos', 'neg'}:
            
            # we look up all adducts we have a method for
            adducts = list(settings.get('ex2ad')[abs(charge)][ionmode].keys())
        
        ad_default = settings.get('adducts_default')[ionmode][abs(charge)]
        ad_constr  = settings.get('adduct_constraints')[ionmode]
        
        exmethods = settings.get('ad2ex')[abs(charge)][ionmode]
        methods = dict((ad, exmethods[ad]) for ad in adducts)
        
        for ad, method in iteritems(methods):
            
            exmz = getattr(mz, method)()
            
            res  = self.lookup_accuracy(exmz, tolerance = tolerance)
            
            if adduct_constraints:
                
                ires = tuple(
                    i for i in xrange(res[0].shape[0])
                    if (
                        (
                            res[1][i].hg not in ad_constr and
                            ad in ad_default
                        ) or
                        (
                            res[1][i].hg in ad_constr and
                            ad in ad_constr[res[1][i].hg]
                        )
                    )
                )
                
                res = (res[0][ires,], res[1][ires,], res[2][ires,])
            
            if len(res[0]):
                
                result[ad] = res
        
        return result
    
    def adduct_lookup_many(
            self,
            mzs,
            adducts = None,
            ionmode = None,
            charge = None,
            adduct_constraints = True,
            tolerance = None,
        ):
        """Performs the lookup on a vector of m/z values.
        Calls the ``adduct_lookup`` method on all m/z's.
        
        Returns array of dicts with lookup results.

        Parameters
        ----------
        mzs :
            
        adducts :
             (Default value = None)
        ionmode :
             (Default value = None)
        charge :
             (Default value = None)
        adduct_constraints :
             (Default value = True)
        tolerance :
             (Default value = None)

        Returns
        -------

        """
        
        result = []
        
        for mz in mzs:
            
            result.append(
                self.adduct_lookup(
                    mz,
                    adducts = adducts,
                    ionmode = ionmode,
                    charge = charge,
                    adduct_constraints = adduct_constraints,
                    tolerance = tolerance,
                )
            )
        
        return np.array(result)
    
    def export_db(self, fname = 'molecule_database.tsv'):
        """

        Parameters
        ----------
        fname :
             (Default value = 'molecule_database.tsv')

        Returns
        -------

        """
        
        hdr = [
            'exact_mass', 'category', 'std_name', 'database_names',
            'database_id', 'database', 'prefix', 'sum_cc', 'sum_unsat',
            'fa1_cc', 'fa1_unsat', 'fa2_cc', 'fa2_unsat',
            'fa3_cc', 'fa3_unsat',
        ]
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('%s\n' % '\t'.join(hdr))
            
            for mass, data in zip(self.masses, self.data):
                
                _ = fp.write('%.12f\t%s\n' % (
                    mass,
                    '\t'.join(str(f) for f in data)
                ))
    
    
    def export_db_lipidblast(self, fname = 'molecule_database.csv'):
        """
        Exports the lipid metabolite database in LipidBlast format.
        
        Parameters
        ----------
        fname : str
            File name to export the database to.
        """
        
        hdr = [
            'Retention Time (min)', 'Neutral Mass', 'Compound ID',
            'Description', 'Formula', 'URL',
        ]
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('%s\n' % ','.join(hdr))
            
            for mass, data in zip(self.masses, self.data):
                
                try:
                    
                    std_name = '%s;' % data.subspecies_str()
                    
                except AttributeError:
                    
                    std_name = ''
                
                _ = fp.write(
                    '%s\n' % ','.join((
                        '', # we don't know RT, leaving it empty
                        '%.012f' % mass,
                        (
                            'lipyd.lipid'
                                if data.lab.db == 'lipyd.lipid' else
                            data.lab.db_id
                        ),
                        '"%s%s"' % (
                            std_name,
                            data.lab.names[0],
                        ),
                        data.lab.formula,
                        (
                            'http://saezlab.github.io/lipyd'
                                if data.lab.db == 'lipyd.lipid' else
                            self.get_url(data.lab.db_id)
                        ),
                    ))
                )
    
    
    @staticmethod
    def get_url(db_id):
        """

        Parameters
        ----------
        db_id :
            

        Returns
        -------

        """
        
        if db_id[:3] == 'SLM':
            
            return 'http://swisslipids.org/#/entity/%s' % db_id
            
        elif db_id[:2] == 'LM':
            
            return 'http://lipidmaps.org/data/LMSDRecord.php?LMID=%s' % db_id
    
    @staticmethod
    def records_string(
            records,
            adducts = None,
            databases = None,
            show_ppm = False,
            show_adduct = False,
            show_db = False,
            use_db_names = False,
        ):
        """

        Parameters
        ----------
        records :
            
        adducts :
             (Default value = None)
        databases :
             (Default value = None)
        show_ppm :
             (Default value = False)
        show_adduct :
             (Default value = False)
        show_db :
             (Default value = False)

        Returns
        -------

        """
        
        result = set()
        
        for add, data in iteritems(records):
            
            if adducts is not None and add not in adducts:
                
                continue
            
            for rec_mz, rec, rec_ppm in zip(*data):
                
                if databases is not None and rec.lab.db not in databases:
                    
                    continue
                
                name = (
                    lipproc.summary_str(rec.hg, rec.chainsum)
                        if not use_db_names and rec.hg else
                    rec.lab.names[0]
                        if rec.lab.names else
                    'Unknown'
                )
                
                details = []
                
                if show_ppm:
                    details.append('%.01fppm' % rec_ppm)
                
                if show_adduct:
                    details.append(add)
                
                if show_db:
                    details.append(rec.lab.db)
                
                result.add('%s%s%s%s' % (
                    name,
                    '[' if details else '',
                    ','.join(details),
                    ']' if details else '',
                ))
        
        return ';'.join(result)
    
    
    def build_names(self):
        """
        Builds a dictionary for names to index or mass lookup.
        """
        
        names = collections.defaultdict(set)
        
        for i, rec in enumerate(self.data):
            
            if not rec.hg:
                
                continue
            
            name = rec.summary_str()
            
            names[name].add(i)
        
        self.names = dict(
            (
                name,
                np.array(sorted(idx))
            )
            for name, idx in iteritems(names)
        )
    
    
    def idx_from_name(
            self,
            name = None,
            hg = None,
            chainsum = None,
            chains = (),
        ):
        """
        For a lipid name returns the indices of the matching records.
        
        Returns
        -------
        Set of indices. If name not in the database returns None.
        """
        
        if name is None and hg is not None:
            
            if chainsum is None:
                
                chainsum = lipproc.sum_chains(chains)
            
            name = lipproc.summary_str(hg = hg, chainsum = chainsum)
        
        if name in self.names:
            
            return self.names[name]
    
    
    def masses_from_name(self, name = None, **kwargs):
        """
        For a lipid name looks up the corresponding exact masses.
        
        Returns
        -------
        Array of masses.
        """
        
        idx = self.idx_from_name(name = name, **kwargs)
        
        if idx is not None:
            
            return self.masses[idx]
    
    
    def mass_from_name(
            self,
            name = None,
            database_preference = None,
            **kwargs,
        ):
        """
        For a lipid name returns one exact mass, preferably the one from
        the database in front of the ``database_preference`` list.
        """
        
        database_preference = database_preference or self.database_preference
        
        idx = self.idx_from_name(name = name, **kwargs)
        
        if idx is not None:
            
            for db in database_preference:
                
                for i in idx:
                    
                    rec = self.data[i]
                    
                    if rec.lab.db == db:
                        
                        return self.masses[i]
    
    
    def mz_from_name(
            self,
            adduct,
            name = None,
            database_preference = None,
            **kwargs,
        ):
        
        exmass = self.mass_from_name(
            name = name,
            database_preference = database_preference,
            **kwargs,
        )
        
        if exmass is not None:
            
            adduct_method = (
                settings.get('ex2ad_all')[adduct]
            )
            
            return getattr(formula.Formula(exmass), adduct_method)()
    
    
    def mz_lowest_error_from_name(
            self,
            measured_mz,
            adduct,
            name = None,
            **kwargs,
        ):
        """
        Regarding a measured m/z and an assumed adduct type and name
        returns the m/z of the record with matching name and lowest error.
        """
        
        exmasses = self.masses_from_name(name = name, **kwargs)
        
        if exmasses is not None:
            
            adduct_method = (
                settings.get('ex2ad_all')[adduct]
            )
            
            addmasses = np.array([
                getattr(formula.Formula(exmass), adduct_method)()
                for exmass in exmasses
            ])
            
            ppms = np.array([
                common.ppm(addmass, measured_mz)
                for addmass in addmasses
            ])
            
            return addmasses[np.argmin(np.abs(ppms))]
    
    
    def name_inconsistent(self, name = None, **kwargs):
        """
        Tells if a name has different corresponding masses in the names
        dictionary. If a name not in the database still returns False.
        If the name inconsistent returns True.
        """
        
        masses = self.masses_from_name(name = name, **kwargs)
        
        return masses is not None and max(masses) - min(masses) > 1e-6
    
    
    def get_inconsistent_names(self):
        """
        Collects all inconsistent names.
        """
        
        inconsistent = []
        
        for name in self.names.keys():
            
            if self.name_inconsistent(name = name):
                
                inconsistent.append(name)
        
        return inconsistent
    
    
    def _dec_lookup_daltons_tolerance(self, func):
        """
        Calls a lookup method with tolerance given in Daltons instead of
        ppms.
        """
        
        


def init_db(**kwargs):
    """Initializes a database.
    
    Args
    ----

    Parameters
    ----------
    **kwargs :
        

    Returns
    -------

    """
    
    mod = sys.modules[__name__]
    
    setattr(mod, 'db', MoleculeDatabaseAggregator(**kwargs))

def get_db():
    """Returns the module's default database.
    Initializes the database with default paremeters if no database
    yet available.

    Parameters
    ----------

    Returns
    -------

    """
    
    mod = sys.modules[__name__]
    
    if not hasattr(mod, 'db'):
        
        init_db()
    
    return getattr(mod, 'db')

def lookup(m, tolerance = None):
    """

    Parameters
    ----------
    m :
        
    tolerance :
         (Default value = None)

    Returns
    -------

    """
    
    db = get_db()
    return db.lookup(m)

def adduct_lookup(
        mz,
        ionmode,
        adduct_constraints = True,
        tolerance = None,
    ):
    """Does a series of lookups in the database assuming various adducts.
    Calculates the exact mass for the m/z for each possible adduct
    and searches these exact masses in the database.
    
    Returns a dict of tuples with 3-3 numpy arrays.
    Keys of the dict are adduct types. The arrays are exact masses,
    database record details and accuracies (ppm).

    Parameters
    ----------
    mz :
        
    ionmode :
        
    adduct_constraints :
         (Default value = True)
    tolerance :
         (Default value = None)

    Returns
    -------

    """
    
    db = get_db()
    
    return db.adduct_lookup(
        mz,
        ionmode = ionmode,
        adduct_constraints = adduct_constraints,
        tolerance = tolerance,
    )

def adduct_lookup_many(
        mzs,
        adducts = None,
        ionmode = None,
        charge = None,
        adduct_constraints = True,
        tolerance = None,
    ):
    """Performs the lookup on a vector of m/z values.
    Calls the ``adduct_lookup`` method on all m/z's.
    
    Returns array of dicts with lookup results.

    Parameters
    ----------
    mzs :
        
    adducts :
         (Default value = None)
    ionmode :
         (Default value = None)
    charge :
         (Default value = None)
    adduct_constraints :
         (Default value = True)
    tolerance :
         (Default value = None)

    Returns
    -------

    """
    
    db = get_db()
    
    return db.adduct_lookup_many(
        mzs = mzs,
        adducts = adducts,
        ionmode = ionmode,
        charge = charge,
        adduct_constraints = adduct_constraints,
        tolerance = tolerance,
    )

def possible_classes(
        mz,
        ionmode,
        adduct_constraints = True,
        main_only = True,
        tolerance = None,
    ):
    """For a m/z returns a set of possible classes.

    Parameters
    ----------
    bool :
        main_only:
        Return only the set of main classes or subclasses.
    mz :
        
    ionmode :
        
    adduct_constraints :
         (Default value = True)
    main_only :
         (Default value = True)
    tolerance :
         (Default value = None)

    Returns
    -------

    """
    
    return set(
        r.hg.main if main_only else r.hg
        for rr in itertools.chain(
            adduct_lookup(
                mz,
                ionmode,
                adduct_constraints = adduct_constraints,
                tolerance = tolerance,
            ).values()
        )
        for r in rr[1]
        if r.hg is not None
    )


def records_string(
        records,
        adducts = None,
        databases = None,
        show_ppm = False,
        show_adduct = False,
        show_db = False,
        use_db_names = False,
    ):
    """

    Parameters
    ----------
    records :
        
    adducts :
         (Default value = None)
    databases :
         (Default value = None)
    show_ppm :
         (Default value = False)
    show_adduct :
         (Default value = False)
    show_db :
         (Default value = False)

    Returns
    -------

    """
    
    return MoleculeDatabaseAggregator.records_string(
        records      = records,
        adducts      = adducts,
        databases    = databases,
        show_ppm     = show_ppm,
        show_adduct  = show_adduct,
        show_db      = show_db,
        use_db_names = use_db_names
    )


def result_summary(result):
    """
    Provides a concise summary of the lookup results.
    
    Parameters
    ----------
    result : dict
        A result dictionary returned by a lookup.

    Returns
    -------
    A set of tuples with adduct type, database name, database ID and the
    string representation of the lipid.
    """
    
    return set(
        (
            adduct, # adduct
            lip.lab.db, # database name
            lip.lab.db_id, # database ID
            (
                lipproc.summary_str(lip.hg, lip.chainsum)
                if lip.hg and lip.chainsum else
                lip.hg.main
                if lip.hg else
                None
            ) # summary string or headgroup class or None
        )
        for adduct, res in result.items()
        for lip in res[1]
    )


def mass_from_name(
        name = None,
        database_preference = None,
        hg = None,
        chainsum = None,
        chains = (),
    ):
    
    db = get_db()
    
    return db.mass_from_name(
        name = name,
        database_preference = database_preference,
        hg = hg,
        chainsum = chainsum,
        chains = chains,
    )


def mz_from_name(
        adduct,
        name = None,
        database_preference = None,
        hg = None,
        chainsum = None,
        chains = (),
    ):
    
    db = get_db()
    
    return db.mz_from_name(
        adduct = adduct,
        name = name,
        database_preference = database_preference,
        hg = hg,
        chainsum = chainsum,
        chains = chains,
    )


def mz_lowest_error_from_name(
        measured_mz,
        adduct,
        name = None,
        **kwargs,
    ):
    """
    Regarding a measured m/z and an assumed adduct type and name
    returns the m/z of the record with matching name and lowest error.
    """
    
    db = get_db()
    
    return db.mz_lowest_error_from_name(
        measured_mz = measured_mz,
        adduct = adduct,
        name = name,
        **kwargs,
    )


