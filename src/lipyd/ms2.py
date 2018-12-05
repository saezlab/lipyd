#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import sys
import imp
import re
import math
import copy
import itertools
import collections
from argparse import Namespace
import numpy as np

from lipyd.common import *
import lipyd.mgf as mgf
import lipyd.mz as mzmod
import lipyd.session as session
import lipyd.settings as settings
import lipyd.lookup as lookup
import lipyd.fragdb as fragdb
import lipyd.moldb as moldb
import lipyd.lipproc as lipproc


ChainFragment = collections.namedtuple(
    'ChainFragment',
    ['c', 'u', 'fragtype', 'chaintype', 'i', 'intensity']
)


class MS2Identity(collections.namedtuple(
        'MS2IdentityBase',
        [
            'score',
            'max_score',
            'score_pct',
            'hg',
            'chainsum',
            'chains',
            'chain_details',
            'scan_details',
            'precursor_details',
        ]
    )):
    
    def __new__(
            cls,
            score = 0,
            max_score = 0,
            score_pct = 0,
            hg = None,
            chainsum = None,
            chains = None,
            chain_details = None,
            scan_details = None,
            precursor_details = None,
        ):
        
        return super(MS2Identity, cls).__new__(
            cls,
            score,
            max_score,
            score_pct,
            hg,
            chainsum = chainsum,
            chains = chains,
            chain_details = chain_details,
            scan_details = scan_details,
            precursor_details = precursor_details,
        )
    
    def __str__(self):
        
        return (
            lipproc.full_str(self.hg, self.chains)
                if self.chains else
            lipproc.summary_str(self.hg, self.chainsum)
        )
    
    def full_str(self):
        
        details = []
        
        details.append('score=%.01f' % self.score_pct)
        
        if self.scan_details:
            
            if self.scan_details.deltart is not None:
                details.append('deltart=%.02f' % self.scan_details.deltart)
            if self.scan_details.sample_id is not None:
                sample_id = self.scan_details.sample_id
                if isinstance(sample_id, tuple):
                    sample_id = ''.join(str(i) for i in sample_id)
                details.append('sample=%s' % sample_id)
            if self.scan_details.scan_id is not None:
                details.append('scan=%u' % self.scan_details.scan_id)
        
        return '%s[%s]' % (
            self.__str__(),
            ','.join(details),
        )
    
    def summary(self):
        
        return self.__str__(), self.score_pct
    
    def __eq__(self, other):
        
        return (
            isinstance(other, MS2Identity) and
            self.hg == other.hg and
            self.chainsum == other.chainsum and
            self.chains == other.chains
        )


ChainIdentificationDetails = collections.namedtuple(
    'ChainIdentificationDetails',
    ['rank', 'i', 'fragtype']
)
ChainIdentificationDetails.__new__.__defaults__ = (None, None, None)


ScanDetails = collections.namedtuple(
    'ScanDetails',
    ['sample_id', 'scan_id', 'source', 'deltart']
)
ChainIdentificationDetails.__new__.__defaults__ = (None, None, None, None)


PrecursorDetails = collections.namedtuple(
    'MS1LookupDetails',
    ['adduct', 'exmass', 'error', 'db', 'db_id', 'charge', 'id']
)
PrecursorDetails.__new__.__defaults__ = (None, None, None, None, None, None)


class mz_sorted(object):
    
    def __init__(self, scan):
        
        self.scan = scan
    
    def __enter__(self):
        
        self.scan.sort_mz()
    
    def __exit__(self, extyp, exval, tb):
        
        self.scan.sort_intensity()


class intensity_sorted(object):
    
    def __init__(self, scan):
        
        self.scan = scan
    
    def __enter__(self):
        
        self.scan.sort_intensity()
    
    def __exit__(self, extyp, exval, tb):
        
        self.scan.sort_mz()


class ScanBase(object):
    
    def __init__(
            self,
            mzs,
            ionmode,
            precursor = None,
            intensities = None,
            tolerance = None,
            scan_id = None,
        ):
        
        self.tolerance = tolerance or settings.get('ms2_tolerance')
        self.sorted_by = None
        self.mzs = mzs
        self.ionmode = ionmode
        self.adducts = {}
        self.intensities = (
            np.array([1.0] * len(self.mzs))
                if intensities is None else
            intensities
        )
        self.precursor = precursor
        self.scan_id = scan_id
        
        if self.mzs is not np.ndarray:
            
            self.mzs = np.array(self.mzs)
        
        if self.intensities is not np.ndarray:
            
            self.intensities = np.array(self.intensities)
        
        self.annotate()
        self.normalize_intensities()
        
        with mz_sorted(self):
            
            self.iisort = np.argsort(self.intensities)[::-1]
        
        self.irank = np.arange(len(self.mzs))
        self.imzsort  = np.argsort(self.mzs)
        self.sorted_by = 'intensities'
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def __len__(self):
        
        return len(self.mzs)
    
    def sort_mz(self):
        """
        Sorts the scan by m/z values ascending.
        """
        
        if self.sorted_by == 'mzs':
            
            return
            
        elif hasattr(self, 'imzsort') and self.sorted_by == 'intensities':
            
            self.sort(self.imzsort)
            
        else:
            
            isort = np.argsort(self.mzs)
            self.sort(isort)
        
        self.sorted_by = 'mzs'
    
    def sort_intensity(self):
        """
        Sorts the scan by intensities descending.
        """
        
        if self.sorted_by == 'intensities':
            
            return
            
        elif hasattr(self, 'iisort') and self.sorted_by == 'mzs':
            
            self.sort(self.iisort)
            
        else:
            
            isort = np.argsort(self.intensities)[::-1]
            self.sort(isort)
        
        self.sorted_by = 'intensities'
    
    def sort(self, isort):
        """
        Applies sorted indices to the scan.
        """
        
        self.intensities = self.intensities[isort]
        self.mzs = self.mzs[isort]
        
        for attr in ('irank', 'annot', 'inorm'):
            
            if hasattr(self, attr):
                
                setattr(self, attr, getattr(self, attr)[isort])
        
        for ad, data in iteritems(self.adducts):
            
            data['annot'] = data['annot'][isort]
    
    def annotate(self):
        """
        Annotates the fragments in the scan with identities provided by
        the fragment database.
        """
        
        self.annot = self.get_annot()
    
    def get_annot(self, precursor = None, tolerance = None):
        """
        Returns array of annotations.
        Makes it possible to use different precursor or tolerance.
        """
        
        precursor = precursor or self.precursor
        tolerance = tolerance or self.tolerance
        
        annotator = fragdb.FragmentAnnotator(
            self.mzs,
            self.ionmode,
            precursor,
            tolerance = tolerance,
        )
        
        return np.array(list(annotator)) # this is array
                                         # only to be sortable
    
    def normalize_intensities(self):
        """
        Creates a vector of normalized intensities i.e. divides intensities
        by their maximum.
        """
        
        self.imax  = self.intensities.max()
        self.inorm = self.intensities / self.imax


class Scan(ScanBase):
    
    method_hg = {
        'fa_neg_1': ('FA', ()),
    }
    
    def __init__(
            self,
            mzs,
            ionmode,
            precursor = None,
            intensities = None,
            ms1_records = None,
            add_precursor_details = False,
            scan_id = None,
            sample_id = None,
            source = None,
            deltart = None,
            logger = None,
            verbose = False,
            tolerance = None,
            ms1_tolerance = None,
            rt = None,
        ):
        
        ScanBase.__init__(
            self,
            mzs,
            ionmode,
            precursor,
            intensities,
            tolerance = tolerance,
        )
        
        # get some settings
        self.ms1_tolerance = ms1_tolerance or settings.get('ms1_tolerance')
        self.check_ratio_g = settings.get(
            'even_chain_fragment_intensity_ratios_gl_gpl'
        )
        self.check_ratio_s = settings.get(
            'even_chain_fragment_intensity_ratios_sl'
        )
        self.iratio_logbase = settings.get(
            'chain_fragment_instensity_ratios_logbase'
        )
        self.chain_details = settings.get('ms2_scan_chain_details')
        
        if ms1_records is None and precursor is not None:
            
            # do the database lookup if not provided,
            # this is not efficient but makes possible
            # to easily use standalone `Scan` instances
            # for testing and such
            self.ms1_records = moldb.adduct_lookup(
                precursor, ionmode, tolerance = self.ms1_tolerance
            )
            
        else:
            
            # even if precursor is None, we end up with an empty dict
            self.ms1_records = ms1_records or {}
        
        self.add_precursor_details = add_precursor_details
        
        self.scan_id   = scan_id
        self.sample_id = sample_id
        self.source    = source
        self.deltart   = deltart
        self.rt        = rt
        self.log       = logger
        self.verbose   = verbose
        
        self.scan_details = ScanDetails(
            sample_id = self.sample_id,
            scan_id   = self.scan_id,
            source    = self.source,
            deltart   = self.deltart,
        )
    
    @classmethod
    def from_mgf(
            cls,
            fname,
            scan_id,
            ionmode,
            sample_id = None,
            precursor = None,
            mgf_charge = None,
            **kwargs
        ):
        
        mgfreader = mgf.MgfReader(fname, charge = mgf_charge)
        sc = mgfreader.scan_by_id(scan_id)
        
        precursor = precursor or mgfreader.precursor_by_id(scan_id)
        
        if sc is not None:
            
            return cls(
                mzs = sc[:,0],
                intensities = sc[:,1],
                ionmode = ionmode,
                precursor = precursor,
                scan_id = scan_id,
                sample_id = sample_id,
                source = fname,
                **kwargs
            )
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def print_scan(self):
        """
        Prints the list of fragments as an annotated table.
        """
        
        if self.log:
            
            self.log.msg(self.scan_str())
    
    def show(self):
        """
        Prints the scan table to standard output.
        """
        
        sys.stdout.write(self.scan_str())
    
    def scan_str(self):
        """
        Returns the scan table as string.
        """
        
        lindent = ' ' * 12
        
        header = '%s\n%s%s\n' % (
            ''.join((
                lindent,
                'Frag. m/z'.ljust(12),
                'Intensity'.ljust(12),
                'Identity'.ljust(36),
                'NL mass'.rjust(12)
            )),
            lindent,
            '=' * 12 * 6
        )
        
        table = '\n'.join((
            ''.join((
                lindent,
                '%12.4f' % self.mz[i],
                '%10u'   % self.intensities[i],
                ann.name,
                (
                    '%12.4f' % self.nl(mz[i])
                        if self.precursor else
                    'NA'.rjust(12)
                )
            ))
            for i in xrange(len(self.mz))
            for ann in (
                self.annot[i]
                    if self.annot[i] else
                (Namespace(name = 'Unknown'),)
            )
        ))
        
        return '%s\n%s\n\n' % (
            self.sample.__str__(),
            header,
            table
        )
    
    def html_table(self):
        
        # TODO
        pass
    
    def nl(self, mz, adduct = None):
        """
        For m/z returns the corresponding neutral loss m/z.
        If precursor ion mass is unknown returns `numpy.nan`.
        """
        
        if adduct is None:
            
            return self.precursor - mz if self.precursor else np.nan
            
        else:
            
            self.adduct(adduct)
            
            return self.adducts[adduct]['fake_precursor'] - mz
    
    def full_list_str(self):
        """
        Returns list of fragments as single string.
        """
        
        return '; '.join(
            '/'.join(
                '%s (%u)' % (
                    ann.name,
                    self.intensities[i]
                )
                    if ann is not None else
                'Unknown (%.03f) (%u)' % (
                    self.mzs[i],
                    self.intensities[i]
                )
            )
            for i in xrange(len(self))
            for ann in (
                self.annot[i]
                    if self.annot[i] else
                (None,)
            )
        )
    
    def most_abundant_mz(self):
        """
        Returns the m/z of the fragment with highest intensity.
        """
        
        result = self.mzs[0]
        
        if self.verbose:
            
            self.log.msg('\t\t  -- Most abundant m/z is %.03f' % result)
        
        return result
    
    def mz_match(self, mz_detected, mz):
        """
        Tests if two m/z's can be accepted to be equal.
        """
        
        return lookup.match(mz_detected, mz, self.tolerance)
    
    def mz_lookup(self, mz):
        """
        Returns the index of the closest m/z value
        detected in the scan if it is within the
        range of tolerance, otherwise None.
        """
        
        self.sort_mz()
        
        imz = lookup.find(self.mzs, mz, self.tolerance)
        i = self.imzsort[imz] if imz else None
        
        self.sort_intensity()
        
        return i
    
    def has_mz(self, mz):
        """
        Tells if an m/z exists in this scan.
        """
        
        result = self.mz_lookup(mz) is not None
        
        if self.verbose:
            
            self.log.msg(
                '\t\t  -- m/z %.03f occures in this scan? -- %s' % (
                    mz, str(result)
                )
            )
        
        return result
    
    def has_nl(self, nl, adduct = None):
        """
        Tells if a neutral loss exists in this scan.
        """
        
        result = self.has_mz(self.nl(nl, adduct = adduct))
        
        if self.verbose:
            
            self.feature.msg(
                '\t\t  -- neutral loss of %.03f occures in '
                'this scan? Looked up m/z %.03f - %.03f = %.03f -- %s' % (
                    nl,
                    self.precursor,
                    nl,
                    self.nl(nl, adduct = adduct),
                    str(result)
                )
            )
        
        return result
    
    def fragment_by_name(self, name, adduct = None):
        """
        Returns the index of a fragment by its name.
        Returns `None` if the fragment does not exist in the scan.
        Returns `False` if the fragment name could not be found in
        the database.
        
        The lookup still goes by m/z, the name first looked up in the
        fragment database and the scan searched for the corresponding m/z.
        The name makes if obvious if this is a charged fragment or a neutral
        loss, hence it is not necessary to provide this information.
        
        Args
        ----
        :param str name:
            Fragment full name as used in the database 2nd column,
            e.g. `PE [P+E] (140.0118)`.
        """
        
        frag = fragdb.by_name(name, self.ionmode)
        
        if frag is not None:
            
            if frag[6] == 0:
                
                return self.nl_lookup(frag[0], adduct = adduct)
                
            else:
                
                return self.mz_lookup(frag[0])
        
        return False
    
    def has_fragment(self, name, adduct = None):
        """
        Tells if a fragment exists in this scan by its name.
        
        Returns bool or `None` if fragment name could not be found
        in the database.
        """
        
        i = self.fragment_by_name(name, adduct = adduct)
        
        return None if i is False else i is not None
    
    def nl_lookup(self, nl, adduct = None):
        """
        Looks up if a neutral loss exists in this scan and returns its index.
        """
        
        return self.mz_lookup(self.nl(nl, adduct = adduct))
    
    def most_abundant_fragment_is(self, name, adduct = None):
        """
        Tells if the fragment name is the highest abundant.
        Returns `None` if the fragment name could not be
        found in the database.
        """
        
        frag = fragdb.by_name(name, self.ionmode)
        
        if frag is not None:
            
            mz = (
                self.nl(frag[0], adduct = adduct)
                    if frag[6] == 0 else
                frag[0]
            )
            
            return self.mz_match(self.mzs[0], mz)
    
    def fragment_among_most_abundant(self, name, n = 2, adduct = None):
        """
        Tells if the fragment is among the top `n`.
        """
        
        frag = fragdb.by_name(name, self.ionmode)
        
        if frag is not None:
            
            mz = (
                self.nl(frag[0], adduct = adduct)
                    if frag[6] == 0 else
                frag[0]
            )
            
            return self.mz_among_most_abundant(mz, n = n)
    
    def fragment_percent_of_most_abundant(
            self,
            name,
            percent = 80.0,
            adduct = None,
        ):
        """
        Tells if a fragment has at least certain percent of intensity
        compared to the highest peak.
        """
        
        frag = fragdb.by_name(name, self.ionmode)
        
        if frag is not None:
            
            mz = (
                self.nl(frag[0], adduct = adduct)
                    if frag[6] == 0 else
                frag[0]
            )
            
            return self.mz_percent_of_most_abundant(mz, percent = percent)
    
    def most_abundant_mz_is(self, mz):
        """
        Tells if the m/z with the highest intensity is `mz`.
        Returns `None` if the fragment name could not be
        found in the database.
        """
        
        result = self.mz_match(self.most_abundant_mz(), mz)
        
        if self.verbose:
            
            self.log.msg(
                '\t\t  -- Is m/z %.03f the most abundant one? -- %s' % (
                    mz,
                    str(result)
                )
            )
        
        return result
    
    def mz_among_most_abundant(self, mz, n = 2):
        """
        Tells if an m/z is among the most aboundant `n` fragments
        in a spectrum.
        
        Args
        ----
        :param float mz:
            The m/z value.
        :param int n:
            The number of most abundant fragments considered.
        """
        
        self.sort_mz()
        
        i = lookup.find(
            self.mzs[self.irank < n], # intensity rank < n
            mz,
            self.tolerance
        )
        
        self.sort_intensity()
        
        if self.verbose:
            
            self.log.msg(
                '\t\t  -- m/z %.03f is among the %u most abundant? -- %s' % (
                    mz, n, str(i is not None)
                )
            )
        
        return i is not None
    
    def nl_among_most_abundant(self, nl, n = 2, adduct = None):
        """
        Tells if a neutral loss corresponds to one of the
        most aboundant `n` fragments in a spectrum.
        
        Args
        ----
        :param float nl:
            The mass of the neutral loss.
        :param int n:
            The number of most abundant fragments considered.
        """
        
        result = self.mz_among_most_abundant(
            self.nl(nl, adduct = adduct),
            n = n,
        )
        
        if self.verbose:
            
            self.log.msg(
                '\t\t  -- neutral loss %.03f is among '
                'the %u most abundant? -- %s' % (
                    nl, n, str(result)
                )
            )
        
        return result
    
    def get_intensity(self, mz):
        """
        Returns the relative intensity of a fragment ion from its m/z.
        Value is `None` if m/z does not present.
        """
        
        i = self.mz_lookup(mz)
        
        if i is not None:
            
            return self.inorm[i]
        
        return None
    
    def get_nl_intensity(self, nl, adduct = None):
        """
        Returns the relative intensity of a neutral loss fragment ion.
        Value is `None` if neutral loss does not present.
        """
        
        return self.get_intensity(self.nl(nl, adduct = adduct))
    
    def mz_percent_of_most_abundant(self, mz, percent = 80.0):
        """
        Tells if an m/z has at least certain percent of intensity
        compared to the most intensive fragment.
        
        Args
        ----
        :param float mz:
            The m/z value.
        :param float percent:
            The threshold in percent of the highest intensity.
        
        """
        
        i = self.get_intensity(mz)
        result = i and i >= percent / 100.
        
        if self.verbose:
            
            self.feature.msg(
                '\t\t  -- m/z %.03f has abundance at least %.01f %% of'
                ' the highest abundance? -- %s\n' % (
                    mz, percent, str(result)
                )
            )
        
        return result
    
    @classmethod
    def match_chtype(cls, value, accepted):
        """
        Matches strings or strings to set of strings, optionally negative.
        Calls `match_chattr` with `typ = basestring`.
        """
        
        return cls.match_chattr(value, accepted, typ = basestring)
    
    @staticmethod
    def match_chattr(value, accepted, typ = int):
        """
        Args
        ----
        :param int value:
            The actual value.
        :param int,set accepted:
            A single value or a set of values to match against.
            Negative match is possible by providing a tuple with `False`
            as it's first element and the set of not acceptable values
            as the second element.
        """
        
        return (
            accepted is None or
            # simple match
            (isinstance(accepted, typ) and value == accepted) or
            (
                not isinstance(accepted, typ) and (
                    # multiple values
                    value in accepted or (
                        # negation
                        hasattr(accepted, '__getitem__') and
                        accepted[0] == False and
                        value not in accepted[1]
                    )
                )
            )
        )
    
    @classmethod
    def match_annot(
            cls,
            annot,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None
        ):
        """
        Tests a fragment annotation against criteria of fragment type,
        chain type, carbon count and unsaturation.
        """
        
        return all((
            cls.match_chattr(annot.fragtype, frag_type, typ = basestring),
            cls.match_chattr(annot.chaintype, chain_type, typ = basestring),
            cls.match_chattr(annot.c, c),
            cls.match_chattr(annot.u, u),
        ))
    
    def highest_fragment_by_chain_type(
            self,
            head = None,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Returns the highest instensity fragment matching a particular
        chain type.
        Returns fragment index or `None` if no such fragment exists.
        
        Arguments passed to `chain_fragment_type_is`.
        """
        
        frags = self.fragments_by_chain_type(
            head = head,
            frag_type = frag_type,
            chain_type = chain_type,
            c = c,
            u = u,
            adduct = adduct,
        )
        
        try:
            
            return next(frags)
            
        except StopIteration:
            
            return None
    
    def fragments_by_chain_type(
            self,
            head = None,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Collects fragments matching a particular chain type.
        Yields indices.
        Arguments passed to `chain_fragment_type_is`.
        """
        
        head = len(self.mzs) if head is None else min(head, len(self.mzs))
        
        for i in xrange(head):
            
            if self.chain_fragment_type_is(
                i,
                frag_type = frag_type,
                chain_type = chain_type,
                c = c,
                u = u,
                return_annot = False,
                adduct = adduct,
            ):
                
                yield i
    
    def chain_fragment_type_among_most_abundant(
            self,
            n = 2,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Tells if a particular type of aliphatic chain fragment can be
        found among the `n` highest intensity fragments.
        
        Arguments passed to `chain_fragment_type_is`.
        """
        
        return bool(len(list(
            self.fragments_by_chain_type(
                head = n,
                frag_type = frag_type,
                chain_type = chain_type,
                c = c,
                u = u,
                adduct = adduct,
            )
        )))
    
    def chain_fragment_type_is(
            self,
            i,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            return_annot = False,
            adduct = None,
        ):
        """
        Tells if an aliphatic chain fragment is a specified type. The type
        should be the string representation of the fragment,
        e.g. `FA-O` for fatty acid minus oxygen fragments.
        
        Returns bool or fragment annotations if `return_annot = True`.
        
        Args
        ----
        :param int i:
            Index of the fragment.
        :param bool return_annot:
            Return iterator with the matching fragment annotations.
        """
        
        if i >= len(self.mzs):
            
            return False
        
        annot = self.annot if adduct is None else self.adduct_annot(adduct)
        
        result = any((
            self.match_annot(an, frag_type, chain_type, c, u)
            for an in annot[i]
        ))
        
        if self.verbose:
            
            criteria = []
            if frag_type is not None:
                criteria.append('of type `%s`' % frag_type)
            if chain_type is not None:
                criteria.append('of chain type `%s`' % chain_type)
            if c is not None:
                criteria.append('with carbon count of %a' % c)
            if u is not None:
                criteria.append('with unsaturation of %a' % u)
            
            self.log.msg(
                '\t\t  -- Fragment #%u (%.03f): '
                'is it a fragment %s? -- %s' % (
                    i,
                    self.mz[i],
                    ' and '.join(criteria),
                    str(result)
                )
            )
        
        if return_annot:
            
            result = (
                an
                for an in annot[i]
                if self.match_annot(an, frag_type, chain_type, c, u)
            )
        
        return result
    
    def chains_of_type(
            self,
            chain_type = None,
            frag_type = None,
            c = None,
            u = None,
            yield_annot = False,
            adduct = None,
        ):
        """
        Iterates chain fragments matching certain criteria.
        Yields fragment indices or indices with annotations.
        
        Args
        ----
        :param bool yield_annot:
            Yield tuples of indices and annotations instead of indices only.
        """
        
        for i in xrange(len(self.mzs)):
            
            if self.chain_fragment_type_is(
                i = i,
                chain_type = chain_type,
                frag_type = frag_type,
                c = c,
                u = u,
                adduct = adduct,
            ):
                
                if yield_annot:
                    
                    for annot in self.chain_fragment_type_is(
                        i = i,
                        chain_type = chain_type,
                        frag_type = frag_type,
                        c = c,
                        u = u,
                        return_annot = True,
                        adduct = adduct,
                    ):
                        
                        yield i, annot
                    
                else:
                    
                    yield i
    
    def has_chain_fragment_type(
            self,
            chain_type = None,
            frag_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Tells if at least one fragment matches certain criteria.
        Arguments passed to `chain_fragment_type_is`.
        """
        
        return self.highest_fragment_by_chain_type(
            chain_type = chain_type,
            frag_type = frag_type,
            c = c,
            u = u,
            adduct = adduct,
        ) is not None
    
    def matching_chain_combinations(
            self,
            record,
            head = None,
            intensity_threshold = None,
            expected_intensities = None,
            no_intensity_check = False,
            chain_param = (),
            adduct = None,
        ):
        """
        Provides a way to see if specific chain combinations exist.
        The database record defines the chain layout of the molecule.
        Further arguments are passed to `chain_combinations`.
        The `chain_param` tuple contains dicts to match against chain
        fragments. All of these dicts must match at least one fragment
        identification. Only combinations matching all criteria yielded.
        
        Args
        ----
        :param lipproc.LipidRecord record:
            A lipid database record matching the MS1 m/z.
        :param int head:
            Consider only the n most intensive fragments.
        :param float intensity_threshold:
            Consider only fragments with intensity higher than threshold.
            Relative to highest fragment, between 0 and 1.
        :param expected_intensities:
            See at `intensity_ratios`.
        :param bool no_intensity_check:
            Completely skip checking intensity ratios.
        :param tuple chain_param:
            Tuple of dicts. Each dict contains criteria for one chain moiety.
            Keys can be `chain_type`, `frag_type`, `c` and `u`.
            These can be single str or int values or sets of multiple
            values. If omitted or `None` any value will pass the filter.
            An empty tuple which is the default value will pass through
            everything, this is equivalent with calling `chain_combinations`.
        """
        
        def match(key, param, value):
            
            return (
                key not in param or
                param[key] is None or (
                    type(param[key]) in {int, str} and
                    value == param[key]
                ) or (
                    type(param[key]) in {set, list, tuple} and
                    value in param[key]
                )
            )
        
        if (
            record.chainsum and
            len(record.chainsum) > 1 and
            len(chain_param) == 1
        ):
            
            chain_param = chain_param * len(record.chainsum)
        
        for chains, details in self.chain_combinations(
            record,
            head = None,
            intensity_threshold = 0,
            expected_intensities = None,
            no_intensity_check = False,
            frag_types = None,
            fragment_details = True,
            adduct = adduct,
        ):
            
            if (
                not chain_param or
                all((
                    not param or
                    any((
                        all((
                            match('chain_type', param, ch.typ),
                            match('frag_type', param, details.fragtype[i]),
                            match('c', param, ch.c),
                            match('u', param, ch.u),
                        ))
                        for i, ch in enumerate(chains)
                    ))
                    for param in chain_param
                ))
            ):
                
                yield chains, details
    
    def has_chain_combination(
            self,
            record,
            head = None,
            intensity_threshold = None,
            expected_intensities = None,
            no_intensity_check = False,
            chain_param = (),
            adduct = None,
        ):
        """
        Tells if a certain chain combination exists in the scan.
        
        Arguments passed to `matching_chain_combinations`.
        """
        
        ccomb = self.matching_chain_combinations(
            record = record,
            head = head,
            intensity_threshold = intensity_threshold,
            expected_intensities = expected_intensities,
            no_intensity_check = no_intensity_check,
            chain_param = chain_param,
            adduct = adduct,
        )
        
        try:
            
            _ = next(ccomb)
            
            return True
            
        except StopIteration:
            
            return False
    
    def _matching_chain_pairs(
            self,
            record,
            chain_type = None,
            frag_type = None,
            c = None,
            u = None,
            partner_chain_types = None,
            partner_frag_types = None,
            count_only = False,
            adduct = None,
        ):
        
        # small caching of constraint matching
        type_pos = {}
        
        def get_positions(self, frag_type):
            
            if frag_type not in type_pos:
                
                type_pos[frag_type] = self.positions_for_frag_type(
                    record, frag_type
                )
            
            return type_pos[frag_type]
        # ##
        
        for i, iannot in self.chains_of_type(
            chain_type = chain_type,
            frag_type = frag_type,
            c = c,
            u = u,
            yield_annot = True,
            adduct = adduct,
        ):
            
            partner_c = record.chainsum.c - annot.c
            partner_u = record.chainsum.u - annot.u
            
            if partner_c < 1 or partner_u < 0:
                
                continue
            
            pos_i = get_positions(iannot.fragtype)
            
            for j, jannot in self.chains_of_type(
                c = partner_c,
                u = partner_u,
                yield_annot = True,
                adduct = adduct,
            ):
                
                if (
                    partner_chain_types is None or
                    jannot.chaintype in partner_chain_types
                ) and (
                    partner_frag_types is None or
                    jannot.fragtype in partner_frag_types
                ):
                    
                    pos_j = get_positions(jannot.fragtype)
                    
                    if (
                        not pos_i or
                        not pos_j or (
                            len(pos_i) == 1 and
                            len(pos_j) == 1 and
                            not pos_i - pos_j
                        )
                    ):
                        
                        continue
                    
                    yield (
                        lipproc.Chain(
                            
                        )
                    )
    
    def positions_for_frag_type(self, record, frag_type):
        """
        Returns the possible chain positions for a record and a fragment type.
        """
        
        # constraints for the fragment type
        constr = fragdb.constraints(frag_type, self.ionmode)
        # set of possible positions of the chain
        # which this fragment originates from
        return lipproc.match_constraints(record, constr)[1]
    
    def is_chain(self, i, adduct = None):
        """
        Examines if a fragment has an aliphatic chain.
        """
        
        annot = self.adduct_annot(adduct)
        
        result = any(not np.isnan(an.c) for an in annot[i])
        
        if self.verbose:
            
            self.log.msg(
                '\t\t -- Fragment #%u (%.03f)'
                'has an aliphatic chain? -- %s' % (
                    i,
                    self.mzs[i],
                    str(result)
                )
            )
        
        return result
    
    def is_chain_type(self, i, typ = 'FA', adduct = None):
        """
        Checks if a fragment might origin from a certain aliphatic
        chain type (e.g. `FA` -- fatty acyl, `FAL` -- fatty alkyl,
        `Sph` -- sphingosin base).
        """
        
        return self.chain_fragment_type_is(
            i, chain_type = typ, adduct = adduct
        )
    
    def is_fa(self, i, adduct = None):
        """
        Tells if a fragment origins from a fatty acyl moiety.
        """
        
        return self.is_chain_type(i, adduct = adduct)
    
    def is_fal(self, i, adduct = None):
        """
        Tells if a fragment origins from a fatty alkyl moiety.
        """
        
        return self.is_chain_type(i, 'FAL', adduct = adduct)
    
    def is_sph(self, i, adduct = None):
        """
        Tells if a fragment origins from a shpingosin backbone.
        """
        
        return self.is_chain_type(i, 'Sph', adduct = adduct)
    
    def is_type(self, i, typ, adduct = None):
        """
        Tells if a fragment is a certain type.
        """
        
        return self.chain_fragment_type_is(
            i, frag_type = typ, adduct = adduct
        )
    
    def annot_by_type(
            self,
            i,
            chain_type = None,
            frag_type = None,
            adduct = None,
        ):
        """
        Returns the annotations matching certain types.
        """
        
        annot = self.adduct_annot(adduct)
        
        return tuple(
            an
            for an in annot[i]
            if (
                self.match_chtype(an.chaintype, chain_type) and
                self.match_chtype(an.fragtype,  frag_type)
            )
        )
    
    def cu_by_type(
            self,
            i,
            chain_type = None,
            frag_type = None,
            adduct = None,
        ):
        """
        Returns `(carbon count, unsaturation)` tuples for fragment `i`
        considering only the the requested chain types and fragment types.
        """
        
        return tuple(
            (a.c, a.u)
            for a in
            self.annot_by_type(
                i,
                chain_type = chain_type,
                frag_type = frag_type,
                adduct = adduct,
            )
        )
    
    def _build_chain_list(self, annot = None):
        """
        Builds a list of chains which facilitates the anlysis of chain
        combinations.
        """
        
        annot = annot if type(annot) is np.ndarray else self.annot
        
        return tuple(
            ChainFragment(
                a.c, a.u, a.fragtype, a.chaintype, i, self.intensities[i]
            )
            for i, aa in enumerate(annot)
            for a in aa
            if a.c and not np.isnan(a.c)
        )
    
    def build_chain_list(self, rebuild = False):
        
        if (
            not rebuild and
            hasattr(self, 'chain_list')
        ):
            return
        
        self.chain_list = self._build_chain_list()
    
    def chain_among_most_abundant(
            self,
            head = 1,
            chain_type = None,
            frag_type = None,
            c = None,
            u = None,
            min_mass = None,
            skip_non_chains = False,
            adduct = None,
        ):
        """
        Returns `True` if the defined type of chain fragment can be found
        among the most abundant fragments.
        """
        
        if self.verbose:
            
            self.log.msg(
                '\t\t -- Checking for certain type of chain among the top '
                '%u fragments.' % head
            )
        
        result = any((
            self.chain_fragment_type_is(
                i,
                frag_type = frag_type,
                chain_type = chain_type,
                c = c,
                u = u,
                adduct = adduct,
            )
            for i in (
                xrange(head)
                    if not skip_non_chains else
                itertools.islice(
                    (
                        i for i in xrange(len(self.mzs))
                            if (
                                not skip_non_chains or self.is_chain(i)
                            ) and (
                                min_mass is None or self.mzs[i] >= min_mass
                            )
                    ),
                    head
                )
            )
        ))
        
        if self.verbose:
            
            self.log.msg(
                '\t\t -- Checked certain type of chain among the top '
                '%u fragments. -- %s' % (head, str(result))
            )
        
        return result
    
    def get_most_abundant_chain(
            self,
            head = 1,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Looks up the most abundant fatty acid fragment of the given type.
        Returns the fragment index.
        """
        
        for i in xrange(len(self)):
            
            if self.chain_fragment_type_is(
                i,
                frag_type = frag_type,
                chain_type = chain_type,
                c = c,
                u = u,
                adduct = adduct,
            ):
                
                return i
    
    def chain_percent_of_most_abundant(
            self,
            percent,
            frag_type = None,
            chain_type = None,
            c = None,
            u = None,
            adduct = None,
        ):
        """
        Tells if a certain chain present with an abundance at least the
        given percent of the most abundant fragment.
        
        Args
        ----
        :param float percent:
            Percentage, between 0 and 100.
        """
        
        result = any((
            self.chain_among_most_abundant(
                i,
                frag_type = frag_type,
                chain_type = chain_type,
                c = c,
                u = u,
                adduct = adduct,
            )
            for i in
            itertools.takewhile(
                lambda i:
                    self.inorm[i] > percent / 100.0,
                xrange(len(self.mzs))
            )
        ))
        
        return result
    
    def mz_most_abundant_fold(self, mz, fold):
        """
        Tells if an m/z is the most abundant fragment
        and it has at least a certain fold higher intensity
        than any other fragment.
        
        Args
        ----
        :param float mz:
            The m/z value.
        :param float fold:
            The m/z must be this times higher than any other.
        """
        
        result = (
            self.most_abundant_mz_is(mz) and (
                len(self.mzs) == 1 or
                self.intensities[1] * fold <= self.intensities[0]
            )
        )
        
        if self.verbose:
            
            self.log.msg(
                '\t\t  -- m/z %.03f is at least %u times higher than '
                'any other? -- %s\n' % (mz, fold, str(result))
            )
        
        return result
    
    def cer_fa_test(self, i_fa, i_sph, adduct = None):
        
        return (
            self.chain_fragment_type_is(
                i_fa,
                frag_type = 'FA+C2+NH2-O',
                adduct = adduct,
            ) and
            self.chain_fragment_type_id(
                i_sph,
                frag_type = 'Sph-C2H4-NH2-H2O',
                adduct = adduct,
            ) and
            self.intensities[i_fa] > self.intensities[i_sph] * 2
        )
    
    def has_chain_combinations(self, rec, adduct = None, **kwargs):
        """
        Calls `chain_combinations` only to check if at least one
        conbination explicitely confirmed.
        """
        
        ccomb = self.chain_combinations(rec, adduct = adduct, **kwargs)
        
        try:
            
            _ = next(ccomb)
            
            return True
            
        except StopIteration:
            
            return False
    
    def chain_combinations(
            self,
            rec,
            head = None,
            intensity_threshold = 0,
            expected_intensities = None,
            no_intensity_check = False,
            frag_types = None,
            fragment_details = None,
            adduct = None,
        ):
        """
        Finds all combinations of chain derived fragments matching the
        total carbon count and unsaturation in a database record.
        
        Yields tuple of chains (`lipproc.Chain` objects).
        
        Arguments not listed here explained at `frags_for_positions`.
        
        Args
        ----
        :param lipproc.LipidRecord rec:
            The database record to match against.
        :param bool no_intensity_check:
            Completely skip checking intensity ratios.
        :param float intensity_threshold:
            Only fragments with intensities above this threshold will be
            considered. Intensities relative to the highest, between 0 and 1.
        :param tuple frag_types:
            See at `frags_for_positions`.
        """
        
        if not rec.chainsum and not rec.chains:
            
            return
        
        self.build_chain_list()
        
        chainsum = rec.chainsum or lipproc.sum_chains(rec.chains)
        
        frags_for_position = self.frags_for_positions(
            rec,
            head = head,
            intensity_threshold = intensity_threshold,
            frag_types = frag_types,
            adduct = adduct,
        )
        
        if len(frags_for_position) != len(rec.chainsum.typ):
            # if one or more chains have no corresponding fragment
            # we do not yield anything;
            # for finding missing those chains `missing_chains`
            # can be used
            return
        
        # iterate all combinations
        for frag_comb in itertools.product(
            *(
                # making a sorted list of lists from the dict
                i[1] for i in
                sorted(frags_for_position.items(), key = lambda i: i[0])
            )
        ):
            
            if (
                sum(frag.c for frag in frag_comb) == chainsum.c and
                sum(frag.u for frag in frag_comb) == chainsum.u
            ):
                
                if (
                    # bypass intensity check
                    no_intensity_check or
                    self._intensity_check(
                        frag_comb, chainsum, expected_intensities
                    )
                ):
                    
                    # now all conditions satisfied:
                    yield self._chains_frag_comb(
                        frag_comb, chainsum, details = fragment_details
                    )
    
    def frags_for_positions(
            self,
            rec,
            head = None,
            intensity_threshold = 0,
            frag_types = None,
            adduct = None,
        ):
        """
        Returns the possible fragments for each positions (sn1, sn2 in
        glycerophospholipids, sphingosine base and N-acyl in sphingolipids,
        etc).
        
        :param int head:
            If `None` or `numpy.inf` all fragment ions will be considered,
            otherwise only the first most aboundant until the number `head`.
        :param float intensity_threshold:
            Only fragments with intensities above this threshold will be
            considered. Intensities relative to the highest, between 0 and 1.
        :param tuple frag_types:
            Limit the query to certain fragment types in addition to
            built in fragment constraints and other criteria.
            A tuple of tuples with fragment type names can be provided
            each for one position with None values where default fragment
            types should be used. E.g. `(('FA_mH', 'Lyso_PA'), None)` means
            the chain in first position might be found as fatty acid minus
            hydrogen fragment or lysophosphatidic acid fragment, while the
            second position could be anything allowed by the built in
            constraints.
        """
        
        frags_for_position = collections.defaultdict(list)
        
        chain_list = self.adduct_chain_list(adduct)
        
        for frag in chain_list:
            
            if (
                (head and frag.i >= head) or
                self.inorm[frag.i] < intensity_threshold
            ):
                
                break
            
            chpos = self.positions_for_frag_type(rec, frag.fragtype)
            
            for ci in chpos:
                
                if (
                    # frag_types constraints
                    not frag_types or
                    not frag_types[ci] or
                    frag.fragtype in frag_types[ci]
                ):
                    
                    frags_for_position[ci].append(frag)
        
        return dict(frags_for_position)
    
    def intensity_ratios(
            self,
            intensities,
            frag_indices = None,
            expected = None,
            logbase = None,
        ):
        """
        Tells if the ratio of a list of intensities fits the one in
        `expected` or is more or less even if `expected` is `None`.
        
        :param list intensities:
            List of intensities.
        :param list expected:
            List with expected intensity proportions. E.g. `[1, 1, 2]`
            means the third ion is twice higher intense than the 2 others.
        :param int logbase:
            The fold difference tolerance when comparing intensities.
            E.g. if this is 2, then an almost twice less or more intense
            ion will considered to have similar intensity.
        """
        
        logbase = settings.get('chain_fragment_instensity_ratios_logbase')
        
        if len(intensities) == 1:
            
            return True
        
        if any(i <= 0.0 for i in intensities):
            
            raise ValueError(
                'Negative intensity value encountered'
                '(sample=%s, ion mode=%s, scan=%u)' % (
                    str(self.sample), self.ionmode, self.scan_id
                )
            )
        
        frag_indices = frag_indices or list(range(len(intensities)))
        # to know if one fragment contributes more than one times;
        # intensities divided by the times the fragment is incident
        cntr = collections.Counter(frag_indices)
        
        # by default expecting more or less equal intensities
        if expected is None:
            
            expected = [1.0] * len(intensities)
        
        # intensities corrected by the expected and the counts
        intcorr = [
            ins / expected[i] / cntr[ind]
            for (i, ins), ind in zip(enumerate(intensities), frag_indices)
        ]
        
        return (
            all((
                math.log(co[0], logbase) - math.log(co[1], logbase) <= 1
                for co in itertools.combinations(intcorr, 2)
            ))
        )
    
    def _intensity_check(
            self,
            frag_comb,
            chainsum,
            expected_intensities = None
        ):
        """
        Performs the chain intensity ratio check according to settings.
        """
        
        return (
            not (
                # need to check intensity ratios
                (chainsum.typ[0] == 'Sph' and self.check_ratio_s) or
                (chainsum.typ[0] != 'Sph' and self.check_ratio_g) or
                expected_intensities
            ) or
            self.intensity_ratios(
                # intensity ratios are ok
                intensities = tuple(f.intensity for f in frag_comb),
                frag_indices = tuple(f.i for f in frag_comb),
                expected = expected_intensities,
                logbase = self.iratio_logbase
            )
        )
    
    def _chains_frag_comb(
            self,
            frag_comb,
            chainsum,
            details = None,
            missing_position = None,
            missing_chain = None,
        ):
        """
        Returns a tuple of chains from a fragment annotation combination
        and a database record chain summary object.
        
        Potentially includes a missing chain which does not yield any
        fragment.
        """
        
        # boolean: whether we provide details or not
        details = self.chain_details if details is None else details
        
        return (
            tuple(
                lipproc.Chain(
                    c = frag_comb[ifrag].c,
                    u = frag_comb[ifrag].u,
                    typ = frag_comb[ifrag].chaintype,
                    attr = lipproc.ChainAttr(
                        # take the sphingosine base type
                        # from the chainsum of the record
                        sph = chainsum.attr[ichain].sph,
                        ether = frag_comb[ifrag].chaintype == 'FAL',
                        oh = chainsum.attr[ichain].oh
                    )
                )
                if ifrag is not None else
                missing_chain
                # chain indices and fragment indices
                for ichain, ifrag in iterator_insert(
                    len(chainsum),
                    missing_position,
                )
            ),
            ChainIdentificationDetails(
                rank     = tuple(
                    frag_comb[ifrag].i
                    if ifrag is not None else None
                    for ichain, ifrag in iterator_insert(
                        len(chainsum),
                        missing_position,
                    )
                ),
                i        = tuple(
                    self.inorm[frag_comb[ifrag].i]
                    if ifrag is not None else None
                    for ichain, ifrag in iterator_insert(
                        len(chainsum),
                        missing_position,
                    )
                ),
                fragtype = tuple(
                    frag_comb[ifrag].fragtype
                    if ifrag is not None else None
                    for ichain, ifrag in iterator_insert(
                        len(chainsum),
                        missing_position,
                    )
                ),
            ) if details else None
        )
    
    def missing_chain(
            self,
            rec,
            missing_position = 1,
            head = None,
            intensity_threshold = 0,
            expected_intensities = None,
            no_intensity_check = False,
            frag_types = None,
            adduct = None,
        ):
        """
        Finds ''missing'' chains i.e. which could complement the chains
        identified among the fragments to fit the species in the record.
        
        Yields tuples with first element a tuple of identified chains and
        as second element the missing chain.
        
        Works a similar way to `chain_combinations`.
        
        Args
        ----
        :param int missing_position:
            Position of the missing chain. 0, 1, 2 are sn1, sn2 and sn3
            positions on glycerol, 0 and 1 are sphingosine base and
            N-acyl in sphingolipids, respectively.
            By default is 1 (sn2 or N-acyl).
        """
        
        chain_list = self.adduct_chain_list(adduct = adduct)
        
        chainsum = rec.chainsum or lipproc.sum_chains(rec.chains)
        
        frags_for_position = self.frags_for_positions(
            rec,
            head = head,
            intensity_threshold = intensity_threshold,
            frag_types = frag_types,
            adduct = adduct,
        )
        
        if missing_position >= len(rec.chainsum.typ):
            
            raise ValueError(
                'No chain known at position %u' % missing_position
            )
        
        if missing_position in frags_for_position:
            
            chains_at_missing = frags_for_position[missing_position]
            del frags_for_position[missing_position]
            
        else:
            
            chains_at_missing = []
        
        # iterate all combinations
        for frag_comb in itertools.product(
            *(
                # making a sorted list of lists from the dict
                i[1] for i in
                sorted(iteritems(frags_for_position), key = lambda i: i[0])
            )
        ):
            
            # if more than one chain missing
            if len(rec.chainsum) - len(frag_comb) > 1:
                
                continue
            
            missing_c = chainsum.c - sum(frag.c for frag in frag_comb)
            missing_u = chainsum.u - sum(frag.u for frag in frag_comb)
            
            # do not yield impossible values
            if missing_c < 1 or missing_u < 0 or missing_u > missing_c - 1:
                
                continue
            
            if (
                # bypass intensity check
                no_intensity_check or
                self._intensity_check(
                    frag_comb, chainsum, expected_intensities
                )
            ):
                
                missing_chain = lipproc.Chain(
                    c = missing_c,
                    u = missing_u,
                    typ = chainsum.typ[missing_position],
                    attr = chainsum.attr[missing_position]
                )
                
                # now all conditions satisfied:
                yield self._chains_frag_comb(
                    frag_comb,
                    chainsum,
                    missing_position = missing_position,
                    missing_chain = missing_chain,
                )
    
    def cu_complete(self, chainsum, chain = None, c = None, u = None):
        """
        Returns the carbon count and unsaturation needed to complete
        the `chain` or `c` and `u` to fit the `chainsum`.
        
        Returns tuple of c and u.
        """
        
        c = c or chain.c
        u = u or chain.u
        
        return chainsum.c - c, chainsum.u - u
    
    def iterrecords(self, adducts = None):
        """
        Iterates MS1 records.
        Yields tuple of adduct type and record.
        """
        
        for add, recs in iteritems(self.ms1_records):
            
            if adducts is None or add in adducts:
                
                for exmass, rec, err in zip(*recs):
                    
                    # if needed we create here the precursor details
                    # it can be included in the identity object
                    if self.add_precursor_details:
                        
                        precursor_details = PrecursorDetails(
                            db = rec.lab.db,
                            db_id = rec.lab.db_id,
                            adduct = add,
                            error = err,
                            exmass = exmass,
                        )
                        
                    else:
                        
                        precursor_details = None
                    
                    yield add, rec, precursor_details
    
    def records_by_type(self, headgroup, sub = (), adducts = None):
        """
        Iterates MS1 database records with a certain headgroup and subtype.
        """
        
        sub = (
            sub
                if type (sub) is set else
            set(sub)
                if type(sub) is list or type(sub) is tuple else
            set([sub])
        )
        
        for add, rec, prec_details in self.iterrecords(adducts = adducts):
            
            if rec.hg and rec.hg.main == headgroup and set(rec.hg.sub) == sub:
                
                yield rec
    
    def first_record(self, headgroup, sub = (), adducts = None):
        """
        Returns the first MS1 database record matching headgroup and subtype.
        """
        
        recbytyp = self.records_by_type(
            headgroup, sub = sub, adducts = adducts
        )
        
        try:
            
            return next(recbytyp)
            
        except StopIteration:
            
            return None
    
    def identify(self, adducts = None):
        
        result = {}
        
        for add, rec, precursor_details in self.iterrecords(adducts):
            
            if rec.hg is None:
                
                continue
            
            rec_str = rec.summary_str()
            
            if rec_str not in result and rec.hg in idmethods[self.ionmode]:
                
                method = idmethods[self.ionmode][rec.hg]
                
                adduct = None if add in {'[M+H]+', '[M-H]-'} else add
                
                result[rec_str] = tuple(
                    method(
                        record = rec,
                        scan = self,
                        adduct = adduct,
                        adduct_str = add,
                        precursor_details = precursor_details,
                    ).identify()
                )
        
        return result
    
    #
    # Sphingolipids
    #
    
    def cer1p_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Ceramide-1-phosphate.

        **Specimen:**
        
        - GLTPD1 - 616.47
      
        **Principle:**
        
        - The most abundant fragment is 78.9591 metaphosphate.
      - If 96.9696 phosphate present adds to the score.
        
        """
        
        score = 0
        fattya = set([])
        if self.most_abundant_mz_is(78.95905658):
            score += 5
            if self.has_mz(96.96962158):
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def hexcer_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Hexosyl-Ceramide.

        **Specimen:**
        
        - GLTP - 744.5627
      
        **Principle:**
        
        - Hexose fragments 71.0115, 89.0220 and 101.0219 must present.
      
        """
        
        score = 0
        fattya = set([])
        
        if all(map(lambda mz: self.mz_among_most_abundant(mz, n = 10),
                   # these are 3 fragments found at GLTP
                   [71.0115000, 89.0220000, 101.021900])):
            
            score += 5
        
        return {'score': score, 'fattya': fattya}
    
    def hexceroh_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Hexosyl-Ceramide-OH
        ('t'). This method is the same as `hexcer_neg_1`.

        **Specimen:**
        
        - GLTP - 760.557
      
        **Principle:**
        
        - Hexose fragments 71.0115, 89.0220 and 101.0219 must present.
      
        """
        
        return self.hexcer_neg_1()
    
    def sm_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Sphingomyeline.

        **Specimen:**
        
        - GLTPD1 - 745.55
      
        **Principle:**
        
        - Must have a neutral loss of CH3+COOH (60.0211).
      - Phosphate+choline-CH3 fragment 168.0431 must be present.
        
        """
        
        score = 0
        fattya = set([])
        
        if self.mz_among_most_abundant(168.0431206) and self.has_nl(60.02113):
            score += 5
        
        return {'score': score, 'fattya': fattya}
    
    def sph1p_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Spingosine-1-phosphate.

        **Specimen:**
        
        - Only observed in standard.
      
        **Principle:**
        
        - Phosphate 78.9590 must be present.
      
        """
        
        score = 0
        fattya = set([])
        
        if self.has_mz(78.95905658):
            score += 5
        
        return {'score': score, 'fattya': fattya}
    
    def cer_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is Ceramide.

        **Specimen:**
        
        - SEC14L1 - 582.509
      
        **Principle:**
        
        - A Ceramide backbone fragment must be among the 2 most abundant.
      - Ceramide backbone fragments lighter by N or C2N but same carbon
          count and unsaturation add to the score.
        
        """
        
        score = 0
        fattya = set([])
        
        if self.fa_among_most_abundant('CerFA', n = 2):
            
            score += 5
            fattya = self.fa_combinations('Cer', sphingo = True)
            fa_h_ccs = self.matching_fa_frags_of_type('Cer', 'CerFA(')
            
            for fa_h_cc in fa_h_ccs:
                
                for fa_other in [
                    '[CerFA-N(C%u:%u)-]-',
                    '[CerFA-C2N(C%u:%u)-]-']:
                    
                    if self.frag_name_present(fa_other % fa_h_cc):
                        
                        score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def cerp_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Ceramide-1-phosphate.
        Gives similar result as Sphingosine-1-phosphate.

        **Specimen:**
        
        - GLTPD1 - 616.47
      
        **Principle:**
        
        - The most abundant fragment must be 78.9591 metaphosphate.
      - Presence of 96.9696 phosphate increase the score.
        
        """
        
        score = 0
        fattya = set([])
        
        if self.most_abundant_mz_is(78.95905658):
            
            score += 5
            
            if self.has_mz(96.96962158):
                
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def cerp_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Ceramide-1-phosphate.

        **Specimen:**
        
        - GLTPD1 + 728.59, 590.45, 702.58, 618.430, 616.415, 640.409
      
        **Principle:**
        
        - A sphingosine fragment with double H2O loss must be among the three
        highest abundant fragments.
        
        """
        
        score = 0
        fattya = set([])
        
        if self.fa_among_most_abundant('-H2O-H2O+]+', n = 3, sphingo = True):
            score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def adduct(self, adduct, reset = False):
        """
        Creates a copy of the current scan assuming that the precursor
        is a certain adduct. The precursor will be converted to [M+H]+
        or [M-H]- adduct and neutral losses will be calculated accordingly.
        """
        
        if adduct in self.adducts and not reset:
            
            return
        
        ad2ex = settings.get('ad2ex')[1][self.ionmode][adduct]
        ex2ad = 'remove_h' if self.ionmode == 'neg' else 'add_h'
        
        fake_precursor = (
            getattr(
                mzmod.Mz(
                    getattr(
                        mzmod.Mz(self.precursor),
                        ad2ex
                    )()
                ),
                ex2ad
            )()
        )
        
        annot = self.get_annot(fake_precursor)
        
        chain_list = self._build_chain_list(annot = annot)
        
        self.adducts[adduct] = {
            'fake_precursor': fake_precursor,
            'annot': annot,
            'chain_list': chain_list,
        }
    
    def adduct_annot(self, adduct = None):
        """
        Gets the annotations for a certain adduct.
        """
        
        return self.adduct_data('annot', adduct = adduct)
    
    def adduct_chain_list(self, adduct = None):
        """
        Gets the chain list for a certain adduct.
        """
        
        return self.adduct_data('chain_list', adduct = adduct)
    
    def adduct_data(self, name, adduct = None):
        
        if adduct is None:
            
            return getattr(self, name)
        
        self.adduct(adduct)
        
        return self.adducts[adduct][name]
    
    def get_ms1_records(
            self,
            hg,
            subtype = None,
            sph = None,
            ether = None,
            oh = None,
            databases = None,
        ):
        """
        Iterates MS1 records for a given type.
        
        Yields tuples of record, adduct type and error in ppm
        """
        
        subtype = subtype or ()
        
        hg = (
            hg
            if isinstance(hg, lipproc.Headgroup) else
            lipproc.Headgroup(main = hg, sub = subtype)
        )
        
        for add, recs in iteritems(self.ms1_records):
            
            for rec_mz, rec, err_ppm in zip(*recs):
                
                if (
                    rec.hg == hg and (
                        databases is None or rec.lab.db in databases
                    ) and (
                        sph is None or rec.chainsum.attr.sph == sph
                    ) and (
                        ether is None or rec.chainsum.attr.ether == ether
                    ) and (
                        oh is None or rec.chainsum.attr.oh == oh
                    )
                ):
                    
                    yield rec, add, err_ppm


class AbstractMS2Identifier(object):
    
    class_methods = {}
    subclass_methods = {}
    
    def __init__(
            self,
            record,
            scan,
            adduct = None,
            adduct_str = None,
            precursor_details = None,
            missing_chains = None,
            explicit_and_implicit = False,
            must_have_chains = True,
            chain_comb_args = {},
            missing_chain_args = {},
        ):
        
        self.score = 0
        self.max_score = 0
        self.rec = record
        self.scn = scan
        self.add = adduct
        self.adduct_str = adduct_str
        self.precursor_details = precursor_details
        self.missing_chains = (
            missing_chains if missing_chains is not None else
            tuple(range(len(record.chainsum))) # any chain can be missing
        )
        self.chain_comb_args = chain_comb_args
        self.missing_chain_args = missing_chain_args or self.chain_comb_args
        self.explicit_and_implicit = explicit_and_implicit
        self.must_have_chains = must_have_chains
        
        self.scores = {}
    
    def identify(self):
        
        if not self.rec.hg:
            
            return
        
        self.confirm_class()
        
        chains_confirmed = False
        
        for chains in self.confirm_chains_explicit():
            
            yield MS2Identity(
                max(self.score, 0),
                self.max_score,
                self.percent_score(),
                self.rec.hg,
                self.rec.chainsum,
                chains = chains[0],
                chain_details = chains[1],
                scan_details = self.scn.scan_details,
                precursor_details = self.precursor_details,
            )
            chains_confirmed = True
        
        if not chains_confirmed or self.explicit_and_implicit:
            
            for chains in self.confirm_chains_implicit():
                
                yield MS2Identity(
                    max(self.score, 0),
                    self.max_score,
                    self.percent_score(),
                    self.rec.hg,
                    self.rec.chainsum,
                    chains = chains[0],
                    chain_details = chains[1],
                    scan_details = self.scn.scan_details,
                    precursor_details = self.precursor_details,
                )
                
                chains_confirmed = True
        
        if not chains_confirmed and not self.must_have_chains and self.score:
            
            yield MS2Identity(
                max(self.score, 0),
                self.max_score,
                self.percent_score(),
                self.rec.hg,
                self.rec.chainsum,
                chains = None,
                chain_details = None,
                scan_details = self.scn.scan_details,
            )
    
    def percent_score(self):
        """
        Returns the score as a percentage of the maximum possible score.
        
        Zero maximum score means something is wrong, then it returns 200.
        """
        
        return (
            max(int(np.round(self.score / self.max_score * 100.)), 0)
            if self.max_score else
            200
        )
    
    def confirm_class(self):
        """
        In this base class pass through everything.
        Most of the subclasses should override this.
        """
        
        self.score = 0
        
        if self.rec.hg is not None and self.rec.hg.main in self.class_methods:
            
            score, max_score = getattr(
                self,
                self.class_methods[self.rec.hg.main]
            )()
            
            self.score += score
            self.max_score += max_score
    
    def confirm_subclass(self):
        
        subclasses = self.rec.hg.sub or ('empty',)
        
        if self.rec.hg is not None:
            
            for sub in subclasses:
                
                if sub not in self.scores and sub in self.subclass_methods:
                    
                    score, max_score = getattr(
                        self,
                        self.subclass_methods[sub]
                    )()
                    
                    self.scores[sub] = score
                    self.score += score
                    self.max_score += max_score
    
    def confirm_chains_explicit(self):
        
        return self.scn.chain_combinations(self.rec, **self.chain_comb_args)
    
    def confirm_chains_implicit(self):
        
        for missing in self.missing_chains:
            
            for chain_comb in self.scn.missing_chain(
                self.rec,
                missing_position = missing,
                **self.missing_chain_args
            ):
                
                yield chain_comb
    
    def matching_chain_combinations(
            self,
            chain_param1,
            chain_param2,
            score_method = lambda ccomb: (min(ccomb, 3) * 2, 6),
        ):
        
        ccomb = len(list(
            self.scn.matching_chain_combinations(
                self.rec,
                chain_param = (chain_param1, chain_param2),
            )
        ))
        
        score, max_score = score_method(ccomb)
        
        self.score += score
        self.max_score += score
    
    def check_lyso(self, score_threshold = 5):
        """
        Checks whether the this mass has been identified in the database
        as a lyso species and calls the corresponding lyso identification
        method.
        
        Returns ``True`` if the score from the lyso is larger than
        ``score_threshold``.
        """
        
        rec_lyso = self.scn.first_record(self.rec.hg.main, sub = ('Lyso',))
        
        if rec_lyso:
            
            lyso_hg = lipproc.Headgroup(
                main = self.rec.hg.main,
                sub = ('Lyso',),
            )
            lyso = idmethods[self.scn.ionmode][lyso_hg](rec_lyso, self.scn)
            lyso.confirm_class()
            
            return lyso.score > score_threshold
        
        return False


#
# Lipid identification methods
#

#
# Fatty acids
#

class FA_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is a fatty acid.
    Here we only check if the most abundant fragment is the
    fatty acid itself.

    **Specimen:**
    
    - in vitro FABP1 -
    
    **Principle:**
    
    - The most abundant fragment must be a fatty acid which matches
      the carbon count and the unsaturation of the whole molecule.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {
                'head': 1,
                'frag_types': {
                    0: {'FA-H'}
                }
            },
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 10
        
        if (
            self.rec.chainsum and
            self.scn.chain_among_most_abundant(
                frag_type = 'FA-H',
                c = self.rec.chainsum.c,
                u = self.rec.chainsum.u,
            )
        ):
            
            self.score = 10


class FA_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a fatty acid.
    Here we only check if the most abundant fragment is the
    fatty acid itself.

    **Specimen:**
    
    - Not known
    
    **Principle:**
    
    - The most abundant fragment must be a fatty acid which matches
      the carbon count and the unsaturation of the whole molecule.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {
                'head': 1
            }
        )
    
    def confirm_class(self):
        
        self.max_score = 10
        
        if (
            self.rec.chainsum and
            self.scn.chain_among_most_abundant(
                frag_type = 'FA+H',
                c = self.rec.chainsum.c,
                u = self.rec.chainsum.u,
            )
        ):
            
            self.score = 10

#
# Glycerolipids
#

class DAG_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a DAG.

    **Specimen:**
    
    - in vivo: SEC14L2 + 584.52
    - in vitro: BNIP2 + 770.67
    
    **Principle:**
    
    - Combination of fatty acid fragments among the 10 most abundant
      fragments must match the expected carbon count and unsaturation.
    - If these are among the 5 highest fragments the score is higher.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 6
        
        if self.scn.has_chain_combinations(self.rec, head = 10):
            
            self.score += 4
        
        if self.scn.has_chain_combinations(self.rec, head = 6):
            
            self.score += 2


class DAG_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is a DAG.

    **Specimen:**
    
    - We don't have yet.
    
    **Principle:**
    
    - Combination of fatty acid fragments among the 10 most abundant
      fragments must match the expected carbon count and unsaturation.
    - If these are among the 5 highest fragments the score is higher.
    
    (Same as in positive ionmode.)
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 6
        
        if self.scn.has_chain_combinations(self.rec, head = 10):
            
            self.score += 4
        
        if self.scn.has_chain_combinations(self.rec, head = 6):
            
            self.score += 2


class TAG_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a TAG.

    **Specimen:**
    
    - STARD11 + 818.7187
    
    **Principle:**
    
    - Combination of fatty acid fragments must match the expected
      carbon count and unsaturation.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 10
        
        if self.scn.has_chain_combinations(self.rec, head = 15):
            
            self.score += 5
        
        if self.scn.has_chain_combinations(self.rec, head = 7):
            
            self.score += 5


class TAG_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is a TAG.

    **Specimen:**
    
    - We don't have yet.
    
    **Principle:**
    
    - Combination of fatty acid fragments must match the
      expected carbon count and unsaturation.
    
    (Same as in positive ionmode.)
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.score = 0
        self.max_score = 5
        
        if self.scn.has_chain_combinations(self.rec):
            
            self.score += 5


class GL_Positive(AbstractMS2Identifier):
    """
    Generic class for identification of glycerolipids.
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    class_methods = {
        'DGTS': 'dgts',
        'DGCC': 'dgcc',
        'DGTA': 'dgts',
        'DGDG': 'dgdg',
        'MGDG': 'mgdg',
        'SQDG': 'sqdg',
    }
    
    def confirm_class(self):
        
        if self.rec.hg.main in self.class_methods:
            
            getattr(self, self.class_methods[self.rec.hg.main])()
    
    def dgts(self):
        
        self.max_score += 30
        
        if self.scn.has_fragment('DGTS [G+TS] (236.1492)'):
            
            self.score += 10
            
            self.score += sum(map(bool,
                (
                    self.scn.has_fragment('DGTS [TS] (144.1019)'),
                    self.scn.has_chain_fragment_type('NL FA-H2O'),
                )
            )) * 10
    
    def dgcc(self):
        
        self.max_score += 20
        
        if self.scn.has_fragment('PC/SM [Ch+H2O] (104.107)'):
            
            self.score += 10
            
            if self.scn.has_fragment('DGCC [C2+Ch] (132.1388)'):
                
                self.score += 10
    
    def sqdg(self):
        
        self.max_score += 10
        
        if self.scn.has_fragment('NL [Hexose+SO3+H2O+H] (NL 261.0280)'):
            
            self.score += 10
    
    def mgdg(self):
        
        self.max_score += 10
        
        if self.scn.has_fragment('[Hexose+H2O-H] (NL 197.07)'):
            
            self.score += 10
    
    def dgdg(self):
        
        self.max_score += 10
        
        if self.scn.has_fragment('NL [2xHexose+H2O-H] (NL 359.1190)'):
            
            self.score += 10


class GL_Negative(AbstractMS2Identifier):
    """
    Generic class for identification of glycerolipids.
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score += 10
        
        if self.scn.has_chain_combination(
            record = self.rec,
            chain_param = (
                {'frag_type': {'FA-H', 'FA-'}},
            )
        ):
            
            self.score += 10

#
# Glycerophospholipids
#

class PE_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is Phosphatidylethanolamine.

    **Specimen:**
    
    - GM2A - 714.507 and 716.523
    
    **Principle:**
    
    - The most abundant fragment is a fatty acid [M-H]- ion.
    - 140.0118 PE headgroup must be present.
    - Other headgroup ions 196.0380 and 178.0275 add to the score.
    - Lyso-PE and [M-H-CO2]- fatty acid fragments complementing the
      highest [M-H]- fatty acid increase the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score += 11
        
        if (
            self.scn.chain_fragment_type_is(
                i = 0,
                chain_type = 'FA',
                frag_type  = 'FA-H'
            ) and
            self.scn.has_fragment('PE [P+E] (140.0118)')
        ):
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('PE [G+P+E-H2O] (196.0380)'),
                self.scn.has_fragment('PE [G+P+E] (178.0275)'),
            ))) * 3
            
            # by default this returns max 6
            self.matching_chain_combinations(
                {'frag_type': 'FA-H'},
                {'frag_type': {
                        'LysoPE',
                        'LysoPEAlkyl',
                        'LysoPEAlkyl-H2O',
                        'FA-H2O-H'
                    }
                }
            )


class PE_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a
    Phosphatidylethanolamine.

    **Specimen:**
    
    - in vivo BPI + 718.536
    - Lyso-PE: in vitro FABP1 + 454.29
    
    **Principle:**
    
    - The PE headgroup neutral loss 141.0191 has the highest intensity.
    - If it is a Lyso-PE score will be zero.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 20
        
        if self.scn.has_fragment('NL PE [P+E] (NL 141.0191)'):
            
            if self.check_lyso():
                
                return
            
            if self.scn.has_fragment('PE [P+E] (142.0264)'):
                
                self.score += 5
            
            self.score += 15


class LysoPE_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a
    Lysophosphatidylethanolamine.

    **Specimen:**
    
    - in vitro FABP1 + 454.29
    
    **Principle:**
    
    - The PE headgroup neutral loss 141.0191 has the highest intensity.
    - A fatty acid-glycerol fragment should match the carbon count and
      unsaturation of the whole molecule.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 15
        
        if self.scn.has_fragment('NL PE [P+E] (NL 141.0191)'):
            
            self.score = 5
            
            if self.scn.has_fragment('PE [P+E] (142.0264)'):
                
                self.score += 5
            
            self.scn.build_chain_list()
            
            if (
                len(self.scn.chain_list) and
                self.scn.chain_fragment_type_is(
                    self.scn.chain_list[0].i,
                    frag_type = 'FA+Glycerol-OH',
                    c = self.rec.chainsum.c,
                    u = self.rec.chainsum.u,
                )
            ):
                
                self.score += 5


class PC_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is a Phosphatidylcholine.

    **Specimen:**
    
    - BPI - 804.57 and 776.545
    
    **Principle:**
    
    - 168.0431 phosphate+choline-CH3 fragment must be present.
    - The highest abundant fragment must be a fatty acid [M-H]- fragment.
    - Lyso-PC fragments complementing the highest [M-H]- fatty acid
      increase the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 17
        
        if (
            self.scn.chain_fragment_type_is(
                i = 0,
                chain_type = 'FA',
                frag_type = 'FA-H'
            ) and
            self.scn.has_fragment('PC/SM PO4+choline-CH3 (168.0431)')
        ):
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('PE [G+P+E-H2O] (196.0380)'),
                self.scn.has_fragment('PE [G+P+E] (178.0275)'),
            ))) * 3
            
            self.matching_chain_combinations(
                {'frag_type': 'FA-H'},
                {'frag_type': 'LysoPC'},
                score_method = lambda ccomb: ((ccomb > 1) * 6, 6),
            )


class PC_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a Phosphatidylcholine.

    **Specimen:**
    
    - BPI + 786.607
    
    **Principle:**
    
    - The most abundant fragment must be choline+phosphate 184.0733.
    - The 86.0964 ethyl-trimetylammonium must be present.
    - The most abundant fatty acid can not have the same carbon count
      and unsaturation as the whole molecule (then it is Lyso-PC).
    - Fragments 104.1069, 124.9998, 60.0808 and 58.0651 increase the
      score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = False,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 13
        
        if (
            self.scn.fragment_percent_of_most_abundant(
                'PC/SM [P+Ch] (184.0733)', 10.0
            ) and
            self.scn.has_fragment('PC/SM [Ch] (86.096)')
        ):
            
            if self.check_lyso():
                
                return
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('PC/SM [Ch+H2O] (104.107)'),
                self.scn.has_fragment('PC/SM [P+Et] (124.9998)'),
                self.scn.has_fragment('PC/SM [N+3xCH3] (60.0808)'),
                self.scn.has_fragment('PC/SM [Ch-Et] (58.0651)'),
            ))) * 2


class LysoPC_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a Lysophosphatidylcholine.
    
    **Specimen:**
    
    - in vitro FABP1 + 522.36
    
    **Principle:**
    
    - Choline-phosphate 184.0733, ethyl-trimethylammonium 86.0964 and
      neutral loss 183.0660 must be present.
    - The latter neutral loss corresponds to a fatty acid+glycerol ion.
    - The carbon count and unsaturation of this fragment should match
      that of the whole molecule.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = False,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 15
        
        if (
            self.scn.most_abundant_fragment_is('PC/SM [P+Ch] (184.0733)') and
            self.scn.has_fragment('PC/SM [Ch] (86.096)')
        ):
            
            self.score += 5
            
            if self.scn.has_fragment('NL PC/SM [P+Ch] (NL 183.066)'):
                
                self.score += 5
            
            if self.scn.has_chain_fragment_type(
                frag_type = {'FA+Glycerol-OH', 'NL FA-H2O'},
                c = self.rec.chainsum.c,
                u = self.rec.chainsum.u,
            ):
                
                self.score += 5


class PI_Negative(AbstractMS2Identifier):
    """
    Examines if a negative MS2 spectrum is Phosphatidylinositol.

    **Specimen:**
    
    - GM2A - 835.52
    
    **Principle:**
    
    - Inositolphosphate-H2O fragment 241.0119, metaphosphate 78.9591 and
      headgroup fragment 152.9958 must be present.
    - Additional headgroup fragments 96.9696, 259.0224 and 297.0381
      increase the score.
    - Presence of Lyso-PI fragments complementing other [M-H]- fatty
      acid fragments increase the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 19
        
        if (
            self.scn.has_fragment('PI [InsP-H2O]- (241.01)') and
            self.scn.has_fragment('PA/PG/PI/PS [G+P] (152.9958)') and
            self.scn.has_fragment('Cer1P/PIP/PL metaphosphate (78.9591)')
        ):
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('Cer1P/PI phosphate (96.9696)'),
                self.scn.has_fragment('PI [InsP-H]- (259.02)'),
                self.scn.has_fragment('PI [G+P+I] (297.04)'),
                self.scn.has_fragment('PI [InsP-2H2O]- (223.00)'),
            ))) * 2
            
            self.matching_chain_combinations(
                {'frag_type': 'FA-H'},
                {'frag_type': {
                        'LysoPI',
                        'LysoPI-H2O',
                    }
                },
                score_method = lambda ccomb: (min(ccomb, 2) * 3, 6),
            )


class PI_Positive(AbstractMS2Identifier):
    """
    Examines if a negative MS2 spectrum is Phosphatidylinositol.

    **Specimen:**
    
    - SEC14L2 + 906.60 and 882.6
    
    **Principle:**
    
    - Combinations of fatty acid fragments must match the expected
      carbon count and unsaturation for PI.
    - Presence of neutral losses 259.0219 and 277.0563 adds to the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 9
        
        if self.scn.has_chain_combinations(self.rec):
            
            self.score += 1
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('NL PI [P+Ins] (NL 259.0219)'),
                self.scn.has_fragment('NL PI [P+Ins+NH3] (NL 277.0563)'),
            ))) * 4


class PS_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is a Phosphatidylserine.

    **Specimen:**
    
    - ORP9 - 788.54
    
    **Principle:**
    
    - The most abundant fragment is an [M-H]- fatty acid fragment.
    - Glycerophosphate fragment 152.9958 must be present.
    - Metaphosphate 78.9591 increases the score.
    - Serine-H2O neutral loss 87.0320 adds to the score.
    - Presence of Lyso-PS and Lyso-PA fragments complementing
      the highest [M-H]- fatty acid fragment increase the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 17
        
        if (
            self.scn.has_chain_combinations(self.rec) and
            self.scn.chain_fragment_type_is(
                0, chain_type = 'FA', frag_type = 'FA-H'
            ) and
            self.scn.fragment_among_most_abundant(
                'PA/PG/PI/PS [G+P] (152.9958)', 5
            )
        ):
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('Cer1P/PIP/PL metaphosphate (78.9591)'),
                self.scn.has_fragment('PS [Ser-H2O] (87.0320)'),
            ))) * 3
            
            self.matching_chain_combinations(
                {'frag_type': 'FA-H'},
                {'frag_type': {'LysoPS', 'LysoPA'}},
                score_method = lambda ccomb: (min(ccomb, 2) * 3, 6),
            )


class PS_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a Phosphatidylserine.

    **Specimen:**
    
    - BPI + 790.56
    
    **Principle:**
    
    - PS headgroup neutral loss 185.0089 must be the highest intensity.
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 5
        
        if self.scn.fragment_among_most_abundant('PS [P+S] (NL 185.0089)', 1):
            
            self.score += 5


class PG_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is Phosphatidylglycerol.
    The result will be the same as `bmp_neg_1`, as in negative
    mode we do not know a way to distinguish these species.
    

    **Specimen:**
    
    - GM2A - 799.54
    - BPIFB2 - 773.5258 (might be BMP)
    
    **Principle:**
    
    - The most abundant fragment is a fatty acid [M-H]- ion.
    - The 152.9958 glycerophosphate fragment must be present.
    - If Lyso-PG fragment present with carbon count complementing
      the [M-H]- fatty acid score is higher.
    - Presence of 171.0064 headgroup fragment adds to the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 14
        
        if (
            self.scn.has_chain_combinations(self.rec) and
            self.scn.chain_fragment_type_is(
                0, chain_type = 'FA', frag_type = 'FA-H'
            ) and
            self.scn.has_fragment('PA/PG/PI/PS [G+P] (152.9958)')
        ):
            
            self.score += 5
            
            if self.scn.has_fragment('PG headgroup (171.0064)'):
                
                self.score += 3
            
            self.matching_chain_combinations(
                {'frag_type': 'FA-H'},
                {'frag_type': {
                        'LysoPG',
                        'LysoPG-H2O',
                    }
                },
                score_method = lambda ccomb: (min(ccomb, 2) * 3, 6),
            )


class PG_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a Phosphatidylglycerol.
    At in vivo observed only in standard.
    
    **Principle:**
    
    - The PG headgroup neutral loss (189.0402) is the fragment ion
      with the highest intensity?
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 5
        
        if (
            self.scn.most_abundant_fragment_is(
                'NL PG [G+P+NH3] (NL 189.0402)'
            )
        ):
            
            self.score += 5
        
        # alternative for lyso
        if (
            self.rec.hg.sub == ('Lyso',) and
            self.scn.chain_fragment_type_among_most_abundant(
                n = 1,
                frag_type = {'FA+Glycerol-OH'},
                c = self.rec.chainsum.c,
                u = self.rec.chainsum.u,
            )
        ):
            
            self.score += 5
            self.max_score += 5


class BMP_Negative(PG_Negative):
    """
    Examines if a negative mode MS2 spectrum is Bismonoacylglycerophosphate.
    The result will be the same as for PG, as in negative
    mode we do not know a way to distinguish these species.
    

    **Specimen:**
    
    - GM2A - 799.54
    - BPIFB2 - 773.5258 (might be BMP)
    
    **Principle:**
    
    - The most abundant fragment is a fatty acid [M-H]- ion.
    - The 152.9958 glycerophosphate fragment must be present.
    - If Lyso-PG fragment present with carbon count complementing
      the [M-H]- fatty acid score is higher.
    - Presence of 171.0064 headgroup fragment adds to the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        PG_Negative.__init__(self, record, scan, **kwargs)


class BMP_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum
    is a Bismonoacylglycerophosphate.

    **Specimen:**
    
    - BPIFB2 + 792.57
    
    **Principle:**
    
    - A glycerol+fatty acid fragment can be found among the 3 highest?
    - The PG headgroup neutral loss (189.0402) is among the fragments?
    - If so, does it have a lower intensity than half of the fatty
      acid+glycerol fragment?
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 10
        
        if (
            self.scn.has_chain_combinations(self.rec, head = 15) and
            self.scn.chain_fragment_type_among_most_abundant(
                chain_type = 'FA', frag_type = 'FA+Glycerol-OH', n = 3
            )
        ):
            
            self.score += 5
            
            i_hg = self.scn.fragment_by_name('NL PG [G+P+NH3] (NL 189.0402)')
            
            if i_hg is not None:
                
                i_gfa = self.scn.highest_fragment_by_chain_type(
                    head = 4, frag_type = 'FA+Glycerol-OH'
                )
                
                if self.scn.intensities[i_gfa] < self.scn.intensities[i_hg]:
                    
                    self.score = 0
                
                else:
                    
                    self.score += 5


class PA_Negative(AbstractMS2Identifier):
    """
    Examines if a negative mode MS2 spectrum is Phosphatidylglycerol.
    The result will be the same as `bmp_neg_1`, as in negative
    mode we do not know a way to distinguish these species.
    

    **Specimen:**
    
    - GM2A - 799.54
    - BPIFB2 - 773.5258 (might be BMP)
    
    **Principle:**
    
    - The most abundant fragment is a fatty acid [M-H]- ion.
    - The 152.9958 glycerophosphate fragment must be present.
    - If Lyso-PG fragment present with carbon count complementing
      the [M-H]- fatty acid score is higher.
    - Presence of 171.0064 headgroup fragment adds to the score.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 25
        
        if (
            self.scn.has_chain_combinations(self.rec) and
            self.scn.chain_fragment_type_is(
                0, chain_type = 'FA', frag_type = 'FA-H'
            ) and
            self.scn.fragment_among_most_abundant(
                'PA/PG/PI/PS [G+P] (152.9958)', 10
            ) and
            self.scn.fragment_among_most_abundant(
                'Cer1P/PIP/PL metaphosphate (78.9591)', 10
            )
        ):
            
            self.score += 20
            
            if self.scn.has_fragment('Cer1P/PI phosphate (96.9696)'):
                
                self.score += 5


class PA_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a Phosphatidylglycerol.
    At in vivo observed only in standard.
    
    **Principle:**
    
    - The PG headgroup neutral loss (189.0402) is the fragment ion
      with the highest intensity?
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 20
        
        if (
            self.scn.fragment_among_most_abundant(
                'NL [P] (NL 97.9769)', 3, adduct = self.add
            )
        ):
            
            self.score += 10
            
            if self.scn.has_chain_combination(
                self.rec,
                chain_param = ({
                    'frag_type': {
                        'FA+Glycerol-OH',
                        'FA-OH',
                        'FA-H2O-OH',
                    }
                },)
            ):
                
                self.score += 10

#
# Vitamins
#

class VA_Positive(AbstractMS2Identifier):
    """
    Examines if a positive MS2 spectrum is vitamin A (retinol).

    **Specimen:**
    
    - in vivo RBP1 + 269.2245
    - in vivo RBP4 + 269.2245
    
    **Principle:**
    
    - The most abundant ion is the whole molecule m/z = 269.224.
    - Presence off 3 other ions adds to the score but not
      mandatory: 213.165, 145.1027, 157.1028.
    
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = False,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 8
        
        if self.scn.fragment_among_most_abundant('Retinol I (269.2264)', 3):
            
            self.score += 5
            
            self.score += sum(map(bool, (
                self.scn.has_fragment('Retinol II (213.1637)'),
                self.scn.has_fragment('Retinol III (157.1012)'),
                self.scn.has_fragment('Retinol IV (145.1012)'),
            )))


class VA_Negative(AbstractMS2Identifier):
    """
    Examines if a positive MS2 spectrum is vitamin A (retinol).

    **Specimen:**
    
    - Standards 141020 negative scan 673
    
    **Principle:**
    
    - 3 fragments seems to always present and be among the most abundant:
      79.055, 119.087 and 255.212; presence of these is the main condition.
    - We also detected 125.061 in our standards which is special because
      contains 2 oxygens; presence of this increase the score.
    """
    
    def __init__(self, record, scan, **kwargs):
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = False,
            **kwargs,
        )
    
    def confirm_class(self):
        
        self.max_score = 8
        
        if all((
            self.scn.fragment_among_most_abundant(fragname, 7)
            for fragname in (
                'Retinoic acid I (79.0553)',
                'Retinoic acid II (119.0866)',
                'Retinoic acid IV (255.2118)',
            )
        )):
            
            self.score += 5
            
            if self.scn.has_fragment('Retinoic acid III (125.0608)'):
                
                self.score += 3

#
# Sphingolipids
#

class Cer_Positive(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a ceramide.
    Identifies ceramide varieties including sphingomyeline,
    ceramide-1-phosphare, ceramide-phosphoethanolamine,
    OH-acyl-ceramide, hexosyl and dihexosyl-ceramides,
    and d, t and DH long chain base varieties.
    
    dCer
    ====
    
    **Specimen:**
    
    - Sphingolipid standards scans 
    - SEC14L1 + 538.52
    - STARD11 + 538.526
    
    **Principle:**
    
    - A sphingosine backbone with two H2O loss must be among the
      10 most abundant fragments.
    - Fatty acid [M+H]+ or [M-O]+ fragments or neutral losses
      complementing the one above increase the score.
    - Sphingosine backbone fragments with same carbon count and
      unsaturation with the one with 2 water loss but [Sph-C-2(H2O)]+
      or [Sph-H2O]+ add to the score.
    - The score increases if the following choline fragments
      can not be found: 58.0651, 104.1070, 124.9998 and 184.0733.
    - The presence of the following fragments increase the score:
      60.0444, 70.0651, 82.0651, 96.0808, 107.0730, 121.0886,
      135.1042 and 149.1199.
    
    DHCer
    =====
    
    **Specimen:**
    
    - Standards 180626 m/z 568.56 scan 2367
    
    **Principle:**
    
    - Same pattern as at dCer but from the sphingosine derived fragments
      it becomes clear if it has no unsaturation.
    
    tCer
    ====
    
    **Specimen:**
    
    - Standards 180628 m/z 584.56 scan 2070
    
    **Principle:**
    
    - Has Sph+H2O-H fragment which does not occur at dCer.
    - Sph-H2O-H and Sph-H are much higher abundant than at dCer.
    - Strong fragment at 60.044 which is missing at dCer.
    - H2O and 2xH2O neutral losses are much higher than at dCer.
    
    dCer-2OH-acyl
    =============
    
    **Specimen:**
    
    - Standards 180615 m/z 584.22 scan 2421
    
    **Principle:**
    
    - Same as other ceramides. It's d and DH forms are isobaric with tCer
      but d and t are clearly distinguishable so this does not cause
      confusion.
    
    dCer-1-P
    ========
    
    **Specimen:**
    
    - in vivo GLTPD1 + 728.59
    
    **Principle:**
    
    - A shpingosine backbone with 2 H2O loss must be among the 3 highest
      intensity fragments.
    - Presence of any of the following fragments increases the score:
      82.0651, 115.9875.
      107.0729, 135.1043, 149.1199.
    
    Hex-dCer
    ========
    
    **Specimen:**
    
    - in vivo GLTP + 810.68
    
    **Principle:**
    
    - Hexose fragments 198.0740, 180.0634 and 162.0528 must present.
      These are neutral losses of hexose, hexose-H2O and hexose+H2O
    
    Hex-tCer
    ========
    
    **Specimen:**
    
    - in vivo GLTP + 826.67
    - in vitro GLTP + 826.67, 800.66,
    
    **Principle:**
    
    - Hexose fragments 198.0740, 180.0634 and 162.0528 must present.
      These are neutral losses of hexose, hexose-H2O and hexose+H2O
    
    Hex2-dCer
    =========
    
    **Specimen:**
    
    - in vivo GLTP + 988.73
    
    **Principle:**
    
    - Loss of double hexose with our without extra water or water loss
      are the characteristic fragments of this class.
    
    SHex-dCer
    =========
    
    **Specimen:**
    
    - in vitro 890.64
    
    dSM & DHSM
    ==========
    
    **Specimen:**
    
    - in vivo GLTPD1 + 703.57
    - in vitro GLTPD1 + 813.68
    
    **Principle:**
    
    - The following choline fragments must be present: 60.0808, 86.0964,
      104.1069, 124.9998 and 184.0733. The last one is the most intensive.
    - If 58.0651 can be found it adds to the score.
    - dSM and DHSM are not distinguishable in our settings. Maybe the
      [Sph-2xH2O+H]+ ion (264 @ 18:1) presents more often at d and only
      eventually at DH.
    
    PE-Cer
    ======
    
    We do not have this in standards or in screens so we can not test this.
    Based on Amiar 2016 and Narayanaswamy 2014.
    
    **Principle:**
    
    - Neutral loss of 141.0191 must be present.
    - 142.0264 phospho-ethanolamine fragment and neutral loss of
      phospho-ethanolamine+water might be present.
    - Sph-2xH2O fragment increases the score.
    
    """
    
    class_methods = {
        'SM': 'sm',
        'Sph': 'sph',
    }
    
    subclass_methods = {
        '1P': 'cer1p',
        'Hex': 'hexcer',
        'Hex2': 'hex2cer',
        'SHex': 'shexcer',
        'SHex2': 'shex2cer',
        'PE': 'pe_cer',
        'M2': 'm2',
        'M1': 'm1',
        'M3': 'm3',
        'PC': 'pc',
        'empty': 'cer',
    }
    
    def __init__(self, record, scan, **kwargs):
        
        self.nacyl = record.chainsum is not None and len(record.chainsum) > 1
        self.oacyl = record.chainsum is not None and len(record.chainsum) > 2
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (1,) if self.nacyl else (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
        
        self.sph_scores     = {}
        self.sph_max_scores = {}
        self.fa_scores      = {}
        self.fa_max_scores  = {}
    
    def confirm_class(self):
        
        AbstractMS2Identifier.confirm_class(self)
        
        self.max_score += 2
        
        self.score += sum(map(
            self.scn.has_fragment,
            (
                'NL [H2O] (NL 18.0106)',
                'NL [2xH2O] (NL 36.0211)',
            )
        ))
        
        self.confirm_subclass()
    
    def confirm_chains_explicit(self):
        """
        Most of the time we don't really have fatty acid derived
        fragments from ceramides in positive mode. However certain
        sphingosin base derived fragments correspond to neutral losses
        of fatty acids as the molecule is composed of a sphingosine
        and a fatty acid but this is redundant. Hence we call both
        explicit and implicit identification as practically the it
        is implicit anyways.
        """
        
        for chains in itertools.chain(
            AbstractMS2Identifier.confirm_chains_explicit(self),
            AbstractMS2Identifier.confirm_chains_implicit(self),
        ):
            
            if chains[0][0].attr.sph == self.rec.chainsum.attr[0].sph:
                
                # the sphingosin base and fatty acyl related part of the
                # score is valid only for the current chain combination
                # hence now we add these to the overall score, yield the
                # identification and then subtract them from the score
                sph_score, sph_max_score = self.sphingosine_base(
                    chains[0][0].attr.sph
                )
                
                self.score += sph_score
                self.max_score += sph_max_score
                
                if self.nacyl:
                    
                    fa_score, fa_max_score = self.fatty_acyl(chains[0][1])
                    self.score += fa_score
                    self.max_score += fa_max_score
                
                yield chains
                
                self.score -= sph_score
                self.max_score -= sph_max_score
                
                if self.nacyl:
                    
                    self.score -= fa_score
                    self.max_score -= fa_max_score
    
    def fatty_acyl(self, fa):
        
        score = 0
        max_score = 0
        
        if fa.attr.oh:
            
            if self.add == '[M-H2O+H]+':
                
                score -= 20
        
        return score, max_score
    
    def cer(self):
        
        
        score = 0
        max_score = 0
        
        non_hex_score, non_hex_max_score = self.non_hex()
        
        score += non_hex_score
        max_score += non_hex_max_score
        
        return score, max_score
    
    def non_hex(self):
        
        score = 0
        max_score = 0
        
        score -= sum(map(bool,
            (
                self.scn.has_fragment('NL [Hexose-H2O] (NL 162.05)'),
                self.scn.has_fragment('NL [Hexose] (NL 180.06)'),
                self.scn.has_fragment('NL [Hexose+H2O] (NL 198.07)'),
                self.scn.has_fragment('NL [2xHexose] (NL 342.1162)'),
                self.scn.has_fragment('NL [2xHexose+H2O] (NL 360.1268)'),
                self.scn.has_fragment('NL [2xHexose-H2O] (NL 324.1056)'),
                self.scn.has_fragment('NL [2xHexose+O] (NL 358.1111)'),
                self.scn.has_fragment('NL [2xHexose+C] (NL 372.1268)'),
                self.scn.has_fragment('NL [S] (NL 79.9568)'),
                self.scn.has_fragment('NL [S+H2O] (97.9674)'),
                self.scn.has_fragment('NL [Hexose+SO3] (NL 242.100)'),
                self.scn.has_fragment('NL [Hexose+SO3+H2O] (NL 260.0202)'),
                self.scn.has_fragment('NL [Hexose+SO3+2xH2O] (NL 278.0308)'),
                self.scn.has_fragment('NL [2xHexose+SO3] (NL 404.0625)'),
                self.scn.has_fragment('NL [2xHexose+SO3+H2O] (NL 422.0730)'),
                self.scn.has_fragment('NL [2xHexose+SO3+2xH2O] (NL 440.0836)'),
            )
        )) * 5
        
        return score, max_score
    
    def sm(self):
        
        score = 0
        max_score = 47
        
        if self.scn.most_abundant_fragment_is('PC/SM [P+Ch] (184.0733)'):
            
            score += 15
            
            score += sum(map(bool,
                (
                    self.scn.has_fragment('PC/SM [N+3xCH3] (60.0808)'),
                    self.scn.has_fragment('PC/SM [Ch] (86.096)'),
                    self.scn.has_fragment('PC/SM [Ch+H2O] (104.107)'),
                    self.scn.has_fragment('PC/SM [P+Et] (124.9998)'),
                    self.scn.has_fragment('PC/SM [Ch-Et] (58.0651)'),
                    self.scn.has_fragment('NL PC/SM [P+Ch] (NL 183.066)'),
                    self.scn.has_fragment('NL SM [P+Ch] (NL 201.0766)'),
                    self.scn.has_fragment('NL SM [N+3xCH3] (77.0841)'),
                    self.scn.has_fragment('NL [H2O] (NL 18.0106)'),
                )
            )) * 3
            
            if self.scn.has_chain_fragment_type(frag_type = 'Sph-2xH2O+H'):
                
                score += 5
            
            self.must_have_chains = False
        
        return score, max_score
    
    def pc(self):
        """
        Lyso-SM aka Sph-PC.
        
        Scherer 2010, Table 1.
        """
        
        score = 0
        max_score = 15
        
        if self.scn.most_abundant_fragment_is('PC/SM [P+Ch] (184.0733)'):
            
            score += 15
        
        return score, max_score
    
    def pe_cer(self):
        
        score = 0
        max_score = 30
        
        if self.scn.has_fragment('NL PE [P+E] (NL 141.0191)'):
            
            score += 15
            
            score += sum(map(bool,
                (
                    self.scn.has_fragment('PE [P+E] (142.0264)'),
                    self.scn.has_fragment('NL PE [P+E+H2O] (NL 159.0297)'),
                    self.scn.has_fragment('NL PE [P+E-H2O] (NL 123.0085)'),
                )
            )) * 5
        
        return score, max_score
    
    def cer1p(self):
        
        score = 0
        max_score = 31
        
        if self.scn.has_fragment('NL [P+H2O] (NL 115.9875)'):
            
            score += 10
        
        if self.scn.chain_among_most_abundant(3, frag_type = 'Sph-2xH2O-H'):
            
            score += 10
        
        if self.scn.has_chain_combination(
            self.rec,
            chain_param = (
                {'frag_type': {
                        'Sph-2xH2O+H',
                        'Sph-H2O+H',
                        'Sph-H2O-H'
                    }
                },
                {'frag_type': 'FA+NH+C2H2-OH'},
            )
        ):
            
            score += 5
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [P] (NL 79.9663)'),
                self.scn.has_fragment('NL [P] (NL 97.9769)'),
            )
        )) * 3
        
        non_hex_score, non_hex_max_score = self.non_hex()
        
        score += non_hex_score
        
        return score, max_score
    
    def hexcer(self):
        
        score = 0
        max_score = 14
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [Hexose-H2O] (NL 162.05)'),
                self.scn.has_fragment('NL [Hexose] (NL 180.06)'),
                self.scn.has_fragment('NL [Hexose+H2O] (NL 198.07)'),
            )
        )) * 3
        
        if self.hexcer_chain_combination():
            
            score += 5
        
        return score, max_score
    
    def hex2cer(self):
        
        score = 0
        max_score = 39
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [2xHexose] (NL 342.1162)'),
                self.scn.has_fragment('NL [2xHexose+H2O] (NL 360.1268)'),
            )
        )) * 10
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [2xHexose-H2O] (NL 324.1056)'),
                self.scn.has_fragment('NL [2xHexose+O] (NL 358.1111)'),
                self.scn.has_fragment('NL [2xHexose+C] (NL 372.1268)'),
            )
        )) * 3
        
        if self.hexcer_chain_combination():
            
            score += 10
        
        return score, max_score
    
    def hexcer_chain_combination(self):
        
        return self.scn.has_chain_combination(
            self.rec,
            chain_param = (
                {'frag_type': {
                        'Sph-2xH2O+H',
                        'Sph-2xH2O-H',
                        'Sph-H2O+H',
                        'Sph-H2O-H',
                        'Sph-C-2xH2O',
                    }
                },
                {'frag_type': {
                        'FA-OH',
                        'NL FA',
                        'FA+NH2-O',
                    }
                },
            )
        )
    
    def shexcer(self):
        
        score = 0
        max_score = 25
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [S] (NL 79.9568)'),
                self.scn.has_fragment('NL [S+H2O] (97.9674)'),
                self.scn.has_fragment('NL [Hexose+SO3] (NL 242.100)'),
                self.scn.has_fragment('NL [Hexose+SO3+H2O] (NL 260.0202)'),
                self.scn.has_fragment('NL [Hexose+SO3+2xH2O] (NL 278.0308)'),
            )
        )) * 5
        
        return score, max_score
    
    def shex2cer(self):
        
        score = 0
        max_score = 25
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('NL [S] (NL 79.9568)'),
                self.scn.has_fragment('NL [S+H2O] (97.9674)'),
                self.scn.has_fragment('NL [2xHexose+SO3] (NL 404.0625)'),
                self.scn.has_fragment('NL [2xHexose+SO3+H2O] (NL 422.0730)'),
                self.scn.has_fragment('NL [2xHexose+SO3+2xH2O] (NL 440.0836)'),
            )
        )) * 5
        
        return score, max_score
    
    def m2(self):
        
        score = 0
        max_score = 46
        
        if self.scn.has_fragment('PC/SM [Ch-Et] (58.0651)'):
            
            score += 10
        
        if (
            self.rec.chainsum.u > 0 and
            self.scn.has_fragment('[C7+NH2] (110.0964)')
        ):
            
            score += 10
        
        score += sum(map(bool, (
            self.scn.has_fragment('[C5+NH2+2H] (84.0808)'),
            self.scn.has_fragment('[C6+NH2] (96.0808)'),
        ))) * 3
        
        if self.scn.has_chain_combination(
                self.rec,
                chain_param = (
                    {
                        'frag_type': {
                            'Sph-2xH2O+CH3',
                            'Sph-O-H2O+CH3+H',
                            'Sph-2xH2O+2xCh3+H',
                            'Sph-H2O+2xCH3+H',
                        }
                    },
                )
            ):
            
            score += 20
        
        return score, max_score
    
    def m1(self):
        
        score = 0
        max_score = 30
        
        if self.scn.has_fragment('PC/SM [Ch-Et] (58.0651)'):
            
            score += 10
        
        if self.scn.has_chain_combination(
                self.rec,
                head = 10,
                chain_param = (
                    {
                        'frag_type': {
                            'Sph-2xH2O+CH3',
                            'Sph-O-H2O+CH3+H',
                            'Sph-H2O+CH3+H',
                        }
                    },
                )
            ):
            
            score += 20
        
        if self.scn.has_chain_fragment_type(
                frag_type = {'Sph-2xH2O+2xCH3+H', 'Sph-H2O+2xCH3+H'},
                c = self.rec.chainsum.c - 1,
                u = self.rec.chainsum.u,
            ):
            
            score -= 20
        
        return score, max_score
    
    def m3(self):
        
        score = 0
        max_score = 20
        
        if self.scn.fragment_among_most_abundant(
            3, 'PC/SM [N+3xCH3] (60.0808)'
        ):
            
            score += 20
        
        return score, max_score
    
    def sph(self):
        
        score = 0
        max_score = 9
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('[C3+NH2] (56.0495)'),
                self.scn.has_fragment('[C2+NH2+O] (60.0444)'),
                self.scn.has_fragment('[C4+NH2+OH] (86.0600)'),
            )
        )) * 3
        
        return score, max_score
    
    def sphingosine_base(self, sph):
        
        if sph not in self.sph_scores:
            
            method = 'sphingosine_%s' % sph.lower()
            
            self.sph_scores[sph], self.sph_max_scores[sph] = (
                getattr(self, method)() if hasattr(self, method) else (0, 0)
            )
        
        return self.sph_scores[sph], self.sph_max_scores[sph]
    
    def sphingosine_d(self):
        
        score = 0
        max_score = 20
        
        if self.rec.chainsum and self.rec.chainsum.u == 0:
            
            return score, max_score
        
        if (
            self.nacyl and self.scn.chain_fragment_type_is(
                0,
                frag_type = 'Sph-2xH2O+H',
                u = (False, {0}),
            )
        ) or (
            not self.nacyl and self.scn.chain_fragment_type_is(
                0,
                frag_type = 'Sph-H2O+H',
                u = (False, {0}),
            )
        ) or (
            self.rec.hg.main == 'SM' and self.scn.chain_among_most_abundant(
                5,
                frag_type = 'Sph-2xH2O+H',
                u = (False, {0}),
                skip_non_chains = True,
            )
        ):
            
            score += 6
            
            score += sum(map(bool,
                (
                    not self.scn.has_fragment('[C2+NH2+O] (60.0444)'),
                    self.scn.has_fragment('NL [C+2xH2O] (NL 48.0211)')
                )
            )) * 2
            
            if (
                self.scn.chain_fragment_type_among_most_abundant(
                    4, frag_type = 'Sph-H2O+H', u = (False, {0})
                ) and
                self.scn.chain_fragment_type_among_most_abundant(
                    4, frag_type = 'Sph-C-O-H2O-H', u = (False, {0})
                ) and
                self.scn.chain_fragment_type_among_most_abundant(
                    4, frag_type = 'Sph-2xH2O+H', u = (False, {0})
                )
            ):
                score += 10
        
        return score, max_score
    
    def sphingosine_dh(self):
        
        score = 0
        max_score = 20
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('[C2+NH2+O] (60.0444)'),
                not self.scn.has_fragment('NL [C+2xH2O] (NL 48.0211)')
            )
        ))
        
        score += sum(map(bool,
            (
                self.scn.chain_fragment_type_among_most_abundant(
                    5, frag_type = 'Sph-H2O+H', u = 0
                ),
                self.scn.chain_fragment_type_among_most_abundant(
                    5, frag_type = 'FA+NH2-O', u = 0
                ),
                self.scn.chain_fragment_type_among_most_abundant(
                    10, frag_type = 'Sph-2xH2O+H', u = 0
                ),
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph-C-O-H2O-H', u = 0
                ),
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph+H', u = 0
                ),
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph-C-O-H2O-NH', u = 0
                )
            )
        )) * 3
        
        return score, max_score
    
    def sphingosine_t(self):
        
        score = 0
        max_score = 20
        
        if all((
            self.scn.chain_fragment_type_among_most_abundant(
                5, frag_type = 'Sph-H2O-H',
            ),
            self.scn.chain_fragment_type_among_most_abundant(
                10, frag_type = 'Sph-2xH2O-H',
            ),
            (
                self.scn.fragment_among_most_abundant(
                    '[C2+NH2+O] (60.0444)'
                ) or
                self.rec.hg.sub
            )
        )):
            
            score = 9
            
            score += sum(map(bool,
                (
                    not self.scn.has_fragment('NL [C+2xH2O] (NL 48.0211)'),
                    self.scn.has_fragment('NL [3xH2O] (NL 54.0317)')
                )
            ))
            
            score += sum(map(bool,
                (
                    self.scn.has_chain_fragment_type(
                        frag_type = 'Sph-C-2xH2O',
                    ),
                    self.scn.has_chain_fragment_type(
                        frag_type = 'Sph+H2O-H',
                    ),
                    self.scn.chain_fragment_type_among_most_abundant(
                        5, frag_type = 'Sph-H',
                    ),
                )
            )) * 3
        
        return score, max_score
    
    def sphingosine_k(self):
        
        score = 0
        max_score = 39
        
        if self.scn.has_chain_fragment_type(frag_type = 'Sph-NH2-H2O-2H'):
            
            score += 15
        
        score += sum(map(bool,
            (
                self.scn.has_fragment('[C2+NH2+O] (60.0444)'),
                self.scn.has_fragment('[C4+NH2+OH] (86.0600)'),
                self.scn.has_fragment('[C6+OH] (99.0804)'),
                self.scn.has_fragment('[C3+NH2] (56.0495)'),
            )
        )) * 3
        
        score += sum(map(bool,
            (
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph-C-2xH2O',
                ),
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph-H2O-H',
                ),
                self.scn.has_chain_fragment_type(
                    frag_type = 'Sph-H',
                ),
                self.scn.chain_fragment_type_among_most_abundant(
                    5, frag_type = 'Sph-H',
                ),
            )
        )) * 3
        
        return score, max_score


class Cer_Negative(AbstractMS2Identifier):
    """
    Examines if a positive mode MS2 spectrum is a ceramide.
    Identifies ceramide varieties including sphingomyeline,
    ceramide-1-phosphare, ceramide-phosphoethanolamine,
    OH-acyl-ceramide, hexosyl and dihexosyl-ceramides,
    and d, t and DH long chain base varieties.
    
    dCer
    ====
    
    **Specimen:**
    
    - in vivo SEC14L1 583, 554, 580 (formiate adduct)
    - in vivo STARD11 583, 554 (formiate adduct)
    - standards
    
    DHCer
    =====
    
    **Specimen:**
    
    - standards
    
    tCer
    ====
    
    **Specimen:**
    
    - standards
    
    """
    
    class_methods = {
        'Cer': 'cer',
        'SM':  'sm',
    }
    
    subclass_methods = {
        '1P': 'cer1p',
        'Hex': 'hexcer',
    }
    
    def __init__(self, record, scan, **kwargs):
        
        self.nacyl = record.chainsum is not None and len(record.chainsum) > 1
        self.oacyl = record.chainsum is not None and len(record.chainsum) > 2
        
        AbstractMS2Identifier.__init__(
            self,
            record,
            scan,
            missing_chains = (),
            chain_comb_args = {},
            must_have_chains = True,
            **kwargs,
        )
        
        self.sph_scores     = {}
        self.fa_scores      = {}
        self.sph_max_scores = {}
        self.fa_max_scores  = {}
    
    def confirm_class(self):
        
        AbstractMS2Identifier.confirm_class(self)
        
        self.confirm_subclass()
    
    def confirm_chains_explicit(self):
        
        for chains in itertools.chain(
            AbstractMS2Identifier.confirm_chains_explicit(self),
            AbstractMS2Identifier.confirm_chains_implicit(self),
        ):
            
            if chains[0][0].attr.sph == self.rec.chainsum.attr[0].sph:
                
                # the sphingosin base and fatty acyl related part of the
                # score is valid only for the current chain combination
                # hence now we add these to the overall score, yield the
                # identification and then subtract them from the score
                sph_score, sph_max_score = self.sphingosine_base(
                    chains[0][0].attr.sph
                )
                self.score += sph_score
                self.max_score += sph_max_score
                
                if self.nacyl:
                    
                    fa_score, fa_max_score = self.fatty_acyl(chains[0][1])
                    self.score += fa_score
                    self.max_score += fa_max_score
                
                yield chains
                
                self.score -= sph_score
                self.max_score -= sph_max_score
                
                if self.nacyl:
                    
                    self.score -= fa_score
                    self.max_score -= fa_max_score
    
    def fatty_acyl(self, fa):
        
        score = 0
        max_score = 0
        
        if len(fa.attr.oh) == 1:
            
            max_score = 30
            
            if self.scn.has_chain_combination(
                self.rec,
                head = 20, # to exclude tCer
                chain_param = (
                    {
                        'frag_type': {
                            'Sph-H', # b1
                            'Sph-C2H4-NH2-H2O', # b5
                        }
                    },
                    {
                        'frag_type': {
                            'FA+C2+NH2+O',  # a5 @ hydroxyacyl
                            'FA+CH2+NH2+O', # a1 @ hydroxyacyl
                        }
                    }
                )
            ):
                
                score += 30
        
        return score, max_score
    
    def cer(self):
        
        max_score = 23
        
        cer_nl = (
            'NL H2O (NL 18.0106)', # Hsu c1
            'NL 2xH2O (NL 36.0211)', # Hsu c4
            'NL C+H2O (NL 30.0106)', # Hsu c2
            'NL CH2+H2O (NL 32.0262)', # Hsu c3
            'NL C+2xH2O (NL 48.0211)', # Hsu c5
            'NL C+3xH2O (66.0455)', # Hsu c6
        )
        
        score = sum(
            self.scn.has_fragment(frag_name, adduct = self.add)
            for frag_name in cer_nl
        ) * 3
        
        if self.scn.has_chain_combinations(self.rec, adduct = self.add):
            
            score += 5
        
        return score, max_score
    
    def sphingosine_base(self, sph):
        
        if sph not in self.sph_scores:
            
            method = 'sphingosine_%s' % sph.lower()
            
            self.sph_scores[sph], self.sph_max_scores[sph] = (
                getattr(self, method)() if hasattr(self, method) else (0, 0)
            )
        
        return self.sph_scores[sph], self.sph_max_scores[sph]
    
    def sphingosine_d_dh(self):
        
        score = 0
        max_score = 20
        
        if self.scn.has_chain_combination(
            self.rec,
            chain_param = (
                {
                    'frag_type': {
                        'Sph-H', # b1
                        'Sph-C2H4-3H', # b2
                        'Sph-CH2-H2O-H', # b3
                        'Sph-H2O-NH2-2H', # b4
                        'Sph-C2H4-NH2-H2O', # b5
                    }
                },
                {
                    'frag_type': {
                        'FA+C2+NH2', # a2
                        'FA+C2+NH2-O', # a3
                    }
                }
            )
        ):
            
            score += 20
        
        return score, max_score
    
    def sphingosine_d(self):
        
        score = 0
        max_score = 20
        
        d_dh_score, d_dh_max_score = self.sphingosine_d_dh()
        
        score += d_dh_score
        max_score += d_dh_max_score
        
        if self.scn.has_fragment('NL C+H2O (NL 30.0106)', adduct = self.add):
            
            score += 20
        
        return score, max_score
    
    def sphingosine_dh(self):
        
        score = 0
        max_score = -20
        
        d_dh_score, d_dh_max_score = self.sphingosine_d_dh()
        
        score += d_dh_score
        max_score += d_dh_max_score
        
        if self.scn.has_fragment('NL C+H2O (NL 30.0106)', adduct = self.add):
            
            score -= 20
        
        return score, max_score
    
    def sphingosine_t(self):
        
        score = 0
        max_score = 28
        
        if self.scn.has_fragment('NL C+3xH2O (66.0455)', adduct = self.add):
            
            score += 5
        
        if self.scn.has_fragment('HexCer identity II'):
            
            score += 3
        
        if self.scn.has_chain_combination(
            self.rec,
            chain_param = (
                {
                    'frag_type': {
                        'Sph-CH2-NH2-4H', # b6
                    }
                },
                {
                    'frag_type': {
                        'FA+C2H2+NH2', # a1
                        'FA+C3H2+NH2', # a10
                    }
                }
            )
        ):
            
            score += 20
        
        # differentiate from hydroxyacyl-dCer
        if self.scn.chain_percent_of_most_abundant(
            frag_type = {'FA+C2H2+NH2+O', 'FA+C2+NH2+O'},
            percent = 5.0,
        ):
            
            score -= 10
        
        return score, max_score
    
    def cer1p(self):
        
        score = 0
        max_score = 70
        
        if any(map(bool, (
            self.scn.has_fragment('Cer1P/PIP/PL metaphosphate (78.9591)'),
            self.scn.has_fragment('Cer1P/PI phosphate (96.9696)'),
        ))):
            
            score += 20
            
            if self.scn.has_fragment(
                'NL H2O (NL 18.0106)', adduct = self.add
            ):
                
                score += 10
            
            if self.scn.has_chain_fragment_type(
                frag_type = {'NLFA_pH2O', 'NLFA_p2xH2O'},
                adduct = self.add
            ):
                
                score += 10
            
            self.must_have_chains = False
        
        if (
            self.rec.hg.main == 'Sph' and
            self.scn.fragment_among_most_abundant(
                'Cer1P/PIP/PL metaphosphate (78.9591)', 3
            )
        ):
            
            score += 20
            self.must_have_chains = False
            
            if self.scn.has_fragment('Cer1P/PI phosphate (96.9696)'):
                
                score += 10
        
        return score, max_score
    
    def sm(self):
        
        score = 0
        max_score = 45
        
        if self.scn.fragment_among_most_abundant(
            'NL CH2 (NL 14.0157)', 3, adduct = self.add
        ) and self.scn.fragment_among_most_abundant(
            'PC/SM PO4+choline-CH3 (168.0431)', 5
        ):
            
            score += 30
            
            score += sum(map(bool, (
                self.scn.has_fragment(
                    'Cer1P/PIP/PL metaphosphate (78.9591)'
                ),
                self.scn.has_fragment(
                    'NL choline+H2O', adduct = self.add
                ),
                self.scn.has_fragment(
                    'NL choline+H2O-CH3', adduct = self.add
                ),
            ))) * 5
            
            self.must_have_chains = False
        
        return score, max_score
    
    def hexcer(self):
        
        score = 0
        max_score = 90
        
        self.score += sum(map(bool, (
            self.scn.fragment_among_most_abundant('HexCer identity I', 10),
            self.scn.fragment_among_most_abundant('HexCer identity II', 10),
            self.scn.fragment_among_most_abundant('HexCer identity III', 10),
            self.scn.has_fragment('[Hexose] (179.0561)'),
            self.scn.has_fragment('[Hexose-H2O] (161.0455)'),
            self.scn.has_fragment('[Hexose-HCHO] (149.0455)'),
            self.scn.has_fragment('NL hexose (162.053)'),
            self.scn.has_fragment('NL hexose+H2O (180.063)'),
        ))) * 10
        
        if self.scn.has_chain_combinations(self.rec):
            
            score += 10
        
        return score, max_score
    
    def hex2cer(self):
        
        score = 0
        max_score = 140
        
        self.score += sum(map(bool, (
            self.scn.has_fragment('HexCer identity I'),
            self.scn.has_fragment('HexCer identity II'),
            self.scn.has_fragment('HexCer identity III'),
            self.scn.has_fragment('[Hexose] (179.0561)'),
            self.scn.has_fragment('[Hexose-H2O] (161.0455)'),
            self.scn.has_fragment('[Hexose-HCHO] (149.0455)'),
            self.scn.has_fragment('NL hexose (162.053)'),
            self.scn.has_fragment('NL hexose+H2O (180.063)'),
            self.scn.has_fragment('[2xHexose-HCHO] (311.0984)'),
            self.scn.has_fragment('[2xHexose-H2O] (323.0984)'),
            self.scn.has_fragment('[2xHexose] (341.1089)'),
            self.scn.has_fragment('NL 2xHexose (324.106)'),
            self.scn.has_fragment('NL 2xHexose+H2O (342.1162)'),
        ))) * 10
        
        if self.scn.has_chain_combinations(self.rec):
            
            score += 10
        
        return score, max_score
    
    def shexcer(self):
        
        score = 0
        max_score = 100
        
        self.missing_chains = (1,)
        
        if self.scn.has_fragment('Sulphate (96.9601)'):
            
            score += 20
        
        self.score += sum(map(bool, (
            self.scn.has_fragment('[Sulfohexose] (259.0129)'),
            self.scn.has_fragment('[Sulfohexose] (256.9972)'),
            self.scn.has_fragment('[Sulfohexose-H2O] (241.0024)'),
            self.scn.has_fragment('[Sulfohexose+Et+N] (300.0395)'),
        ))) * 10
        
        if self.scn.has_chain_fragment_type(
            frag_type = {
                'Sph+C6O5H8+SO3+H2O',
                'Sph+C6O5H8+SO3+CO+H2O',
            }
        ):
            
            score += 20
            
            if self.scn.has_chain_combination(
                self.rec,
                chain_param = (
                    {
                        'frag_type': {
                            'Sph+C6O5H8+SO3',
                            'Sph+C6O5H8+SO3+H2O',
                            'Sph+C6O5H8+SO3+CO+H2O',
                        }
                    },
                    {
                        'frag_type': {
                            'NLFA',
                            'NLFA_mH2O',
                        }
                    }
                )
            ):
                
                score += 20
        
        return score, max_score
    
    def shex2cer(self):
        
        score = 0
        max_score = 140
        
        self.missing_chains = (1,)
        
        if self.scn.has_fragment('Sulphate (96.9601)'):
            
            score += 20
        
        score += sum(map(bool, (
            self.scn.has_fragment('[Sulfohexose] (259.0129)'),
            self.scn.has_fragment('[Sulfohexose] (256.9972)'),
            self.scn.has_fragment('[Sulfohexose-H2O] (241.0024)'),
            self.scn.has_fragment('[Sulfohexose+Et+N] (300.0395)'),
            self.scn.has_fragment('[2xHexose-H2O+SO3] (403.0552)'),
            self.scn.has_fragment('[2xHexose+SO3] (419.0501)'),
            self.scn.has_fragment('[2xHexose+SO3] (421.0658)'),
            self.scn.has_fragment('[2xHexose+SO3+Et+N] (462.0923)'),
        ))) * 10
        
        if self.scn.has_chain_fragment_type(
            frag_type = {
                'Sph+C12O10H18+SO3',
                'Sph+C12O10H18+SO3+H2O',
                'Sph+C12O10H18+SO3+CO+H2O',
            }
        ):
            
            score += 20
            
            if self.scn.has_chain_combination(
                self.rec,
                chain_param = (
                    {
                        'frag_type': {
                            'Sph+C12O10H18+SO3',
                            'Sph+C12O10H18+SO3+H2O',
                            'Sph+C12O10H18+SO3+CO+H2O',
                        }
                    },
                    {
                        'frag_type': {
                            'NLFA',
                            'NLFA_mH2O',
                        }
                    }
                )
            ):
                
                score += 20
        
        return score, max_score
    
    def pe_cer(self):
        
        score = 0
        max_score = 30
        
        score += sum(map(bool, (
            self.scn.has_fragment('PE [P+E] (140.0118)'),
            self.scn.has_fragment('NL PE [P+E] (141.0191)'),
            self.scn.has_fragment('PE [P+E-H2O] (122.0013)'),
        ))) * 10
        
        return score, max_score

#
# Scan.identify() dispatches identification methods as below
#

idmethods = {
    'neg': {
        lipproc.Headgroup(main = 'FA'):  FA_Negative,
        lipproc.Headgroup(main = 'DAG'): DAG_Negative,
        lipproc.Headgroup(main = 'TAG'): TAG_Negative,
        lipproc.Headgroup(main = 'DGTA'): GL_Negative,
        lipproc.Headgroup(main = 'DGTS'): GL_Negative,
        lipproc.Headgroup(main = 'DGCC'): GL_Negative,
        lipproc.Headgroup(main = 'SQDG'): GL_Negative,
        lipproc.Headgroup(main = 'MGDG'): GL_Negative,
        lipproc.Headgroup(main = 'DGDG'): GL_Negative,
        lipproc.Headgroup(main = 'DGTA', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'DGTS', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'DGCC', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'SQDG', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'MGDG', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'DGDG', sub = ('Lyso',)): GL_Negative,
        lipproc.Headgroup(main = 'PE'):  PE_Negative,
        lipproc.Headgroup(main = 'PE', sub = ('Lyso',)): PE_Negative,
        lipproc.Headgroup(main = 'PC'):  PC_Negative,
        lipproc.Headgroup(main = 'PC', sub = ('Lyso',)): PC_Negative,
        lipproc.Headgroup(main = 'PI'):  PI_Negative,
        lipproc.Headgroup(main = 'PI', sub = ('Lyso',)):  PI_Negative,
        lipproc.Headgroup(main = 'PS'):  PS_Negative,
        lipproc.Headgroup(main = 'PS', sub = ('Lyso',)): PS_Negative,
        lipproc.Headgroup(main = 'PG'):  PG_Negative,
        lipproc.Headgroup(main = 'PG', sub = ('Lyso',)):  PG_Negative,
        lipproc.Headgroup(main = 'BMP'): BMP_Negative,
        lipproc.Headgroup(main = 'PA'):  PA_Negative,
        lipproc.Headgroup(main = 'PA', sub = ('Lyso',)):  PA_Negative,
        lipproc.Headgroup(main = 'VA'): VA_Negative,
        lipproc.Headgroup(main = 'Cer'): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('1P',)): Cer_Negative,
        lipproc.Headgroup(main = 'SM'): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex',)): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex2',)): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex',)): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex2',)): Cer_Negative,
        lipproc.Headgroup(main = 'Cer', sub = ('PE',)): Cer_Negative,
        lipproc.Headgroup(main = 'Sph'): Cer_Negative,
        lipproc.Headgroup(main = 'Sph', sub = ('1P',)): Cer_Negative,
    },
    'pos': {
        lipproc.Headgroup(main = 'FA'):  FA_Positive,
        lipproc.Headgroup(main = 'DAG'): DAG_Positive,
        lipproc.Headgroup(main = 'DGTA'): GL_Positive,
        lipproc.Headgroup(main = 'DGTS'): GL_Positive,
        lipproc.Headgroup(main = 'DGCC'): GL_Positive,
        lipproc.Headgroup(main = 'SQDG'): GL_Positive,
        lipproc.Headgroup(main = 'MGDG'): GL_Positive,
        lipproc.Headgroup(main = 'DGDG'): GL_Positive,
        lipproc.Headgroup(main = 'DGTA', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'DGTS', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'DGCC', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'SQDG', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'MGDG', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'DGDG', sub = ('Lyso',)): GL_Positive,
        lipproc.Headgroup(main = 'TAG'): TAG_Positive,
        lipproc.Headgroup(main = 'PE'):  PE_Positive,
        lipproc.Headgroup(main = 'PE', sub = ('Lyso',)):  LysoPE_Positive,
        lipproc.Headgroup(main = 'PC'):  PC_Positive,
        lipproc.Headgroup(main = 'PC', sub = ('Lyso',)):  LysoPC_Positive,
        lipproc.Headgroup(main = 'PI'):  PI_Positive,
        lipproc.Headgroup(main = 'PI', sub = ('Lyso',)):  PI_Positive,
        lipproc.Headgroup(main = 'PS'):  PS_Positive,
        lipproc.Headgroup(main = 'PS', sub = ('Lyso',)):  PS_Positive,
        lipproc.Headgroup(main = 'PG'):  PG_Positive,
        lipproc.Headgroup(main = 'PG', sub = ('Lyso',)):  PG_Positive,
        lipproc.Headgroup(main = 'BMP'): BMP_Positive,
        lipproc.Headgroup(main = 'PA'):  PA_Positive,
        lipproc.Headgroup(main = 'PA', sub = ('Lyso',)):  PA_Positive,
        lipproc.Headgroup(main = 'VA'): VA_Positive,
        lipproc.Headgroup(main = 'Cer'): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('1P',)): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex',)): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex2',)): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex',)): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex2',)): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('PE',)): Cer_Positive,
        lipproc.Headgroup(main = 'SM'): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('1P', 'Lyso')): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex', 'Lyso')): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('Hex2', 'Lyso')): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex', 'Lyso')): Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('SHex2', 'Lyso')):
            Cer_Positive,
        lipproc.Headgroup(main = 'Cer', sub = ('PE', 'Lyso')): Cer_Positive,
        lipproc.Headgroup(main = 'SM',  sub = ('Lyso',)): Cer_Positive,
        lipproc.Headgroup(main = 'Sph'): Cer_Positive,
        lipproc.Headgroup(main = 'Sph', sub = ('1P',)): Cer_Positive,
        lipproc.Headgroup(main = 'Sph', sub = ('M1',)): Cer_Positive,
        lipproc.Headgroup(main = 'Sph', sub = ('M2',)): Cer_Positive,
        lipproc.Headgroup(main = 'Sph', sub = ('M3',)): Cer_Positive,
    }
}


class MS2Feature(object):
    
    scan_methods = {
        'mgf':  'mgf_iterscans',
        'mzml': 'mzml_iterscans',
    }
    
    def __init__(
            self,
            mz,
            ionmode,
            resources,
            rt,
            ms1_records = None,
            rt_range = .5,
            check_rt = True,
            add_precursor_details = False,
        ):
        """
        Collects the MS2 scans from the provided resources for a single
        feature. Calls identification methods on all scans collected.
        
        :param float mz:
            m/z value of the precursor ion.
        :param str ionmode:
            Ion mode of the experiment. Either ``pos`` or ``neg``.
        :param dict resources:
            ``dict`` of MS2 scan resources. These are either ``mgf.MgfReader``
            objects or paths to MGF files. Later more resource types
            will be available, for example MzML format. Keys of the ``dict``
            are used as sample labels. Thes can be strings or tuples.
        :param dict ms1_records:
            A data structure resulted by ``moldb.adduct_lookup``. If ``None``
            the lookup will be done here.
        :param float rt_range:
            If a single retention time value provided this is the largest
            accepted difference between an MS2 scan's RT and the precursor's
            RT. E.g. if ``rt = 8.3`` and ``rt_range = 0.5``, scans between
            7.8 and 8.8 will be considered. If a tuple of floats provided
            for RT, scans between these two values will be considered.
        :param bool check_rt:
            Check if the retention time of the scan is enough close to the
            precursor's RT. If ``False``, scans will be matched only by the
            m/z value of the precursor and scans with any large RT difference
            will be analysed.
        """
        
        self.mz = mz
        self.ionmode = ionmode
        self.ms1_records = ms1_records or moldb.adduct_lookup(mz, ionmode)
        self.add_precursor_details = add_precursor_details
        self.resources = resources
        self.rt = (rt - rt_range, rt + rt_range) if type(rt) is float else rt
        self.rtmean = sum(self.rt) / 2.0
        self.rt_range = rt_range
        self.check_rt = check_rt
    
    def main(self):
        
        self.ms1_lookup()
        self.build_scans()
        self.identify()
    
    def iterscans(self):
        
        for sample_id, resource in iteritems(self.resources):
            
            res_type = self.guess_resouce_type(resource)
            
            if res_type not in self.scan_methods:
                
                raise ValueError(
                    'Unknown MS2 resource type: %s' % str(resource)
                )
            
            scan_method = getattr(self, self.scan_methods[res_type])
            
            for scan in scan_method(resource, sample_id):
                
                yield scan
    
    def mgf_iterscans(self, mgf_resource, sample_id = None):
        
        if isinstance(mgf_resource, basestring):
            
            mgffile = mgf.MgfReader(mgfname, charge = None)
            
        elif isinstance(mgf_resource, mgf.MgfReader):
            
            mgffile = mgf_resource
            
        else:
            
            raise ValueError(
                'Mgf files should be lipyd.mgf.MgfReader '
                'instances of file names.'
            )
        
        idx, rtdiff = mgffile.lookup(self.mz, rt = self.rtmean)
        
        for i, rtd in zip(idx, rtdiff):
            
            if self.check_rt:
                
                scan_rt = self.rtmean + rtd
                
                if scan_rt < self.rt[0] or scan_rt > self.rt[1]:
                    
                    continue
            
            sc = mgffile.get_scan(i)
            
            yield Scan(
                mzs = sc[:,0],
                intensities = sc[:,1],
                ionmode = self.ionmode,
                precursor = self.mz,
                ms1_records = self.ms1_records,
                add_precursor_details = self.add_precursor_details,
                scan_id = mgffile.mgfindex[i,3],
                sample_id = sample_id,
                source = mgffile.fname,
                deltart = rtd,
                rt = mgffile.mgfindex[i,2],
            )
    
    def mzml_iterscans(self, mzml_resource, sample_id = None):
        
        raise NotImplementedError
    
    @staticmethod
    def guess_resouce_type(res):
        
        if isinstance(res, basestring) and os.path.exists(res):
            
            if res[-3:].lower() == 'mgf':
                
                return 'mgf'
            
        elif isinstance(res, mgf.MgfReader):
            
            return 'mgf'
    
    def build_scans(self):
        
        self.scans = np.array(list(self.iterscans()))
        self.deltart = np.array([sc.rt - self.rtmean for sc in self.scans])
        
        rtsort = [
            it[0]
            for it in sorted(
                (it for it in enumerate(self.deltart)),
                key = lambda it: abs(it[1])
            )
        ]
        
        self.scans = self.scans[rtsort]
        self.deltart = self.deltart[rtsort]
    
    def identify(self):
        
        self.identities = []
        
        for scan in self.scans:
            
            identity = scan.identify()
            
            if identity:
                
                self.identities.append(identity)
    
    @staticmethod
    def identities_sort(ids):
        """
        Sorts identities by score and deltart.
        """
        
        return sorted(
            ids,
            key = lambda i:
                (
                    100 - i.score_pct,
                    # if we don't have scan details don't consider
                    i.scan_details.deltart if i.scan_details else 0,
                )
        )
    
    def identities_group_by(self, by = 'subspecies'):
        """
        Parameters
        ----------
        species : str
            Either ``species`` or ``subspecies`` or ``subclass``.
        """
        
        identities = {}
        
        for scan in self.identities:
            
            for sum_str, varieties in iteritems(scan):
                
                for var in varieties:
                    
                    key = (
                        sum_str
                        if by == 'species' else
                        lipproc.full_str(var.hg, var.chains)
                        if by == 'subspecies' else
                        lipproc.subclass_str(var.hg)
                        if by == 'subclass' else
                        lipproc.class_str(var.hg)
                    )
                    
                    if key not in identities:
                        
                        identities[key] = []
                    
                    identities[key].append(var)
        
        for k, v in iteritems(identities):
            
            identities[k] = self.identities_sort(v)
        
        return identities
    
    def identities_group_by_species(self):
        
        return self.identities_group_by(by = 'species')
    
    def identities_group_by_subspecies(self):
        
        return self.identities_group_by(by = 'subspecies')
    
    def identities_group_by_subclass(self):
        
        return self.identities_group_by(by = 'subclass')
    
    def identities_group_by_class(self):
        
        return self.identities_group_by(by = 'class')
    
    def identities_str(
            self,
            only_best = True,
            only_top = True,
            group_by = 'subspecies'
        ):
        """
        only_best : bool
            Only the ones with highest score.
        only_top : bool
            Only the first one for each group (the one with highest score
            and lowest delta RT), ignoring the repeated identifications.
        group_by : str
            Passed to ``identities_group_by``. Group by lipid identification
            level ``subspecies``, ``species``, ``subclass`` or ``class``.
        """
        
        identities = self.identities_group_by(by = group_by)
        
        if only_best:
            
            max_score = max(
                (ids[0].score_pct for ids in identities.values()),
                default = 0,
            )
            
        else:
            
            max_score = 0
        
        return ';'.join(
            (
                self._identities_str(
                    ids if only_top else (ids[0],),
                    sum_str = group_by != 'subspecies'
                )
                if ids else ''
            )
            for ids in identities.values()
            if (
                not only_best or (
                    ids[0].score_pct > 0 and
                    ids[0].score_pct == max_score
                )
            )
        )
    
    def identities_str_all(self):
        
        return self.identities_str(
            only_best = False,
            only_top  = False,
            group_by  = 'subspecies',
        )
    
    def identities_str_best_sum(self):
        
        return self.identities_str(
            only_best = True,
            only_top  = True,
            group_by  = 'species',
        )
    
    def identities_str_full_toplist(self):
        
        return self.identities_str(
            only_best = False,
            only_top  = True,
            group_by  = 'subclass',
        )
    
    def _identities_str(self, ids, sum_str = False):
        
        return (
            ';'.join(
                sorted(self.format_ms2id(i, sum_str = sum_str) for i in ids)
            )
        )
    
    def format_ms2id(self, i, sum_str = False):
        # TODO: more generic handling of sample IDs!!!
        
        adduct = (
            ',adduct=%s' % i.precursor_details.adduct
                if self.add_precursor_details else
            ''
        )
        error = (
            ',ms1ppm=%0.1f' % i.precursor_details.error
                if self.add_precursor_details else
            ''
        )
        
        return '%s[score=%u,deltart=%.02f,fraction=%s%u,scan=%u%s%s]' % (
            (
                lipproc.summary_str(i.hg, i.chainsum)
                if sum_str and i.chainsum else
                i.__str__()
            ),
            i.score_pct,
            i.scan_details.deltart,
            i.scan_details.sample_id[0],
            i.scan_details.sample_id[1],
            i.scan_details.scan_id,
            adduct,
            error,
        )
    
    def identity_summary(
            self,
            chains = True,
            scores = True,
            drt = True,
            sample_ids = False,
            scan_ids = False,
        ):
        
        identities = {}
        
        for i, scan_i in enumerate(self.identities):
            
            for sum_str, varieties in iteritems(scan_i):
                
                for var in varieties:
                    
                    key = [sum_str]
                    
                    if chains:
                        
                        key.append(var.__str__())
                        
                    if scores:
                        
                        key.append(var.score_pct)
                        
                    else:
                        
                        if var.score == 0:
                            
                            continue
                    
                    if drt:
                        
                        key.append(self.deltart[i])
                    
                    if sample_ids:
                        
                        key.append(self.scans[i].sample_id)
                    
                    if scan_ids:
                        
                        key.append(self.scans[i].scan_id)
                    
                    key = tuple(key)
                    
                    if key not in identities:
                        
                        identities[key] = []
                    
                    identities[key].append(var)
        
        return identities
    
    def ms1_lookup(self):
        
        if self.ms1_records is None:
            
            self.ms1_records = moldb.adduct_lookup(self.mz)
