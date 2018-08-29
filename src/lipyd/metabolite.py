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

from __future__ import generator_stop

from future.utils import iteritems
from past.builtins import xrange, range

import itertools
import collections
import operator
import functools
import copy

import numpy as np

import lipyd.formula as formula
import lipyd.lipproc as lipproc


class AbstractMetaboliteComponent(formula.Formula):
    
    def __init__(self,
                 core = 0.0,
                 charge = 0,
                 isotope = 0,
                 name = 'Unknown',
                 getname = lambda parent: parent.name,
                 **kwargs):
        """
        Represents a component of a molecule. It can be the basis of a core
        carrying more substituents. Or can be the basis of a substituent
        which might provide a custom set of options or a homolog series
        with different carbon counts and unsaturations on the aliphatic
        chain.
        
        :param core: You can provide the core (unchanged part) 3 ways:
            As a `float` which considered to be an exact mass.
            As a formula, e.g. `C2H5OH`.
            As a `dict` of atom counts, e.g. `{'C': 2, 'H': 6, 'O': 1}`.
        """
        
        formula.Formula.__init__(
            self,
            core if type(core) is not float else None,
            charge,
            isotope,
            **kwargs
        )
        
        if not self.has_mass():
            
            if type(core) is float:
                
                self.mass = core
                
            else:
                
                raise ValueError('Please provide either formula or '
                                'atom counts or mass.')
        
        self.name = name
        self.getname = getname


class AbstractMetabolite(AbstractMetaboliteComponent):
    
    def __init__(self,
            core = 0.0,
            charge = 0,
            isotope = 0,
            name = 'Unknown',
            hg = None,
            getname = lambda parent, subs:
                '%s(%s)' % (
                    parent.name,
                    '/'.join(
                        s.cc_unsat_str
                        for s in subs
                        if hasattr(s, 'cc_unsat_str') and s.cc_unsat_str
                    )
                ),
            subs = None,
            sum_only = False,
            **kwargs
        ):
        """
        Represents a metabolite with an unchanged core and a set of variable
        substituents.
        
        :param bool sum_only: Do not iterate all aliphatic chains
            independently but consider only sum of chain lengths and
            unsaturations.
        """
        
        AbstractMetaboliteComponent.__init__(
            self,
            core,
            charge = charge,
            isotope = isotope,
            name = name,
            getname = getname,
            **kwargs
        )
        
        self.subs = subs or ()
        self.subs = tuple(self.get_substituent(s) for s in self.subs)
        self.sum_only = sum_only
        self.sub0 = None
        self.hg = hg
    
    def __iter__(self):
        
        for subs, inst in self.subsproduct():
            
            yield inst
    
    @staticmethod
    def get_substituent(sub):
        """
        Creates a `Formula` object if the substituent is not
        already an instance of `Formula` or `Substituent`.
        """
        
        return (
                formula.Formula(sub)
            if hasattr(sub, 'lower') or type(sub) is float
                else formula.Formula(**sub)
            if type(sub) is dict
                else sub
        )
    
    def itersubs(self):
        """
        Iterates all combinations of all substituents.
        Yields tuples of substituents.
        """
        
        self._restore_sub0()
        
        for subs in itertools.product(*self.subs):
            
            yield subs
    
    def subsproduct(self):
        """
        Iterates instances and substituents in parallel.
        Yields tuples of two elements.
        First element is a tuple of all molecule parts (substituents).
        Second element is the actual instance, i.e. the whole molecule.
        """
        
        iterator = self.itersum() if self.sum_only else self.itersubs()
        
        for subs in iterator:
            
            self.inst_name = self.getname(self, subs)
            
            self.inst = functools.reduce(
                operator.add,
                itertools.chain((self,), (s for s in subs))
            )
            self.inst.name = self.inst_name
            
            yield subs, self.inst
    
    def iterlines(self):
        """
        Iterates standard lines.
        """
        
        for subs, inst in self.subsproduct():
            
            full_name = inst.getname()
            
            chains = tuple(
                s.attrs.chain for s in subs if hasattr(s.attrs, 'chain')
            )
            chainsum = lipproc.sum_chains(chains)
            name = (
                (lipproc.summary_str(self.hg, chainsum),)
                    if self.hg and self.sum_only else
                (lipproc.full_str(self.hg, chains),)
                    if self.hg and not self.sum_only else
                ()
            )
            lab = lipproc.LipidLabel(
                db_id = None,
                db = 'lipyd.lipid',
                names = name
            )
            rec = lipproc.LipidRecord(
                lab = lab,
                hg  = self.hg,
                chainsum = chainsum,
                chains = () if self.sum_only else chains
            )
            
            yield inst.mass, rec
    
    def itersum(self):
        """
        Iterates by considering only the sum of chain lengths and
        unsaturations.
        """
        
        self._restore_sub0()
        
        chains = [
            (i, s.chlens, s.unsats)
            for i, s in enumerate(self.subs)
            if self.has_variable_aliphatic_chain(s)
        ]
        
        if len(chains) <= 1:
            
            for subs in self.itersubs():
                
                yield subs
            
            return
        
        min_chlens = sum(min(c[1]) for c in chains[1:])
        min_unsats = sum(min(c[2]) for c in chains[1:])
        sum_chlens = list(set(
            sum(cc) - min_chlens for cc in
            itertools.product(*(c[1] for c in chains))
        ))
        sum_unsats = list(set(
            sum(uu) - min_unsats for uu in
            itertools.product(*(c[2] for c in chains))
        ))
        
        self.sub0 = chains[0]
        isub0 = self.sub0[0]
        sub0 = self.subs[isub0]
        
        sub0.chlens = sum_chlens
        sub0.unsats = sum_unsats
        
        for subs in itertools.product(*(
                s.__iter__()
                if i == isub0 else
                # other substituents iterated by their cores only
                s.__iter__(cores_only = True)
                for i, s in enumerate(self.subs)
            )):
            
            yield subs
        
        # at the end restore the real values
        self._restore_sub0()
    
    def _restore_sub0(self):
        """
        For iterating with considering only total carbon count and
        unsaturation accross all aliphatic chains, we pretend that
        all extra carbons and unsaturations are added to the first
        substituent with variable aliphatic chain, e.g. it may
        have C54 while in reality it has max C18 and the other 36
        carbons are in the other chains. After setting these fake
        numbers on the first substituent we need to restore the
        real numbers otherwise it could mess up things later.
        """
        
        if self.sub0:
            
            self.subs[self.sub0[0]].chlens = self.sub0[1]
            self.subs[self.sub0[0]].unsats = self.sub0[2]
    
    def has_variable_aliphatic_chain(self, sub):
        """
        Tells if the substituent really has more than one possible
        chain length or unsaturation variation.
        """
        
        return (
            hasattr(sub, 'variable_aliphatic_chain') and
            sub.variable_aliphatic_chain
        )


class AbstractSubstituent(AbstractMetaboliteComponent):
    
    def __init__(
            self,
            cores = 0.0,
            c = (14, 20),
            u = (0, 1),
            counts = None,
            charges = 0,
            isotopes = 0,
            names = 'Unknown',
            chain_attr = None,
            chain_type = None,
            getname = lambda parent: '%u:%u' % (parent.c, parent.u),
            c_u_diff = lambda c, u: c > u + 1 or (c == 0 and u == 0),
            prefix = '',
            even = False,
            valence = 1,
            **kwargs
        ):
        """
        Represents a set of distinct substituent parts in a molecule or
        a homolog series spanning accross a range of aliphatic chain
        length and unsaturated bonds.
        
        Args:
        -----
        :param list cores:
            List of core variations. Same kind of definitions
            are possible like at `formula.Formula`: formula as `str`, `dict` of
            atoms or exact mass as `float`.
        :param tuple c:
            Tuple of 2 integers: range of chain lengths.
        :param tuple u:
            Tuple of 2 integers: range of unsaturations.
        :param list counts:
            Dictionary with extra atom counts.
            If you have one or more extra oxygen, nitrogen, phosphorous or
            any other atoms in the compound you can include here.
            Alternatively you can also include them in the core.
            Also accounts for the valences of the aliphatic chain
            not occupied by hydrogens. E.g. for a fatty acyl you need to
            remove 3 hydrogens as 3 valences are occupied by the oxygens.
            If you have a secondary amine you need to remove one more
            hydrogen. Otherwise, as this data structure has no information
            about constitution, we could not guess the number of hydrogens.
            Similarly, for oxo groups removal of 2 hydrogens necessary.
        :param list charges:
            List of integers: charges for each core
            variation.
        :param list isotopes:
            List of integers: extra neutrons for each core
            variation.
        :param getname:
            Method (callable) to create a name from the chain
            length and the unsaturation.
        :param c_u_diff:
            Method (callable) to decide if the chain length and
            unsaturation are compatible. E.g. if chain length is only C2, an
            unsaturation of 4 double bonds is not possible. By default, chain
            length must be greater by 2 than unsaturation, it means at acyl
            chains we avoid to assume double bond right next to the carboxyl
            group which is clearly impossible.
        :param bool even:
            If true only even chain lengths are considered.
        :param int valence:
            Total valence of the substituent. Number of bonds connecting
            to other parts of the molecule.
        """
        
        self.prefix   = prefix
        self.c_u_diff = c_u_diff
        self.even     = even
        self.valence  = valence
        self.getname  = getname
        self.cores    = (
            cores if type(cores) in {list, set, tuple} else [cores]
        )
        self.set_attr(charges, 'charges')
        self.set_attr(names, 'names')
        self.set_attr(isotopes, 'isotopes')
        self.counts   = collections.defaultdict(lambda: 0)
        self.counts.update(counts or {})
        
        AbstractMetaboliteComponent.__init__(
            self,
            self.cores[0],
            charge = self.charges[0],
            isotope = self.isotopes[0],
            name = self.names[0],
            getname = getname,
            **kwargs
        )
        
        # range of possible lengths and unsats
        self.set_chlens(c)
        self.set_unsats(u)
        
        # current value of length and unsat
        self.c = self.chlens[0]
        self.u = self.unsats[0]
        
        self.total = len(self.chlens) * len(self.unsats)
        self.variable_aliphatic_chain = (
            not (
                len(self.chlens) == 1 and
                self.chlens[0] == 0
            )
        )
        
        self.chain_attr = chain_attr
        self.chain_type = chain_type
    
    def __iter__(self, cores_only = False):
        
        for i in range(len(self.cores)):
            
            self.update_core(i)
            
            for c in self.chlens:
                
                self.c = c
                
                for u in self.unsats:
                    
                    self.u = u
                    
                    if not self.c_u_diff(self.c, self.u):
                        
                        continue
                    
                    # implicit hydrogens
                    h = c * 2 + 2 - self.valence - 2 * u
                    p = self.get_prefix()
                    new_counts = self.counts.copy()
                    new_counts['C'] += c
                    new_counts['H'] += h
                    new_attrs = copy.deepcopy(self.attrs)
                    # `attrs` might contain methods which are called
                    # with the present instance passed and their returned
                    # value will be the attribute value of the
                    # iteration products
                    for k, v in iteritems(new_attrs.__dict__):
                        if hasattr(v, '__call__'):
                            setattr(new_attrs, k, v(self))
                    
                    if self.chain_type and self.chain_attr and c > 0:
                        
                        new_attrs.chain = self.get_chain()
                    
                    new = self + formula.Formula(**new_counts)
                    
                    new.attrs = new_attrs
                    new.c = c
                    new.u = u
                    new.get_prefix = lambda: p
                    new.variable_aliphatic_chain = (
                        self.variable_aliphatic_chain
                    )
                    # new.attrs.chain = self.get_chain()
                    
                    yield new
                    
                    if cores_only:
                        
                        break
                
                if cores_only:
                    
                    break
    
    def set_attr(self, val, name):
        
        setattr(
            self,
            name,
            val
                if type(val) in {list, set, tuple} else
            (val,) * len(self.cores)
        )
    
    def set_chlens(self, c):
        
        if type(c) is int: c = [c]
        if type(c) is tuple and len(c) != 2: c = list(c)
        
        self.chlens = (
            list(c) if type(c) in {list, set, range}
            else [
                i for i in
                range(c[0], c[1] + 1)
                if not self.even or i % 2 == 0
            ]
        )
    
    def set_unsats(self, u):
        
        if type(u) is int: u = [u]
        if type(u) is tuple and len(u) != 2: u = list(u)
        
        self.unsats = (
            list(u) if type(u) in {list, set, range}
            else list(range(u[0], u[1] + 1))
        )
    
    def update_core(self, i = 0):
        
        try:
            new = self.cores[i]
        
        except:
            raise IndexError('Number of cores is less than %u' % i + 1)
        
        if type(new) is float:
            
            self.mass = new
            self.formula = ''
            
        elif type(new) is dict:
            
            self.formula_from_dict(new)
            
        elif hasattr(new, 'lower'):
            
            self.formula = new
            
        else:
            
            raise ValueError('Wrong mass or formula: `%s`' % str(new))
        
        self.charge = (
            self.charges[i]
            if type(self.charges) in {list, set, tuple} else
            self.charges
        )
        self.isotope = (
            self.isotopes[i]
            if type(self.isotopes) in {list, set, tuple} else
            self.isotypes
        )
        self.name = (
            self.names[i]
            if type(self.names) in {list, set, tuple} else
            self.names
        )
        
        self.reset_atoms()
        self.calc_mass()
    
    def get_prefix(self):
        
        return self.prefix
    
    def get_chain(self):
        
        return lipproc.Chain(
            c = self.c,
            u = self.u,
            typ = self.chain_type,
            attr = self.chain_attr,
        )
