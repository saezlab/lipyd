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

from future.utils import iteritems
from past.builtins import xrange, range

import itertools
import collections
import copy

import lipyd.formula as formula


class AbstractMetaboliteComponent(formula.Formula):
    
    def __init__(self,
                 core = 0.0,
                 charge = 0,
                 isotope = 0,
                 name = 'Unknown',
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
        
        formula.Mass.__init__(self,
            core if type(core) is not float else None,
            charge,
            isotope,
            **kwargs)
        
        if not self.has_mass():
            
            if type(core) is float:
                
                self.mass = core
                
            else:
                
                raise ValueError('Please provide either formula or '
                                'atom counts or mass.')
        
        self.name = name


class AbstractMetabolite(AbstractMetaboliteComponent):
    
    def __init__(self, core = 0.0,
                 charge = 0,
                 isotope = 0,
                 name = 'Unknown',
                 syn = {},
                 subs = {},
                 **kwargs):
        """
        Represents a metabolite with an unchanged core and a set of variable
        substituents.
        """
        
        AbstractMetaboliteComponent.__init__(
            self,
            core,
            charge = charge,
            isotope = isotope,
            name = name,
            **kwargs
        )
        
        self.syn  = syn
        self.subs = subs
    
    def __iter__(self):
        
        for sub in itertools.product(self.subs.values()):
            
            yield sum(itertools.chain([self], (s for s in sub)))


class AbstractSubstituent(AbstractMetaboliteComponent):
    
    def __init__(
            self,
            cores = [0.0],
            c = (0, 1),
            u = (0, 1),
            counts = {},
            charges = [0],
            isotopes = [0],
            names = ['Unknown'],
            getname = lambda c, u: '%u:%u' % (c, u),
            c_u_diff = lambda c, u: c > u + 1,
            **kwargs
        ):
        """
        Represents a set of distinct substituent parts in a molecule or
        a homolog series spanning accross a range of aliphatic chain
        length and unsaturated bonds.
        
        :param list cores: List of core variations. Same kind of definitions
            are possible like at `formula.Formula`: formula as `str`, `dict` of
            atoms or exact mass as `float`.
        :param tuple c: Tuple of 2 integers: range of chain lengths.
        :param tuple u: Tuple of 2 integers: range of unsaturations.
        :param list counts: Dictionary with extra atom counts.
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
        :param list charges: List of integers: charges for each core
            variation.
        :param list isotopes: List of integers: extra neutrons for each core
            variation.
        :param getname: Method (callable) to create a name from the chain
            length and the unsaturation.
        :param c_u_diff: Method (callable) to decide if the chain length and
            unsaturation are compatible. E.g. if chain length is only C2, an
            unsaturation of 4 double bonds is not possible. By default, chain
            length must be greater by 2 than unsaturation, it means at acyl
            chains we avoid to assume double bond right next to the carboxyl
            group which is clearly impossible.
        """
        
        self.c_u_diff = c_u_diff
        self.getname  = getname
        self.cores    = cores if type(cores) is list else [cores]
        self.charges  = charges if type(charges) is list else [charges]
        self.isotopes = isotopes if type(isotopes) is list else [isotopes]
        self.names    = names if type(names) is list else [names]
        self.counts   = collections.defaultdict(lambda: 0)
        self.counts.update(counts)
        
        AbstractMetaboliteComponent.__init__(
            self,
            self.cores[0],
            charge = self.charges[0],
            isotope = self.isotopes[0],
            name = self.names[0],
            **kwargs
        )
        
        # range of possible lengths and unsats
        self.chlens = list(range(*c))
        self.unsats = list(range(*u))
        # current value of length and unsat
        self.c = self.chlens[0]
        self.u = self.unsats[0]
    
    def __iter__(self):
        
        for i in range(len(self.cores)):
            
            self.update_core(i)
            
            for c in self.chlens:
                
                self.c = c
                
                for u in self.unsats:
                    
                    self.u = u
                    
                    if not self.c_u_diff(self.c, self.u):
                        
                        continue
                    
                    # implicit hydrogens
                    h = c * 2 + 1 - 2 * u
                    new_counts = self.counts.copy()
                    new_counts['C'] += c
                    new_counts['H'] += h
                    
                    new =  self + formula.Formula(**new_counts)
                    new.name = self.getname(self.c, self.u)
                    
                    yield new
    
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
        
        self.charge = self.charges[i]
        self.isotope = self.isotopes[i]
        self.name = self.names[i]
        
        self.reset_atoms()
        self.calc_mass()
