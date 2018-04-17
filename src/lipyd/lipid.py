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

import sys
import itertools

import lipyd.metabolite as metabolite
import lipyd.substituent as substituent
import lipyd.formula as formula


sphingolipids = [
    'Sphingosine',
    'SphingosinePhosphate',
    'HydroxySphingosine',
    'MethylSphingosine',
    'DiMethylSphingosine',
    'TriMethylSphingosine',
    'CeramideD',
    'CeramideT',
    'DihydroCeramide',
    'HydroxyacylCeramideD',
    'HydroxyacylCeramideT',
    'HydroxyacylDihydroCeramide'
]

#
# Glycerolipids
#

class AbstractGlycerol(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            sn1 = 'H',
            sn2 = 'H',
            sn3 = 'H',
            pcharge = 0,
            ncharge = 0,
            name = 'Glycerol',
            typ  = 'G',
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = self.pcharge - self.ncharge
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'C3H5O3',
            subs = [
                self.get_substituent(sn1),
                self.get_substituent(sn2),
                self.get_substituent(sn3)
            ],
            charge = self.netcharge,
            **kwargs
        )


#
# Glycerophospholipids
#

class AbstractGPL(AbstractGlycerol):
    
    def __init__(
            self,
            headgroup = 'H',
            lyso = False,
            ether = False,
            fa_args = {},
            name = 'GPL',
            typ  = 'GPL',
            lyso_sn1_fa = True,
            sn2_fa_args = None,
            sn1_ether = False,
            sn2_ether = False,
            **kwargs
        ):
        """
        Represents a generic glycerophospholipid.
        
        :param str headgroup: Formula of the moiety attached to the phosphate.
        :param bool lyso: Whether it is a lyso form or not.
        :param bool ether: Whether it is an ether i.e. having fatty alcohol
            ether on one or both of the sn1 and sn2 positions.
        :param dict fa_args: Arguments for the `substituent.FattyAcyl()`.
        :param str name: Name stem of the lipid class.
        :param str typ: Name of the lipid family, here should be GPL.
        :param bool lyso_sn1_fa: In case of lyso form is the alkyl ester/ether
            in sn1 position?
        :param dict sn2_fa_args: If the sn2 fatty acyl or alcohol has
            different parameters; if `None` it defaults to `fa_args`.
        :param bool sn1_ether: If `ether` is `True`, is there ether in sn1
            position? If `ether` is `True` but both this and `sn2_ether`
            are `False`, sn1 ether position assumed.
        :param bool sn2_ether: If `ether` is `True`, is there ether in sn2
            position?
        :param **kwargs: Passed to `AbstractGlycerol` and finally to
            `metabolite.AbstractMetabolite`.
        """
        
        def get_cls(lyso, fa, ether):
            """
            Returns `substituent.FattyAcyl` or `substituent.FattyAlkoxy`
            class for sn1 and sn2 positions and or `None` if the hydroxyl
            group is free.
            """
            
            return (
                None
                if lyso and not fa else
                substituent.FattyAlkoxy
                if ether else
                substituent.FattyAcyl
            )
        
        self.lyso = lyso
        self.ether = ether
        self.sn2_ether = sn2_ether
        self.sn1_ether = sn1_ether or (ether and not sn2_ether)
        self.lyso_sn1_fa = lyso_sn1_fa
        
        self.sn1_fa_args = fa_args
        self.sn2_fa_args = sn2_fa_args or fa_args
        self.sn1_cls = get_cls(lyso, lyso_sn1_fa, self.sn1_ether)
        self.sn1_cls = get_cls(lyso, not lyso_sn1_fa, self.sn2_ether)
        
        AbstractGlycerol.__init__(
            self,
            sn1  = 'H' if sn1_cls is None else sn1_cls(**sn1_fa_args),
            sn2  = 'H' if sn2_cls is None else sn2_cls(**sn2_fa_args),
            sn3  = 'PO3H%s' % headgroup,
            name = name,
            typ  = typ,
            **kwargs
        )


class Phosphatidylethanolamine(AbstractGPL):
    
    def __init__(self, **kwargs):
        
        AbstractGPL.__init__(
            self,
            headgroup = 'C2H4NH2',
            name = 'PE',
            typ  = 'PE',
            **kwargs
        )


class GPLFactory(object):
    
    def __init__(
            self,
            fa1_args = {'c': (4, 24), 'u': (0, 9)},
            fa2_args = {'c': (4, 24), 'u': (0, 9)},
            double_ester = True,
            ether_ester = True,
            lyso = True,
            lyso_ether = True,
            **kwargs
        ):
        
        l_classes = [
            ('C2H4NH2', 'PE'),
            ('C2H4NC3H9', 'PC'),
            ('C2H4NH2COOH', 'PS'),
            ('C3O2H5', 'PG'),
            ('H', 'PA'),
            ('C3O2H5PO3', 'PGP'),
            ('C6O5H11', 'PI'),
            ('C6O5H11PO3', 'PIP'),
            ('C6O5H11P3O6', 'PIP2'),
            ('C6O5H11P3O9', 'PIP3')
        ]
        l_dihydro = [True, False]
        
        docs = {
            'CeramideDPhosphoethanolamine':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000398516/
                    
                    [(m.name, m.mass) for m in
                        lipid.CeramideDPhosphoethanolamine(
                            sph_args = {'c': (16, 16), 'u': (0, 0)},
                            fa_args = {'c': (30, 30), 'u': (6, 6)}
                        )
                    ]
                    
                    exact mass = 818.63017553472
                """
        }
        
        mod = sys.modules[__name__]
        
        for t, fa_hydroxy, dihydro, (o, name) in itertools.product(
                l_t, l_fa_hydroxy, l_dihydro, l_classes
            ):
            
            if (t and dihydro):
                
                continue
            
            parent, child = self.class_name(name, t, dihydro, fa_hydroxy, o)
            
            exec(
                (
                    'def __init__(self, %s**kwargs):\n'
                    '    \n%s'
                    '    %s.__init__(\n'
                    '        self,\n'
                    '        o = %s,\n'
                    '        name = \'%s\',\n'
                    '        **kwargs\n'
                    '        )\n'
                ) % (
                    # `fa_args_1o` is an argument for
                    # 1-O-acyl ceramides
                    (
                        '\nfa_args_1o = %s,\n' % fa_args_1o.__str__()
                    )
                    if o is None
                    else '',
                    # the 1O substituent is a fatty acyl
                    # if o is None
                    (
                        '\n    fa1o = substituent.FattyAcyl('
                        '**fa_args_1o)\n'
                    )
                    if o is None
                    else '',
                    parent,
                    'fa1o' if o is None else '\'%s\'' % o,
                    name
                ),
                mod.__dict__,
                mod.__dict__
            )
            
            if child in docs:
                
                mod.__dict__['__init__'].__doc__ = docs[child]
            
            cls = type(
                child,
                (getattr(mod, parent), ),
                {'__init__': mod.__dict__['__init__']}
            )
            
            setattr(mod, child, cls)
            
            sphingolipids.append(child)
        
        delattr(mod, '__init__')
    
    def class_name(self, name, t, dihydro, hydroxyacyl, o):
        
        dt = 'T' if t else 'D' if not dihydro else ''
        
        maintype = '%s%s' % (
            'Ceramide' if 'Cer' in name else 'Sphingomyelin',
            dt
        )
        
        parmaintype = 'Ceramide%s' % dt
        
        subtype = '%s%s' % (
            'Hydroxyacyl' if hydroxyacyl else '',
            'Dihydro' if dihydro else ''
        )
        
        parent = '%s%s' % (subtype, parmaintype)
        
        child = '%s%s%s%s%s' % (
                subtype,
                'Sulfo' if 'SHex' in name else '',
                'DiHexosyl'
                    if 'Hex2' in name
                    else 'Hexosyl' if 'Hex' in name
                    else '',
                maintype,
                'Phosphoethanolamine'
                    if 'CerPE' in name
                    else 'Phosphate' if 'CerP' in name
                    else '1OAcyl' if o is None
                    else ''
            )
        
        return parent, child

#
# Shpingolipids
#

class AbstractSphingolipid(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            core = '',
            o = 'H',
            n = 'H',
            sph = None,
            sph_args = {},
            pcharge = 0,
            ncharge = 0,
            name = 'Sphingolipid',
            typ  = 'SL',
            getname = None,
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = self.pcharge - self.ncharge
        
        if not sph:
            
            sph = substituent.Sphingosine(**sph_args)
        
        def _getname(parent, subs):
            
            sep = '-' if self.sum_only else '/'
            chains = [
                s for s in subs
                if parent.has_variable_aliphatic_chain(s)
            ]
            
            return (
                '%s(%s%s)' % (
                    parent.name,
                    sep.join(
                        n for n in (
                            '%s%s' % (
                                s.get_prefix(),
                                '' if self.sum_only else s.name
                            )
                            for s in chains
                        ) if n
                    ),
                    (
                        '%u:%u' % (
                            sum(s.c for s in chains),
                            sum(s.u for s in chains)
                        ) if self.sum_only else ''
                    )
                )
            )
        
        metabolite.AbstractMetabolite.__init__(
            core = '',
            self,
            subs = [
                sph,
                self.get_substituent(n),
                self.get_substituent(o)
            ],
            name = name,
            charge = self.netcharge,
            getname = getname or _getname,
            **kwargs
        )

#
# Sphingoid bases
#

class Sphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        AbstractSphingolipid.__init__(
            self,
            sph_args = sph_args,
            name = 'Sph',
            **kwargs
        )


class SphingosinePhosphate(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        """
        Example:
            http://www.swisslipids.org/#/entity/SLM:000000438/
            
            [(m.name, m.mass) for m in
                lipid.SphingosinePhosphate(
                    sph_args = {'c': (18, 18), 'u': (1, 1)}
                )
            ]
            
            exact mass = 379.24876032958
        
        """
        
        AbstractSphingolipid.__init__(
            self,
            sph_args = sph_args,
            name = 'SphP',
            o = 'PO3H2',
            **kwargs
        )


class HydroxySphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        sph = substituent.HydroxySphingosine(**sph_args)
        
        AbstractSphingolipid.__init__(
            self,
            sph = sph,
            name = 'SphOH',
            **kwargs
        )



class MethylSphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        AbstractSphingolipid.__init__(
            self,
            n = 'CH3',
            sph_args = sph_args,
            name = 'SphM',
            **kwargs
        )


class DiMethylSphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        AbstractSphingolipid.__init__(
            self,
            n = 'C2H6',
            sph_args = sph_args,
            name = 'SphM2',
            **kwargs
        )


class TriMethylSphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        AbstractSphingolipid.__init__(
            self,
            n = 'C3H9',
            sph_args = sph_args,
            name = 'SphM3',
            **kwargs
        )

#
# Ceramides and Sphingomyelins
#

class AbstractCeramide(AbstractSphingolipid):
    
    def __init__(
            self,
            core = '',
            sph_args = {},
            fa_args = {},
            o = 'H',
            fa_hydroxy = False,
            t = False,
            dihydro = False,
            name = 'Cer',
            getname = None,
            **kwargs
        ):
        
        if t:
            
            sph = substituent.HydroxySphingosine(**sph_args)
        elif dihydro:
            sph = substituent.DihydroSphingosine(**sph_args)
        else:
            sph = substituent.Sphingosine(**sph_args)
        
        if fa_hydroxy:
            fa = substituent.HydroxyFattyAcyl(**fa_args)
        else:
            fa = substituent.FattyAcyl(**fa_args)
        
        self.t = t
        self.dihydro = dihydro
        self.fa_hydroxy = fa_hydroxy
        
        #def getname(parent, subs):
            
            #return (
                #'%s(%s%s/%s%s)' % (
                    #parent.name,
                    #subs[0].get_prefix(),
                    #subs[0].name,
                    #subs[1].get_prefix(),
                    #subs[1].name
                #)
            #)
        
        AbstractSphingolipid.__init__(
            self,
            o = o,
            n = fa,
            sph = sph,
            name = name,
            getname = getname,
            **kwargs
        )


class CeramideD(AbstractCeramide):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        """
        Example:
            http://www.swisslipids.org/#/entity/SLM:000397236/
            
            [(m.name, m.mass) for m in
                lipid.CeramideD(
                    sph_args = {'c': (14, 14), 'u': (1, 1)},
                    fa_args = {'c': (20, 20), 'u': (1, 1)}
                )
            ]
            
            exact mass = 535.4964449617
        """
        
        AbstractCeramide.__init__(
            self,
            o = o,
            **kwargs
        )


class CeramideT(AbstractCeramide):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        """
        Example:
            http://www.swisslipids.org/#/entity/SLM:000395636/
            
            [(m.name, m.mass) for m in
                lipid.CeramideT(
                    sph_args = {'c': (18, 18), 'u': (0, 0)},
                    fa_args = {'c': (10, 10), 'u': (0, 0)}
                )
            ]
            
            exact mass = 471.42875932382003
        """
        
        AbstractCeramide.__init__(
            self,
            o = o,
            t = True,
            **kwargs
        )


class DihydroCeramide(AbstractCeramide):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        """
        This is isobaric with dCer.
        """
        
        AbstractCeramide.__init__(
            self,
            dihydro = True,
            o = o,
            **kwargs
        )


class HydroxyacylCeramideD(CeramideD):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        """
        Example:
            http://www.swisslipids.org/#/entity/SLM:000397236/
            
            [(m.name, m.mass) for m in
                lipid.HydroxyacylCeramideD(
                    sph_args = {'c': (14, 14), 'u': (1, 1)},
                    fa_args = {'c': (20, 20), 'u': (1, 1)}
                )
            ]
            
            exact mass = 551.4913595819
            
            551.4913595819 - formula.Formula('O').mass = 535.4964449617
            
            I.e. this species is always one oxygen heavier than the
            corresponding dCer.
        """
        
        CeramideD.__init__(
            self,
            fa_hydroxy = True,
            o = o,
            **kwargs
        )


class HydroxyacylCeramideT(CeramideT):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            fa_hydroxy = True,
            o = o,
            **kwargs
        )


class HydroxyacylDihydroCeramide(CeramideD):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        
        CeramideD.__init__(
            self,
            dihydro = True,
            fa_hydroxy = True,
            o = o,
            **kwargs
        )



class CeramideFactory(object):
    
    def __init__(self, fa_args_1o = {'c': (4, 24), 'u': (0, 9)}, **kwargs):
        
        l_t = [True, False]
        l_fa_hydroxy = [True, False]
        l_classes = [
            ('PO3H2', 'CerP'),
            ('C12H21O10SO3', 'SHex2Cer'),
            ('C6H11O5SO3', 'SHexCer'),
            ('C12H21O10', 'Hex2Cer'),
            ('C6H11O5', 'HexCer'),
            ('PO3C2H4NC3H9', 'SM'),
            ('PO3C2H4NH3', 'CerPE'),
            (None, 'CerA')
            
        ]
        l_dihydro = [True, False]
        
        docs = {
            'CeramideDPhosphoethanolamine':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000398516/
                    
                    [(m.name, m.mass) for m in
                        lipid.CeramideDPhosphoethanolamine(
                            sph_args = {'c': (16, 16), 'u': (0, 0)},
                            fa_args = {'c': (30, 30), 'u': (6, 6)}
                        )
                    ]
                    
                    exact mass = 818.63017553472
                """,
            'SulfoDiHexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000396884/
                    
                    [(m.name, m.mass) for m in
                        lipid.SulfoDiHexosylCeramideD(
                            sph_args = {'c': (18, 18), 'u': (1, 1)},
                            fa_args = {'c': (18, 18), 'u': (1, 1)}
                        )
                    ]
                    
                    exact mass = 967.59020697344
                """,
            'DiHexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000395342/
                    
                    [(m.name, m.mass) for m in
                        lipid.DiHexosylCeramideD(
                            sph_args = {'c': (18, 18), 'u': (1, 1)},
                            fa_args = {'c': (26, 26), 'u': (0, 0)}
                        )
                    ]
                    
                    exact mass = 1001.77424251862
                """,
            'SulfoHexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000396804/
                    
                    [(m.name, m.mass) for m in
                        lipid.SulfoHexosylCeramideD(
                            sph_args = {'c': (18, 18), 'u': (1, 1)},
                            fa_args = {'c': (18, 18), 'u': (2, 2)}
                        )
                    ]
                    
                    exact mass = 803.5217334853199
                """,
            'HexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000395423/
                    
                    [(m.name, m.mass) for m in
                        lipid.HexosylCeramideD(
                            sph_args = {'c': (18, 18), 'u': (1, 1)},
                            fa_args = {'c': (16, 16), 'u': (0, 0)}
                        )
                    ]
                    
                    exact mass = 699.56491844982
                """,
            'SphingomyelinT':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000485623/
                    
                    [(m.name, m.mass) for m in
                        lipid.SphingomyelinT(
                            sph_args = {'c': (17, 17), 'u': (0, 0)},
                            fa_args = {'c': (32, 32), 'u': (5, 5)}
                        )
                    ]
                    
                    exact mass = 920.73464060656
                """,
            'SphingomyelinD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000397988/
                    
                    [(m.name, m.mass) for m in
                        lipid.SphingomyelinD(
                            sph_args = {'c': (16, 16), 'u': (1, 1)},
                            fa_args = {'c': (18, 18), 'u': (1, 1)}
                        )
                    ]
                    
                    exact mass = 700.5519252121201
                """,
            'HydroxyacylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000397236/
                    
                    [(m.name, m.mass) for m in
                        lipid.HydroxyacylCeramideD(
                            sph_args = {'c': (14, 14), 'u': (1, 1)},
                            fa_args = {'c': (20, 20), 'u': (1, 1)}
                        )
                    ]
                    
                    exact mass = 551.4913595819
                    
                    551.4913595819 - formula.Formula('O').mass = 535.4964449617
                    
                    I.e. this species is always one oxygen heavier than the
                    corresponding dCer.
                """,
            'DihydroCeramide':
                """
                This is isobaric with dCer.
                """,
            'CeramideT':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000395636/
                    
                    [(m.name, m.mass) for m in
                        lipid.CeramideT(
                            sph_args = {'c': (18, 18), 'u': (0, 0)},
                            fa_args = {'c': (10, 10), 'u': (0, 0)}
                        )
                    ]
                    
                    exact mass = 471.42875932382003
                """,
            'CeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000397236/
                    
                    [(m.name, m.mass) for m in
                        lipid.CeramideD(
                            sph_args = {'c': (14, 14), 'u': (1, 1)},
                            fa_args = {'c': (20, 20), 'u': (1, 1)}
                        )
                    ]
                    
                    exact mass = 535.4964449617
                """
        }
        
        mod = sys.modules[__name__]
        
        for t, fa_hydroxy, dihydro, (o, name) in itertools.product(
                l_t, l_fa_hydroxy, l_dihydro, l_classes
            ):
            
            if (t and dihydro):
                
                continue
            
            parent, child = self.class_name(name, t, dihydro, fa_hydroxy, o)
            
            exec(
                (
                    'def __init__(self, %s**kwargs):\n'
                    '    \n%s'
                    '    %s.__init__(\n'
                    '        self,\n'
                    '        o = %s,\n'
                    '        name = \'%s\',\n'
                    '        **kwargs\n'
                    '        )\n'
                ) % (
                    # `fa_args_1o` is an argument for
                    # 1-O-acyl ceramides
                    (
                        '\nfa_args_1o = %s,\n' % fa_args_1o.__str__()
                    )
                    if o is None
                    else '',
                    # the 1O substituent is a fatty acyl
                    # if o is None
                    (
                        '\n    fa1o = substituent.FattyAcyl('
                        '**fa_args_1o)\n'
                    )
                    if o is None
                    else '',
                    parent,
                    'fa1o' if o is None else '\'%s\'' % o,
                    name
                ),
                mod.__dict__,
                mod.__dict__
            )
            
            if child in docs:
                
                mod.__dict__['__init__'].__doc__ = docs[child]
            
            cls = type(
                child,
                (getattr(mod, parent), ),
                {'__init__': mod.__dict__['__init__']}
            )
            
            setattr(mod, child, cls)
            
            sphingolipids.append(child)
        
        delattr(mod, '__init__')
    
    def class_name(self, name, t, dihydro, hydroxyacyl, o):
        
        dt = 'T' if t else 'D' if not dihydro else ''
        
        maintype = '%s%s' % (
            'Ceramide' if 'Cer' in name else 'Sphingomyelin',
            dt
        )
        
        parmaintype = 'Ceramide%s' % dt
        
        subtype = '%s%s' % (
            'Hydroxyacyl' if hydroxyacyl else '',
            'Dihydro' if dihydro else ''
        )
        
        parent = '%s%s' % (subtype, parmaintype)
        
        child = '%s%s%s%s%s' % (
                subtype,
                'Sulfo' if 'SHex' in name else '',
                'DiHexosyl'
                    if 'Hex2' in name
                    else 'Hexosyl' if 'Hex' in name
                    else '',
                maintype,
                'Phosphoethanolamine'
                    if 'CerPE' in name
                    else 'Phosphate' if 'CerP' in name
                    else '1OAcyl' if o is None
                    else ''
            )
        
        return parent, child


# creating further Ceramide derived classes:
_factory = CeramideFactory()
del _factory
