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
import collections

import lipyd.metabolite as metabolite
import lipyd.substituent as substituent
import lipyd.formula as formula
import lipyd.common as common
import lipyd.lipproc as lipproc


fattyacids = [
    'FattyAcid'
]

misc = [
    'VitaminA'
]

# will be populated by the factory
glycerolipids = []
glycerophospholipids = []
sphingolipids = []

#
# Glycerolipids
#

class AbstractGlycerol(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            sn1 = 'H',
            sn2 = 'H',
            sn3 = 'H',
            charge = None,
            pcharge = 0,
            ncharge = 0,
            hg = None,
            name = 'Glycerol',
            typ  = 'G',
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = (
            charge
                if charge is not None else
            self.pcharge - self.ncharge
        )
        self.typ = typ
        hg = hg or lipproc.Headgroup(main = name)
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'C3H5O3',
            subs = [
                self.get_substituent(sn1),
                self.get_substituent(sn2),
                self.get_substituent(sn3)
            ],
            charge = self.netcharge,
            name = name,
            hg = hg,
            **kwargs
        )


#
# Glycerophospholipids
#

class AbstractGlycerolipid(AbstractGlycerol):
    
    def __init__(
            self,
            headgroup = 'H',
            phospho = True,
            lyso = False,
            ether = False,
            fa_args = None,
            name = 'GL',
            getname = None,
            hg = None,
            typ  = None,
            sn3_fa = False,
            lyso_sn1_fa = True,
            sn2_fa_args = None,
            sn3_fa_args = None,
            sn1_ether = False,
            sn2_ether = False,
            sn3_ether = False,
            **kwargs
        ):
        """
        Represents a generic glycerolipid.
        
        Args
        ----
        :param str headgroup:
            Formula of the moiety attached to the phosphate.
        :param bool phospho:
            Is this a glycerophospholipid? If `True` implicitely
            adds `PO3H` to sn3 position hence only the remaining
            part of the headgroup needs to be defined, e.g. choline.
        :param bool lyso:
            Whether it is a lyso form or not.
        :param bool ether:
            Whether it is an ether i.e. having fatty alcohol
            ether on one or both of the sn1 and sn2 positions.
        :param dict fa_args:
            Arguments for the `substituent.FattyAcyl()`.
        :param str name:
            Name stem of the lipid class.
        :param str typ:
            Name of the lipid family, here should be GL or GPL,
            which is set automatically depending on `phospho`.
        :param bool sn3_fa:
            Is there a 3rd fatty acyl or alkyl moiety at the sn3 position?
        :param bool lyso_sn1_fa:
            In case of lyso form is the alkyl ester/ether
            in sn1 position?
        :param dict sn2_fa_args:
            If the sn2 fatty acyl or alkyl has
            different parameters; if `None` it defaults to `fa_args`.
        :param dict sn3_fa_args:
            If the sn3 fatty acyl or alkyl has
            different parameters; if `None` it defaults to `fa_args`.
        :param bool sn1_ether:
            If `ether` is `True`, is there ether in sn1
            position? If `ether` is `True` but both this and `sn2_ether`
            are `False`, sn1 ether position assumed.
        :param bool sn2_ether:
            If `ether` is `True`, is there ether in sn2
            position?
        :param bool sn3_ether:
            If `ether` is `True`, is there ether in sn3
            position? For example to hace trialkyl-glycerol
            set `sn1_ether`, `sn2_ether` and `sn3_ether` all
            to `True`.
        :param **kwargs:
            Passed to `AbstractGlycerol` and finally to
            `metabolite.AbstractMetabolite`.
        """
        
        def get_cls(fa, ether):
            """
            Returns `substituent.FattyAcyl` or `substituent.FattyAlkoxy`
            class for sn1 and sn2 positions and or `None` if the hydroxyl
            group is free.
            """
            
            return (
                None
                if self.lyso and not fa else
                substituent.FattyAlkoxy
                if ether else
                substituent.FattyAcyl
            )
        
        self.fa_args = fa_args or {}
        self.lyso = lyso
        self.ether = ether
        # ether is default on sn1
        self.sn2_ether = sn2_ether
        self.sn3_ether = sn3_ether
        # both sn1 and sn2 can be ether if set explicitely
        # otherwise sn1 is ether sn2 is ester
        self.sn1_ether = sn1_ether or (ether and not sn2_ether)
        # in lyso forms the alkyl ester/ether is on the sn1 by default
        self.lyso_sn1_fa = lyso_sn1_fa
        
        # by default all fatty acyls have the same parameters
        # except if arguments given explicitely
        self.sn1_fa_args = self.fa_args
        self.sn2_fa_args = sn2_fa_args or self.fa_args
        self.sn3_fa_args = sn3_fa_args or self.fa_args
        self.sn1_cls = get_cls(lyso_sn1_fa, self.sn1_ether)
        self.sn2_cls = get_cls(not lyso_sn1_fa, self.sn2_ether)
        self.sn3_cls = None if not sn3_fa else get_cls(True, self.sn3_ether)
        
        def _getname(parent, subs):
            
            return (
                lipproc.summary_str(
                    self.hg,
                    lipproc.sum_chains(
                        tuple(
                            s.attrs.chain
                            for s in subs
                            if hasattr(s.attrs, 'chain')
                        )
                    )
                )
                if self.sum_only else
                lipproc.full_str(
                    self.hg,
                    tuple(
                        s.attrs.chain
                        for s in subs
                        if hasattr(s.attrs, 'chain')
                    )
                )
            )
        
        AbstractGlycerol.__init__(
            self,
            sn1  = (
                'H'
                    if self.sn1_cls is None else
                self.sn1_cls(**self.sn1_fa_args)
            ),
            sn2  = (
                'H'
                    if self.sn2_cls is None else
                self.sn2_cls(**self.sn2_fa_args)
            ),
            sn3  = (
                ('PO3H%s' % headgroup)
                    if phospho else
                self.sn3_cls(**self.sn3_fa_args)
                    if self.sn3_cls is not None else
                headgroup
            ),
            name = name,
            typ  = typ or ('GPL' if phospho else 'GL'),
            getname = getname or _getname,
            hg = hg,
            **kwargs
        )


# just an exemple (or reference) how to create a specific
# glycerolipid on top of the abstract one
class Phosphatidylethanolamine(AbstractGlycerolipid):
    
    def __init__(self, **kwargs):
        
        AbstractGPL.__init__(
            self,
            headgroup = 'C2H4NH2',
            name = 'PE',
            typ  = 'PE',
            hg = lipproc.Headgroup('PE'),
            pcharge = 1,
            ncharge = 1,
            **kwargs
        )


class GlycerolipidFactory(object):
    
    def __init__(
            self,
            fa1_args = None,
            fa2_args = None,
            double_ester = True,
            ether_ester = True,
            lyso_ester = True,
            lyso_ether = True,
            **kwargs
        ):
        """
        Populates glycerolipid classes.
        """
        
        fa1_args = fa1_args or {'c': (4, 24), 'u': (0, 9)}
        fa2_args = fa2_args or {'c': (4, 24), 'u': (0, 9)}
        
        l_classes = [
            ('C2H3NH3', 'PE', True),
            ('C2H3NC3H9', 'PC', True),
            ('C2H3NH2COOH', 'PS', True),
            ('C3H5O2H2', 'PG', True),
            ('C3H5O2H2', 'BMP', True),
            ('H', 'PA', True),
            ('C3O2H8PO3', 'PGP', True),
            ('C6O5H11', 'PI', True),
            ('C6O5H12PO3', 'PIP', True),
            ('C6O5H13P2O6', 'PIP2', True),
            ('C6O5H14P3O9', 'PIP3', True),
            ((True,), 'MAG', False),
            ((False,), 'MAG', False),
            ((True, True), 'DAG', False),
            ((False, True), 'DAG', False),
            ((True, True, True), 'TAG', False),
            ((False, True, True), 'TAG', False),
            ('C7H14NO2', 'DGTS', False),
            ('C7H14NO3', 'DGCC', False),
            ('C7H14NO2', 'DGTA', False),
            ('C6O7SH11', 'SQDG', False),
            ('C6O5H11', 'MGDG', False),
            ('C12O10H21', 'DGDG', False),
        ]
        l_ether = [False, True]
        l_lyso  = [False, True]
        
        docs = {
            'DGDG':
                """
                Example:
                    https://pubchem.ncbi.nlm.nih.gov/compound/52922098
                """,
            'PE':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000096355/
                    
                    [(m.name, m.mass) for m in
                        lipid.PE(
                            fa_args = {'c': 26, 'u': 0},
                            sn2_fa_args = {'c': 22, 'u': 6}
                        )
                    ]
                    
                    exact mass = 903.6717059967399
                """,
            'EtherPE':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000029613/
                    
                    [(m.name, m.mass) for m in
                        lipid.EtherPE(
                            fa_args = {'c': 16, 'u': 0},
                            sn2_fa_args = {'c': 28, 'u': 0}
                        )
                    ]
                    
                    exact mass = 845.7237415701001
                """,
            'LysoPE':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000055415/
                    
                    [(m.name, m.mass) for m in
                        lipid.EtherPE(
                            fa_args = {'c': 18, 'u': 1}
                        )
                    ]
                    
                    exact mass = 479.30118982806
                """,
            'MonoalkylGlycerol':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000308014/
                    
                    [(m.name, m.mass) for m in
                        lipid.MonoalkylGlycerol(
                            fa_args = {'c': 18, 'u': 2}
                        )
                    ]
                    
                    exact mass = 340.297745151
                """,
            'MonoalkylMonoacylGlycerol':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000321943/
                    
                    [(m.name, m.mass) for m in
                        lipid.MonoalkylMonoacylGlycerol(
                            fa_args = {'c': 20, 'u': 1},
                            sn2_fa_args = {'c': 24, 'u': 6}
                        )
                    ]
                    
                    exact mass = 708.6056610616
                """,
            'MonoalkylDiacylGlycerol':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000198889/
                    
                    [(m.name, m.mass) for m in
                        lipid.MonoalkylDiacylGlycerol(
                            fa_args = {'c': 18, 'u': 0},
                            sn2_fa_args = {'c': 17, 'u': 0},
                            sn3_fa_args = {'c': 22, 'u': 0}
                        )
                    ]
                    
                    exact mass = 918.89792690768
                """,
            'TriacylGlycerol':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000237541/
                    
                    [(m.name, m.mass) for m in
                        lipid.TriacylGlycerol(
                            fa_args = {'c': 18, 'u': 1},
                            sn2_fa_args = {'c': 26, 'u': 0},
                            sn3_fa_args = {'c': 16, 'u': 1}
                        )
                    ]
                    
                    exact mass = 970.89284152788
                """,
            'DiacylGlycerol':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000124274/
                    
                    [(m.name, m.mass) for m in
                        lipid.DiacylGlycerol(
                            fa_args = {'c': 18, 'u': 4},
                            sn2_fa_args = {'c': 24, 'u': 5}
                        )
                    ]
                    
                    exact mass = 690.5223253592001
                """,
            'PI':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000016800/
                    
                    [(m.name, m.mass) for m in
                        lipid.PI(
                            fa_args = {'c': 38, 'u': 5},
                            sn2_fa_args = {'c': 18, 'u': 3}
                        )
                    ]
                    
                    exact mass = 1130.7762306419602
                """,
            'PGP':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000412170/
                    
                    [(m.name, m.mass) for m in
                        lipid.PGP(
                            fa_args = {'c': 16, 'u': 2},
                            sn2_fa_args = {'c': 26, 'u': 5}
                        )
                    ]
                    
                    exact mass = 928.5230667049199
                """,
            'PC':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000011977/
                    
                    [(m.name, m.mass) for m in
                        lipid.PC(
                            fa_args = {'c': 18, 'u': 3},
                            sn2_fa_args = {'c': 34, 'u': 4}
                        )
                    ]
                    
                    exact mass = 999.76560638386
                """,
            'EtherPA':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000085779/
                    
                    [(m.name, m.mass) for m in
                        lipid.EtherPA(
                            fa_args = {'c': 15, 'u': 0},
                            sn2_fa_args = {'c': 30, 'u': 4}
                        )
                    ]
                    
                    exact mass = 808.6345922110401
                """,
            'PS':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000002856/
                    
                    [(m.name, m.mass) for m in
                        lipid.PS(
                            fa_args = {'c': 34, 'u': 6},
                            sn2_fa_args = {'c': 28, 'u': 5}
                        )
                    ]
                    
                    exact mass = 1133.80238581782
                """,
            'PG':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000072349/
                    
                    [(m.name, m.mass) for m in
                        lipid.PG(
                            fa_args = {'c': 20, 'u': 4},
                            sn2_fa_args = {'c': 34, 'u': 5}
                        )
                    ]
                    
                    exact mass = 1012.71323645876
                """,
            'PIP':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000495110/
                    
                    [(m.name, m.mass) for m in
                        lipid.PIP(
                            fa_args = {'c': 32, 'u': 5},
                            sn2_fa_args = {'c': 13, 'u': 0}
                        )
                    ]
                    
                    exact mass = 1062.61736101716
                """,
            'PIP2':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000491707/
                    
                    [(m.name, m.mass) for m in
                        lipid.PIP2(
                            fa_args = {'c': 16, 'u': 1},
                            sn2_fa_args = {'c': 38, 'u': 5}
                        )
                    ]
                    
                    exact mass = 1266.70889242468
                """,
            'PIP3':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000492994/
                    
                    [(m.name, m.mass) for m in
                        lipid.PIP3(
                            fa_args = {'c': 18, 'u': 1},
                            sn2_fa_args = {'c': 36, 'u': 5}
                        )
                    ]
                    
                    exact mass = 1346.6790633754044
                """,
        }
        
        mod = sys.modules[__name__]
        
        for (headgroup, name, phospho), ether, lyso in itertools.product(
                l_classes, l_ether, l_lyso
            ):
            
            if not lyso and not ether and not double_ester:
                continue
            if ether and not lyso and not ether_ester:
                continue
            if lyso and ether and not lyso_ether:
                continue
            if lyso and not ether and not lyso_ester:
                continue
            if not phospho and (ether or lyso):
                continue
            
            if not phospho:
                
                if type(headgroup) is tuple:
                    
                    # neutral glycerolipids
                    
                    faa = headgroup
                    sn1_ether = not faa[0]
                    sn2_ether = not faa[1] if len(faa) > 1 else False
                    sn3_ether = not faa[2] if len(faa) > 2 else False
                    lyso = len(faa) == 1
                    sn3_fa = len(faa) == 3
                    headgroup = 'H'
                    child = self.gl_class_name(faa)
                    
                else:
                    
                    sn1_ether = sn2_ether = sn3_ether = sn3_fa = False
                    child = name
                
            else:
                sn1_ether = sn2_ether = sn3_ether = sn3_fa = False
                child = self.gpl_class_name(name, lyso, ether)
            
            exec(
                (
                    'def __init__(self, **kwargs):\n'
                    '    \n'
                    '    AbstractGlycerolipid.__init__(\n'
                    '        self,\n'
                    '        headgroup = \'%s\',\n'
                    '        name = \'%s\',\n'
                    '        phospho = %s,\n'
                    '        lyso = %s,\n'
                    '        ether = %s,\n'
                    '        sn1_ether = %s,\n'
                    '        sn2_ether = %s,\n'
                    '        sn3_ether = %s,\n'
                    '        sn3_fa = %s,\n'
                    '        hg = lipproc.Headgroup(\n'
                    '            main = \'%s\',\n'
                    '            sub = (%s)\n'
                    '        ),\n'
                    '        **kwargs\n'
                    '        )\n'
                ) % (
                    headgroup,
                    '%s%s' % ('Lyso' if phospho and lyso else '', name),
                    str(phospho),
                    str(lyso),
                    str(ether),
                    str(sn1_ether),
                    str(sn2_ether),
                    str(sn3_ether),
                    str(sn3_fa),
                    name,
                    '\'Lyso\',' if lyso else ''
                ),
                mod.__dict__,
                mod.__dict__
            )
            
            if child in docs:
                
                mod.__dict__['__init__'].__doc__ = docs[child]
            
            cls = type(
                child,
                (AbstractGlycerolipid, ),
                {'__init__': mod.__dict__['__init__']}
            )
            
            setattr(mod, child, cls)
            
            if phospho:
                
                glycerophospholipids.append(child)
                
            else:
                
                glycerolipids.append(child)
        
        delattr(mod, '__init__')
    
    @staticmethod
    def gpl_class_name(name, lyso, ether):
        
        return '%s%s%s' % (
            'Lyso' if lyso else '',
            'Ether' if ether else '',
            name
        )
    
    @staticmethod
    def gl_class_name(faa):
        
        cnt = collections.Counter(faa)
        
        return (
            '%s%sGlycerol' % (
                '%salkyl' % common.count_prefix[cnt[False]].capitalize()
                    if False in cnt else
                '',
                '%sacyl' % common.count_prefix[cnt[True]].capitalize()
                    if True in cnt else
                ''
            )
        )

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
            sph_args  = None,
            charge    = None,
            pcharge   = None,
            ncharge   = None,
            lcb_type  = 'd',
            name = 'Sphingolipid',
            subtype = (),
            hg = None,
            getname = None,
            **kwargs
        ):
        
        sph_args = sph_args or {}
        
        self.pcharge = pcharge or 0
        self.ncharge = ncharge or 0
        self.charge  = (
                charge
            if type(charge) is int else
                self.pcharge - self.ncharge
        )
        
        self.lcb_type     = lcb_type
        self.hg = hg or lipproc.Headgroup(main = name, sub = subtype)
        
        if not sph:
            
            sph = substituent.Sphingosine(lcb_type = lcb_type, **sph_args)
        
        def _getname(parent, subs):
            
            return (
                lipproc.summary_str(
                    self.hg,
                    lipproc.sum_chains(
                        tuple(
                            s.attrs.chain
                            for s in subs
                            if hasattr(s.attrs, 'chain')
                        )
                    )
                )
                if self.sum_only else
                lipproc.full_str(
                    self.hg,
                    tuple(
                        s.attrs.chain
                        for s in subs
                        if hasattr(s.attrs, 'chain')
                    )
                )
            )
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = '',
            subs = [
                sph,
                self.get_substituent(n),
                self.get_substituent(o)
            ],
            name = name,
            hg = hg,
            charge = self.charge,
            getname = getname or _getname,
            **kwargs
        )

#
# Ceramides and Sphingomyelins
#

class AbstractCeramide(AbstractSphingolipid):
    
    def __init__(
            self,
            core = '',
            sph_args = None,
            fa_args = None,
            o = 'H',
            fa_hydroxy = False,
            lcb_type = 'd',
            name = 'Cer',
            hg = None,
            getname = None,
            **kwargs
        ):
        
        sph_args = sph_args or {}
        fa_args  = fa_args  or {}
        
        if fa_hydroxy:
            fa = substituent.HydroxyFattyAcyl(**fa_args)
        else:
            fa = substituent.FattyAcyl(**fa_args)
        
        self.fa_hydroxy = fa_hydroxy
        
        hg = hg or lipproc.Headgroup(main = 'Cer')
        
        AbstractSphingolipid.__init__(
            self,
            o = o,
            n = fa,
            sph_args = sph_args,
            name = name,
            getname = getname,
            lcb_type = lcb_type,
            hg = hg,
            **kwargs
        )


class SphingolipidFactory(object):
    
    lcb_name = {
        'd':  'D',
        't':  'T',
        'DH': 'Dihydro',
        'k':  'Keto',
    }
    
    pre_subclass_name = {
        'M1':    'Methyl',
        'M2':    'Dimethyl',
        'M3':    'Trimethyl',
        'Hex':   'Hexosyl',
        'Hex2':  'Dihexosyl',
        'SHex':  'SulfoHexosyl',
        'SHex2': 'SulfoDihexosyl',
    }
    
    post_subclass_name = {
        '1P': 'Phosphate',
        'PE': 'PE',
        'PI': 'PI',
    }
    
    def __init__(self, fa_args_1o = None, **kwargs):
        
        fa_args_1o = fa_args_1o or {'c': (4, 24), 'u': (0, 9)}
        
        l_lcb = ['d', 't', 'DH', 'k']
        l_fa_hydroxy = [True, False]
        l_classes = [
            # 1-O substituent, headgroup main class,
            # headgroup subclass, N-acyl
            ('H', 'Cer', (), True),
            ('PO3H2', 'Cer', ('1P',), True),
            ('C12H21O10SO3', 'Cer', ('SHex2',), True),
            ('C6H11O5SO3', 'Cer', ('SHex',), True),
            ('C12H21O10', 'Cer', ('Hex2',), True),
            ('C6H11O5', 'Cer', ('Hex',), True),
            ('PO3C2H4NC3H9', 'SM', (), True),
            ('PO3C2H4NH3', 'Cer', ('PE',), True),
            (None, 'Cer', ('1OAcyl',), True),
            ('H', 'Sph', (), False),
            ('PO3H2', 'Sph', ('1P',), False),
            ('HCH2', 'Sph', ('M1',), False),
            ('HC2H4', 'Sph', ('M2',), False),
            ('HC3H5', 'Sph', ('M3',), False),
        ]
        
        docs = {
            'CeramideDPE':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000398516/
                    
                    [(m.name, m.mass) for m in
                        lipid.CeramideDPE(
                            sph_args = {'c': (16, 16), 'u': (0, 0)},
                            fa_args = {'c': (30, 30), 'u': (6, 6)}
                        )
                    ]
                    
                    exact mass = 818.63017553472
                """,
            'SulfoDihexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000396884/
                    
                    [(m.name, m.mass) for m in
                        lipid.SulfoDihexosylCeramideD(
                            sph_args = {'c': (18, 18), 'u': (1, 1)},
                            fa_args = {'c': (18, 18), 'u': (1, 1)}
                        )
                    ]
                    
                    exact mass = 967.59020697344
                """,
            'DihexosylCeramideD':
                """
                Example:
                    http://www.swisslipids.org/#/entity/SLM:000395342/
                    
                    [(m.name, m.mass) for m in
                        lipid.DihexosylCeramideD(
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
        
        for lcb, fa_hydroxy, (o, name, subtype, nacyl) in (
            itertools.product(l_lcb, l_fa_hydroxy, l_classes)
        ):
            
            if (
                (lcb == 'k' and (nacyl or subtype)) or
                (fa_hydroxy and not nacyl)
            ):
                
                continue
            
            parent, child = self.class_name(
                name, subtype, lcb, fa_hydroxy, o, nacyl
            )
            
            exec(
                (
                    'def __init__(self, %s**kwargs):\n'
                    '    \n%s'
                    '    %s.__init__(\n'
                    '        self,\n'
                    '        o = %s,\n'
                    '        name = \'%s\',\n'
                    '        hg = lipproc.Headgroup(\n'
                    '            main = \'%s\',\n'
                    '            sub = %s,\n'
                    '        ),\n'
                    '        lcb_type = \'%s\',\n'
                    '        %s'
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
                    name,
                    name,
                    str(subtype),
                    lcb,
                    'fa_hydroxy = True,\n' if fa_hydroxy else '',
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
    
    def class_name(self, name, subtype, lcb, hydroxyacyl, o, nacyl):
        
        lcbtype = self.lcb_name[lcb] if lcb in self.lcb_name else ''
        
        main = '%s%s' % (
            'Ceramide'
                if name == 'Cer' else
            'Sphingomyelin'
                if nacyl else
            'Sphingosine',
            lcbtype
        )
        
        parent = 'AbstractCeramide' if nacyl else 'AbstractSphingolipid'
        
        child = '%s%s%s%s%s' % (
                'Hydroxyacyl' if hydroxyacyl else '',
                ''.join(
                    self.pre_subclass_name[s]
                    for s in subtype
                    if s in self.pre_subclass_name
                ),
                main,
                ''.join(
                    self.post_subclass_name[s]
                    for s in subtype
                    if s in self.post_subclass_name
                ),
                '1OAcyl' if o is None else '',
            )
        
        return parent, child

#
# Fatty acids
#

class FattyAcid(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            c = None,
            u = None,
            fa_counts = None,
            getname = None,
            sub = (),
            **kwargs
        ):
        
        def _getname(parent, subs):
            
            return '%s(%u:%u)' % (parent.name, subs[0].c, subs[0].u)
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'OH',
            subs = [
                substituent.FattyAcyl(c = c, u = u, counts = fa_counts)
            ],
            name = 'FA',
            hg = lipproc.Headgroup(main = 'FA', sub = ()),
            getname = getname or _getname,
            **kwargs
        )

#
# Vitamins
#

class VitaminA(metabolite.AbstractMetabolite):
    
    def __init__(
            self,
            c = None,
            u = None,
            fa_counts = None,
            getname = None,
            sub = (),
            **kwargs
        ):
        
        def _getname(parent, subs):
            
            return 'VA'
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = 'C19H26',
            subs = (
                metabolite.AbstractSubstituent(
                    cores = ('COOH', 'CH2OH'),
                    c = (0,),
                    u = (0,),
                ),
            ),
            name = 'VA',
            hg = lipproc.Headgroup(main = 'VA', sub = ()),
            getname = getname or _getname,
            **kwargs
        )


# creating further Sphingolipid derived classes:
_factory = SphingolipidFactory()
del _factory

# creating further Glycerolipid derived classes:
_factory = GlycerolipidFactory()
del _factory
