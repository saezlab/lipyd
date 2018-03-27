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
            name = name,
            charge = self.netcharge,
            **kwargs
        )


class AbstractGPL(AbstractGlycerol):
    
    def __init__(
            self,
            headgroup = 'H',
            lyso = False,
            ether = False,
            fa_args = {},
            name = 'GPL',
            typ  = 'GPL',
            **kwargs
        ):
        
        AbstractGlycerol.__init__(
            self,
            sn1  = substituent.FattyAcyl(**fa_args),
            sn2  = 'H' if lyso else substituent.FattyAcyl(**fa_args),
            sn3  = 'PO3H%s' % headgroup,
            name = name,
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
            **kwargs
        ):
        
        self.pcharge = pcharge
        self.ncharge = ncharge
        self.netcharge = self.pcharge - self.ncharge
        
        if not sph:
            
            sph = substituent.Sphingosine(**sph_args)
        
        metabolite.AbstractMetabolite.__init__(
            self,
            core = '',
            subs = [
                sph,
                self.get_substituent(n),
                self.get_substituent(o)
            ],
            name = name,
            charge = self.netcharge,
            **kwargs
        )


class Sphingosine(AbstractSphingolipid):
    
    def __init__(self, sph_args = {}, **kwargs):
        
        AbstractSphingolipid.__init__(
            self,
            sph_args = sph_args,
            name = 'Sph',
            **kwargs
        )


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
        
        def getname(parent, subs):
            
            return (
                '%s%s(%s%s/%s%s)' % (
                    'DH'
                        if parent.dihydro or (
                            not parent.t and subs[0].u == 0
                        )
                        else '',
                    parent.name,
                    't' if parent.t else 'd',
                    subs[0].getname(),
                    subs[1].getname(),
                    '-2OH' if parent.fa_hydroxy else ''
                )
            )
        
        AbstractSphingolipid.__init__(
            self,
            o = o,
            n = fa,
            sph = sph,
            name = name,
            getname = getname,
            **kwargs
        )

# Ceramides and Sphingomyelins

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


class SphingomyelinD(CeramideD):
    
    def __init__(self, **kwargs):
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
        """
        
        CeramideD.__init__(self, o = 'PO3C2H4NC3H9', name = 'SM', **kwargs)


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


class SphingomyelinT(CeramideT):
    
    def __init__(self, **kwargs):
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
        """
        
        CeramideT.__init__(self, o = 'PO3C2H4NC3H9', name = 'SM', **kwargs)


class DihydroCeramide(AbstractCeramide):
    
    def __init__(
            self,
            o = 'H',
            **kwargs
        ):
        """
        This is perfectly isobaric with dCer.
        """
        
        AbstractCeramide.__init__(
            self,
            dihydro = True,
            **kwargs
        )


class DihydroSphingomyelin(DihydroCeramide):
    
    def __init__(self, **kwargs):
        
        DihydroCeramide.__init__(
            self,
            o = 'PO3C2H4NC3H9',
            name = 'SM',
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
            **kwargs
        )


class HydroxyacylSphingomyelinD(HydroxyacylCeramideD):
    
    def __init__(self, **kwargs):
        
        HydroxyacylCeramideD.__init__(
            self,
            o = 'PO3C2H4NC3H9',
            name = 'SM',
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
            **kwargs
        )


class HydroxyacylSphingomyelinT(HydroxyacylCeramideT):
    
    def __init__(self, **kwargs):
        
        HydroxyacylCeramideT.__init__(
            self,
            o = 'PO3C2H4NC3H9',
            name = 'SM',
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
            **kwargs
        )


class HydroxyacylDihydroSphingomyelinT(HydroxyacylDihydroCeramide):
    
    def __init__(self, **kwargs):
        
        HydroxyacylDihydroCeramide.__init__(
            self,
            o = 'PO3C2H4NC3H9',
            name = 'SM',
            **kwargs
        )


class CeramideDPhosphate(CeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideD.__init__(
            self,
            o = 'PO3H2',
            name = 'CerP',
            **kwargs
        )


class CeramideTPhosphate(CeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            o = 'PO3H2',
            name = 'CerP',
            **kwargs
        )

class HydroxyacylCeramideDPhosphate(HydroxyacylCeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HidroxyacylCeramideD.__init__(
            self,
            o = 'PO3H2',
            name = 'CerP',
            **kwargs
        )


class HydroxyacylCeramideTPhosphate(HydroxyacylCeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HidroxyacylCeramideT.__init__(
            self,
            o = 'PO3H2',
            name = 'CerP',
            **kwargs
        )

class DihydroCeramidePhosphate(DihydroCeramide):
    
    def __init__(
            self,
            **kwargs
        ):
        
        DihydroCeramide.__init__(
            self,
            o = 'PO3H2',
            name = 'CerP',
            **kwargs
        )

# Hexosyl Ceramides

class HexosylCeramideD(CeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
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
        """
        
        CeramideD.__init__(
            self,
            o = 'C6H11O5',
            name = 'HexCer',
            **kwargs
        )


class HexosylCeramideT(CeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            o = 'C6H11O5',
            name = 'HexCer',
            **kwargs
        )


class HydroxyacylHexosylCeramideD(HydroxyacylCeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideD.__init__(
            self,
            o = 'C6H11O5',
            name = 'HexCer',
            **kwargs
        )


class HydroxyacylHexosylCeramideT(HydroxyacylCeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideT.__init__(
            self,
            o = 'C6H11O5',
            name = 'HexCer',
            **kwargs
        )


class DihydroHexosylCeramide(DihydroCeramide):
    
    def __init__(
            self,
            **kwargs
        ):
        
        DihydroCeramide.__init__(
            self,
            o = 'C6H11O5',
            name = 'HexCer',
            **kwargs
        )

# Dihexosyl Ceramides

class DiHexosylCeramideD(CeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
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
        """
        
        CeramideD.__init__(
            self,
            o = 'C12H21O10',
            name = 'Hex2Cer',
            **kwargs
        )


class DiHexosylCeramideT(CeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            o = 'C12H21O10',
            name = 'Hex2Cer',
            **kwargs
        )


class HydroxyacylDiHexosylCeramideD(HydroxyacylCeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideD.__init__(
            self,
            o = 'C12H21O10',
            name = 'Hex2Cer',
            **kwargs
        )


class HydroxyacylDiHexosylCeramideT(HydroxyacylCeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideT.__init__(
            self,
            o = 'C12H21O10',
            name = 'Hex2Cer',
            **kwargs
        )


class DihydroDiHexosylCeramide(DihydroCeramide):
    
    def __init__(
            self,
            **kwargs
        ):
        
        DihydroCeramide.__init__(
            self,
            o = 'C12H21O10',
            name = 'Hex2Cer',
            **kwargs
        )


# Hexosyl Ceramides

class SulfoHexosylCeramideD(CeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
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
        """
        
        CeramideD.__init__(
            self,
            o = 'C6H11O5SO3',
            name = 'SHexCer',
            **kwargs
        )


class SulfoHexosylCeramideT(CeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            o = 'C6H11O5SO3',
            name = 'SHexCer',
            **kwargs
        )


class HydroxyacylSulfoHexosylCeramideD(HydroxyacylCeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideD.__init__(
            self,
            o = 'C6H11O5SO3',
            name = 'SHexCer',
            **kwargs
        )


class HydroxyacylSulfoHexosylCeramideT(HydroxyacylCeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideT.__init__(
            self,
            o = 'C6H11O5SO3',
            name = 'SHexCer',
            **kwargs
        )


class DihydroSulfoHexosylCeramide(DihydroCeramide):
    
    def __init__(
            self,
            **kwargs
        ):
        
        DihydroCeramide.__init__(
            self,
            o = 'C6H11O5SO3',
            name = 'SHexCer',
            **kwargs
        )

# Sulfodihexosyl Ceramides

class SulfoDiHexosylCeramideD(CeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
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
        """
        
        CeramideD.__init__(
            self,
            o = 'C12H21O10SO3',
            name = 'SHex2Cer',
            **kwargs
        )


class SulfoDiHexosylCeramideT(CeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        CeramideT.__init__(
            self,
            o = 'C12H21O10SO3',
            name = 'SHex2Cer',
            **kwargs
        )


class HydroxyacylSulfoDiHexosylCeramideD(HydroxyacylCeramideD):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideD.__init__(
            self,
            o = 'C12H21O10SO3',
            name = 'SHex2Cer',
            **kwargs
        )


class HydroxyacylSulfoDiHexosylCeramideT(HydroxyacylCeramideT):
    
    def __init__(
            self,
            **kwargs
        ):
        
        HydroxyacylCeramideT.__init__(
            self,
            o = 'C12H21O10SO3',
            name = 'SHex2Cer',
            **kwargs
        )


class DihydroSulfoDiHexosylCeramide(DihydroCeramide):
    
    def __init__(
            self,
            **kwargs
        ):
        
        DihydroCeramide.__init__(
            self,
            o = 'C12H21O10SO3',
            name = 'SHex2Cer',
            **kwargs
        )


class CeramideFactory(object):
    
    def __init__(self):
        
        l_t = [True, False]
        l_fa_hydroxy = [True, False]
        l_classes = [
            ('PO3H2', 'CerP')
        ]
        l_dihydro = [True, False]
        
        mod = sys.modules[__name__]
        
        for t, fa_hydroxy, dihydro, (o, name) in itertools.product(
            l_t, l_fa_hydroxy, l_dihydro, l_classes
        ):
            
            def init():
                
                pass
            
            cls = type()
            #setattr(mod)
