
# coding: utf-8

# # **lipyd** – A Python module for  
# # lipidomics LC MS/MS data analysis

# *Author:* Dénes Türei  
# *Contact:* turei.denes@gmail.com  
# *Git repositories (mirrored):*  
# - https://git.embl.de/turei/lipyd
# - https://bitbucket.org/deeenes/lipyd
# - https://github.com/saezlab/lipyd

# ## 1: Chemical calculator

# In[46]:


import imp
import itertools
import pprint
from lipyd import pprint_namedtuple
pprint.PrettyPrinter = pprint_namedtuple.PrettyPrinter
from lipyd import mass
from lipyd import formula
imp.reload(mass)
imp.reload(formula)


# The `mass` module knows exact masses of isotopes, isotopic abundances, weights, etc. The `MassBase` class is able to process chemical formula, calculate masses and do arithmetics.

# In[3]:


mass.get_mass('Na')


# In[6]:


mass.MassBase('C2H6') - mass.MassBase('H') + mass.MassBase('OH')


# Expression evaluation with `mass.expr`:

# In[8]:


mass.expr('C2H6 - H + OH')


# In[9]:


mass.expr('C6H12O6 - water')


# Make a deuterium:

# In[11]:


mass.expr('H + n')


# Hydrogensulphate ion:

# In[12]:


mass.expr('HSO4 + e')


# Additional attributes can be provided in keyword arguments to carry metadata.

# In[40]:


lactic_acid = formula.Formula('CH3CHOHCOOH', name = 'lactic acid')
lactic_acid.attrs.name


# A galactose:

# In[15]:


((2 * formula.Formula('C6H12O6')) - 'H2O').formula


# ## 2: Calculations with adducts

# In[4]:


from lipyd import mz


# This is an oleic acid. What is the mass of the [M-H]- adduct?

# In[17]:


formula.Formula('C18H34O2').remove_h()


# We have seen a mass and wondering what is the exact mass if it is an [M+NH4]+ adduct:

# In[22]:


mz.Mz(854.576482597791).remove_nh4()


# Calculate the [M+Li]+ adduct for the same molecule:

# In[24]:


mz.Mz(836.5426570444603).adduct(mass.MassBase('Li', charge = 1))


# ## 3: Metabolite model

# Metabolites consist of a core and optionally substituents. Substituents might be formulas or moieties with aliphatic chains.

# In[5]:


from lipyd import metabolite
from lipyd import substituent
imp.reload(metabolite)
imp.reload(substituent)


# Make all combinations of halogenated methanes:

# In[6]:


halo_methanes = metabolite.AbstractMetabolite(
    core = 'C',
    subs = [('H', 'F', 'Cl', 'Br', 'I')] * 4
)


# Check the first 3:

# In[7]:


[(m.formula, m.mass) for m in halo_methanes][:3]


# Do the same with all alcohols up to 1-8 carbon count with 0-2 unsaturated bonds:

# In[34]:


chain = metabolite.AbstractSubstituent(c = (1, 8), u = (0, 2))
alcohols = metabolite.AbstractMetabolite(subs = (chain, ('OH',)))
[(m.formula, m.mass) for m in alcohols][:3]


# Make some ceramides:

# In[9]:


# fatty acyls of length 16, 18 or 20 and one or no unsaturation:
fattyacyl = substituent.FattyAcyl(c = (16, 18, 20), u = (0, 1))
lcb = substituent.Sphingosine(c = 18, u = 1)
ceramides = metabolite.AbstractMetabolite(core = 'H', subs = (lcb, fattyacyl), name = 'Cer')
# name, formula, mass, [M+H]+, [M+NH4]+
[(cer.name, cer.formula, cer.mass, cer.add_nh4()) for cer in ceramides]


# ## 4: Lipid definitions

# In the `lipyd.lipid` module more than 150 lipid varieties are predefined.

# In[11]:


from lipyd import lipid


# In[16]:


d_ceramides = lipid.CeramideD(
    sph_args = {'c': 18, 'u': (0, 1)},
    fa_args  = {'c': (14, 24), 'u': (0, 4), 'even': True},
)


# In[30]:


[(cer.name, cer.mass) for cer in d_ceramides][:3]


# ## 5: External databases

# The `lipyd.moldb` module provides access to *SwissLipis* and *LipidMaps*. Automatically downloads and processes the databases which you can search by various names and identifiers, you can also access structures as OpenBabel objects, InChI and SMILE strings.

# In[39]:


from lipyd import moldb
from lipyd import name
from lipyd.lipproc import *


# ### 5.1: SwissLipids

# In[40]:


swl = moldb.SwissLipids(levels = {'species'})


# Get phosphatidylethanolamines as OpenBabel objects:

# In[48]:


swl.reload()
pe = swl.get_hg_obmol('PE')
pe0 = next(pe)
pe0.draw()


# In[49]:


[m.title for m in itertools.islice(pe, 0, 3)]


# ### 5.2: LipidMaps

# In[52]:


lm = moldb.LipidMaps()


# In[53]:


gibberellin = list(lm.get_record('LMPR0104170034', typ = 'mainkey'))[0]
gibberellin['name']


# LipidMaps too is able to yield OpenBabel objects:

# In[55]:


tag = list(lm.get_obmol('TAG(15:0_20:4_20:5)', 'synonym'))[0]
tag.draw()


# ## 6: Lipid name parser

# In order to make the databases computationally useful and use them as a combined database, we need to process their nomenclature. The `lipyd.name` module is able to recognize dozens of lipid names used in SwissLipids and LipidMaps.

# In[67]:


from lipyd import name


# In[68]:


nameproc = name.LipidNameProcessor(database = 'swisslipids', iso = True)


# In[70]:


processed_name = nameproc.process(
    ['Phosphatidylethanolamine (16:0/20:4(5Z,8Z,11Z,14Z))']
)


# In[71]:


pprint.pprint(processed_name)


# It understands even greek names:

# In[76]:


nameproc.process(['eicosapentaenoate'])


# ## 7: Combined molecule database

# In[56]:


db = moldb.MoleculeDatabaseAggregator()


# Either exact masses or adducts can be searched in the database by `lookup` and `adduct_lookup` methods, respectively.

# In[60]:


result = db.adduct_lookup(757.549011, ionmode = 'pos')


# Take a closer look at one of the resulted records:

# In[62]:


pprint.pprint(result['[M+NH4]+'][1][2])


# The exact masses and errors for all hits are also provided. Errors in ppm:

# In[65]:


result['[M+NH4]+'][2]


# Repeat the lookup with lower tolerance, and the items with high ppm disappear:

# In[66]:


result = db.adduct_lookup(757.549011, ionmode = 'pos', tolerance = 5)
result['[M+NH4]+'][2]


# ## 8: MS2 fragment definitions

# The fragment database provided by `lipyd.fragment` and `lipyd.fragdb` modules works similar way as `lipyd.lipid` and `lipyd.moldb`. `lipyd.fragment` contains near 100 predefined aliphatic chain derived fragments. In addition 140 headgroup derived fragments are included like for example 184 for choline.

# In[79]:


from lipyd import fragment
from lipyd import fragdb


# As an example take a look at a [Sph-NH2-OH]- fragment:

# In[85]:


sphfrag = fragment.Sph_mNH2_mOH(c = 18, u = 1)


# At this fragment type the constraints tell us which lipid varieties this fragment can be observed. In this case *dCer* and *DHCer*.

# In[86]:


sphfrag.constraints


# In[87]:


list(sphfrag)[0].charge, list(sphfrag)[0].mass


# ## 9: MS2 fragment database

# Look up a negative mode fragment m/z in the database. It results an array with mass, fragment name, fragment type, aliphatic chain type, carbon count, unsaturation and charge in each row. At neutral losses the charge is 0.

# In[89]:


fragdb.lookup_neg(283.26)


# Now let's annotate an MS2 scan with possible fragment identifications. To do this we open an example MGF file included in the module. The `lipyd.mgf` module serves MS2 scans from MGF files on demand. Btw the `lipyd.settings` module gives easy access for and control over near 100 customizable parameters.

# In[104]:


from lipyd import mgf
from lipyd import settings
mgffile = settings.get('mgf_example')
mgfreader = mgf.MgfReader(mgffile)
precursor = 590.45536 # this is a Cer-1P
idx, rtdiff = mgfreader.lookup_scan_ids(precursor)


# We found the following scans for precursor 590.455:

# In[105]:


idx


# Select a scan from the ones above and annotate its fragments:

# In[106]:


scan = mgfreader.scan_by_id(1941)
annot = fragdb.FragmentAnnotator(
    mzs = scan[:,0],
    ionmode = 'pos',
    precursor = precursor
)


# One example of the annotations. This fragment ranks 25 by intensity.

# In[112]:


pprint.pprint(list(annot)[24])


# ## 10: MS2 spectrum analysis

# The `lipyd.ms2.Scan` class is able to perform the entire identification workflow. By an alternative constructor method it can be initialized by providing and MGF file and scan ID.

# In[114]:


from lipyd import ms2


# In[138]:


mgfname = settings.get('mgf_pos_examples')
scan_id = 3626
scan = ms2.Scan.from_mgf(mgfname, scan_id, 'pos')


# If not provided the `Scan` instance performs the database lookup of the precursor ion. Here are the results:

# In[143]:


pprint.pprint(scan.ms1_records['[M+H]+'][1][0])


# The `identify` method attempts to confirm each of the records by analysing the MS2 spectrum.

# In[145]:


identity = scan.identify()


# The results are grouped by lipid species and come with a score. Hex2-Cer(t42:2) got a score of 45, which is the highest at this scan:

# In[148]:


pprint.pprint(identity['Hex2-Cer(t42:2)'][0])


# At the same time there were attempts to confirm for example Hex-Cer(d53:9-2OH) but it resulted a score of 0.

# In[151]:


pprint.pprint(identity['Hex-Cer(d53:9-2OH)'][0])


# Let's see one more example.

# In[152]:


mgfname = settings.get('mgf_neg_examples')
scan_id = 2516
scan = ms2.Scan.from_mgf(mgfname, scan_id, 'neg')
identity = scan.identify()


# We see that this is a PI(34:1) with score 11 and both fatty acyl fragments are confirmed by [FA-H]- ions (see the `ChainIdentificationDetails` object). These fragments are the 1st and 2nd most abundant with relative intensities of 100% and 99%.

# In[154]:


pprint.pprint(identity['PI(34:1)'][0])

