
Tutorial
========

1: Chemical calculator
----------------------

.. code:: python

    import imp
    import itertools
    import pprint
    from lipyd import pprint_namedtuple
    pprint.PrettyPrinter = pprint_namedtuple.PrettyPrinter
    from lipyd import mass
    from lipyd import formula

The ``mass`` module knows exact masses of isotopes, isotopic abundances,
weights, etc. The ``MassBase`` class is able to process chemical
formula, calculate masses and do arithmetics.

.. code:: python

    mass.get_mass('Na')




.. parsed-literal::

    22.989769282



.. code:: python

    mass.MassBase('C2H6') - mass.MassBase('H') + mass.MassBase('OH')




.. parsed-literal::

    46.04186481376



Expression evaluation with ``mass.expr``:

.. code:: python

    mass.expr('C2H6 - H + OH')




.. parsed-literal::

    46.04186481376



.. code:: python

    mass.expr('C6H12O6 - water')




.. parsed-literal::

    162.0528234236



Make a deuterium:

.. code:: python

    mass.expr('H + n')




.. parsed-literal::

    2.0164899481400003



Hydrogensulphate ion:

.. code:: python

    mass.expr('HSO4 + e')




.. parsed-literal::

    96.96010326786924



Additional attributes can be provided in keyword arguments to carry
metadata.

.. code:: python

    lactic_acid = formula.Formula('CH3CHOHCOOH', name = 'lactic acid')
    lactic_acid.attrs.name




.. parsed-literal::

    'lactic acid'



A galactose:

.. code:: python

    ((2 * formula.Formula('C6H12O6')) - 'H2O').formula




.. parsed-literal::

    'C12H22O11'



2: Calculations with adducts
----------------------------

.. code:: python

    from lipyd import mz

This is an oleic acid. What is the mass of the [M-H]- adduct?

.. code:: python

    formula.Formula('C18H34O2').remove_h()




.. parsed-literal::

    281.24860387047



We have seen a mass and wondering what is the exact mass if it is an
[M+NH4]+ adduct:

.. code:: python

    mz.Mz(854.576482597791).remove_nh4()




.. parsed-literal::

    836.5426570444603



Calculate the [M+Li]+ adduct for the same molecule:

.. code:: python

    mz.Mz(836.5426570444603).adduct(mass.MassBase('Li', charge = 1))




.. parsed-literal::

    843.558111907551



3: Metabolite model
-------------------

Metabolites consist of a core and optionally substituents. Substituents
might be formulas or moieties with aliphatic chains.

.. code:: python

    from lipyd import metabolite
    from lipyd import substituent

Make all combinations of halogenated methanes:

.. code:: python

    halo_methanes = metabolite.AbstractMetabolite(
        core = 'C',
        subs = [('H', 'F', 'Cl', 'Br', 'I')] * 4
    )

Check the first 3:

.. code:: python

    [(m.formula, m.mass) for m in halo_methanes][:3]




.. parsed-literal::

    [('C1H4', 16.031300129039998),
     ('C1F1H3', 34.02187826038),
     ('C1Cl1H3', 49.99232782678)]



Do the same with all alcohols up to 1-8 carbon count with 0-2
unsaturated bonds:

.. code:: python

    chain = metabolite.AbstractSubstituent(c = (1, 8), u = (0, 2))
    alcohols = metabolite.AbstractMetabolite(subs = (chain, ('OH',)))
    [(m.formula, m.mass) for m in alcohols][:3]




.. parsed-literal::

    [('C1H4O1', 32.02621474924),
     ('C2H6O1', 46.04186481376),
     ('C3H8O1', 60.057514878279996)]



Make some ceramides:

.. code:: python

    # fatty acyls of length 16, 18 or 20 and one or no unsaturation:
    fattyacyl = substituent.FattyAcyl(c = (16, 18, 20), u = (0, 1))
    lcb = substituent.Sphingosine(c = 18, u = 1)
    ceramides = metabolite.AbstractMetabolite(core = 'H', subs = (lcb, fattyacyl), name = 'Cer')
    # name, formula, mass, [M+H]+, [M+NH4]+
    [(cer.name, cer.formula, cer.mass, cer.add_nh4()) for cer in ceramides]




.. parsed-literal::

    [('Cer(18:1/16:0)', 'C34H67N1O3', 537.51209502622, 555.5459205795507),
     ('Cer(18:1/16:1)', 'C34H65N1O3', 535.4964449617, 553.5302705150308),
     ('Cer(18:1/18:0)', 'C36H71N1O3', 565.54339515526, 583.5772207085907),
     ('Cer(18:1/18:1)', 'C36H69N1O3', 563.52774509074, 581.5615706440708),
     ('Cer(18:1/20:0)', 'C38H75N1O3', 593.5746952843, 611.6085208376307),
     ('Cer(18:1/20:1)', 'C38H73N1O3', 591.55904521978, 609.5928707731108)]



4: Lipid definitions
--------------------

In the ``lipyd.lipid`` module more than 150 lipid varieties are
predefined.

.. code:: python

    from lipyd import lipid

.. code:: python

    d_ceramides = lipid.CeramideD(
        sph_args = {'c': 18, 'u': (0, 1)},
        fa_args  = {'c': (14, 24), 'u': (0, 4), 'even': True},
    )

.. code:: python

    [(cer.name, cer.mass) for cer in d_ceramides][:3]




.. parsed-literal::

    [('Cer(DH18:0/14:0)', 511.49644496170004),
     ('Cer(DH18:0/14:1)', 509.48079489718003),
     ('Cer(DH18:0/14:2)', 507.46514483266003)]



5: External databases
---------------------

The ``lipyd.moldb`` module provides access to *SwissLipis* and
*LipidMaps*. Automatically downloads and processes the databases which
you can search by various names and identifiers, you can also access
structures as OpenBabel objects, InChI and SMILE strings.

.. code:: python

    from lipyd import moldb
    from lipyd import name
    from lipyd.lipproc import *

5.1: SwissLipids
~~~~~~~~~~~~~~~~

.. code:: python

    swl = moldb.SwissLipids(levels = {'species'})


.. parsed-literal::

            Indexing SwissLipids -- finished: 100%|██████████| 623M/623M [00:13<00:00, 47.4Mit/s]


Get phosphatidylethanolamines as OpenBabel objects:

.. code:: python

    swl.reload()
    pe = swl.get_hg_obmol('PE')
    pe0 = next(pe)
    pe0.draw()

.. code:: python

    [m.title for m in itertools.islice(pe, 0, 3)]




.. parsed-literal::

    ['Phosphatidylethanolamine (28:1)',
     'Phosphatidylethanolamine (39:3)',
     'Phosphatidylethanolamine (O-29:1)']



5.2: LipidMaps
~~~~~~~~~~~~~~

.. code:: python

    lm = moldb.LipidMaps()


.. parsed-literal::

    	:: Indexed 42556 records from `cache/LMSDFDownload12Dec17FinalAll.sdf`.


.. code:: python

    gibberellin = list(lm.get_record('LMPR0104170034', typ = 'mainkey'))[0]
    gibberellin['name']




.. parsed-literal::

    {'COMMON_NAME': 'gibberellin A17',
     'SYNONYMS': 'GA17; gibberellin 17',
     'PUBCHEM_CID': '5460657',
     'CHEBI_ID': '24236',
     'INCHI': 'InChI=1S/C20H26O7/c1-10-8-18-9-19(10,27)7-4-11(18)20(16(25)26)6-3-5-17(2,15(23)24)13(20)12(18)14(21)22/h11-13,27H,1,3-9H2,2H3,(H,21,22)(H,23,24)(H,25,26)/t11-,12-,13-,17-,18+,19+,20-/m1/s1'}



LipidMaps too is able to yield OpenBabel objects:

.. code:: python

    tag = list(lm.get_obmol('TAG(15:0_20:4_20:5)', 'synonym'))[0]
    tag.draw()

6: Lipid name parser
--------------------

In order to make the databases computationally useful and use them as a
combined database, we need to process their nomenclature. The
``lipyd.name`` module is able to recognize dozens of lipid names used in
SwissLipids and LipidMaps.

.. code:: python

    from lipyd import name

.. code:: python

    nameproc = name.LipidNameProcessor(database = 'swisslipids', iso = True)

.. code:: python

    processed_name = nameproc.process(
        ['Phosphatidylethanolamine (16:0/20:4(5Z,8Z,11Z,14Z))']
    )

.. code:: python

    pprint.pprint(processed_name)


.. parsed-literal::

    (
      Headgroup(main='PE', sub=()),
      ChainSummary(
        c = 36,
        u = 4,
        typ = ('FA', 'FA'),
        attr = (
            ChainAttr(sph='', ether=False, oh=()),
         ChainAttr(sph='', ether=False, oh=())
            ),
        iso = None
      ),
      (
        Chain(
          c = 16,
          u = 0,
          typ = 'FA',
          attr = ChainAttr(sph='', ether=False, oh=()),
          iso = ()
        ),
       Chain(
          c = 20,
          u = 4,
          typ = 'FA',
          attr = ChainAttr(sph='', ether=False, oh=()),
          iso = ('5Z', '8Z', '11Z', '14Z')
        )
        )
      )


It understands even greek names:

.. code:: python

    nameproc.process(['eicosapentaenoate'])




.. parsed-literal::

    (Headgroup(main='FA', sub=()),
     ChainSummary(c=20, u=5, typ=('FA',), attr=(ChainAttr(sph='', ether=False, oh=()),), iso=None),
     [Chain(c=20, u=5, typ='FA', attr=ChainAttr(sph='', ether=False, oh=()), iso=())])



7: Combined molecule database
-----------------------------

.. code:: python

    db = moldb.MoleculeDatabaseAggregator()


.. parsed-literal::

            Indexing SwissLipids -- finished: 100%|██████████| 623M/623M [00:13<00:00, 16.0Mit/s]


.. parsed-literal::

    	:: Indexed 42556 records from `cache/LMSDFDownload12Dec17FinalAll.sdf`.


.. parsed-literal::

            Generating metabolites -- finished: 100%|██████████| 44.0/44.0 [00:09<00:00, 4.50it/s]
            Generating metabolites -- finished: 100%|██████████| 18.0/18.0 [00:05<00:00, 4.10it/s]
            Generating metabolites -- finished: 100%|██████████| 106/106 [00:09<00:00, 11.4it/s] 
            Generating metabolites -- finished: 100%|██████████| 1.00/1.00 [00:00<00:00, 11.9it/s]
            Generating metabolites -- finished: 100%|██████████| 1.00/1.00 [00:00<00:00, 48.9it/s]


Either exact masses or adducts can be searched in the database by
``lookup`` and ``adduct_lookup`` methods, respectively.

.. code:: python

    result = db.adduct_lookup(757.549011, ionmode = 'pos')

Take a closer look at one of the resulted records:

.. code:: python

    pprint.pprint(result['[M+NH4]+'][1][2])


.. parsed-literal::

    LipidRecord(
      lab = LipidLabel(db_id=None, db='lipyd.lipid', names=('PC(33:4)',)),
      hg = Headgroup(main='PC', sub=()),
      chainsum = ChainSummary(
          c = 33,
          u = 4,
          typ = ('FA', 'FA'),
          attr = (
              ChainAttr(sph='', ether=False, oh=()),
          ChainAttr(sph='', ether=False, oh=())
              ),
          iso = None
        ),
      chains = ()
    )


The exact masses and errors for all hits are also provided. Errors in
ppm:

.. code:: python

    result['[M+NH4]+'][2]




.. parsed-literal::

    array([-1.69750819e-02, -1.69750819e-02, -2.69161082e-02, -2.69161082e-02,
           -2.69161082e-02, -2.69161082e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.04974547e-02,
           -3.04974547e-02, -3.04974547e-02, -3.04974547e-02, -3.23033863e+00,
           -3.23033863e+00, -3.23033863e+00, -1.11723392e+01,  1.73995818e+01,
            1.73995818e+01,  1.73995818e+01])



Repeat the lookup with lower tolerance, and the items with high ppm
disappear:

.. code:: python

    result = db.adduct_lookup(757.549011, ionmode = 'pos', tolerance = 5)
    result['[M+NH4]+'][2]




.. parsed-literal::

    array([-0.01697508, -0.01697508, -0.02691611, -0.02691611, -0.02691611,
           -0.02691611, -0.03049745, -0.03049745, -0.03049745, -0.03049745,
           -0.03049745, -0.03049745, -0.03049745, -0.03049745, -0.03049745,
           -0.03049745, -0.03049745, -0.03049745, -0.03049745, -0.03049745,
           -0.03049745, -0.03049745, -0.03049745, -0.03049745, -0.03049745,
           -0.03049745, -0.03049745, -0.03049745, -0.03049745, -0.03049745,
           -0.03049745, -3.23033863, -3.23033863, -3.23033863])



8: MS2 fragment definitions
---------------------------

The fragment database provided by ``lipyd.fragment`` and
``lipyd.fragdb`` modules works similar way as ``lipyd.lipid`` and
``lipyd.moldb``. ``lipyd.fragment`` contains near 100 predefined
aliphatic chain derived fragments. In addition 140 headgroup derived
fragments are included like for example 184 for choline.

.. code:: python

    from lipyd import fragment
    from lipyd import fragdb

As an example take a look at a [Sph-NH2-OH]- fragment:

.. code:: python

    sphfrag = fragment.Sph_mNH2_mOH(c = 18, u = 1)

At this fragment type the constraints tell us which lipid varieties this
fragment can be observed. In this case *dCer* and *DHCer*.

.. code:: python

    sphfrag.constraints




.. parsed-literal::

    (FragConstraint(hg='Cer', family=None, sub=None, sph='d', oh=0, chaintype='Sph'),
     FragConstraint(hg='Cer', family=None, sub=None, sph='DH', oh=0, chaintype='Sph'))



.. code:: python

    list(sphfrag)[0].charge, list(sphfrag)[0].mass




.. parsed-literal::

    (-1, 263.2380392001693)



9: MS2 fragment database
------------------------

Look up a negative mode fragment m/z in the database. It results an
array with mass, fragment name, fragment type, aliphatic chain type,
carbon count, unsaturation and charge in each row. At neutral losses the
charge is 0.

.. code:: python

    fragdb.lookup_neg(283.26)




.. parsed-literal::

    array([[283.2642539494093, '[FA(18:0)-H]-', 'FA-H', 'FA', 18, 0, -1],
           [283.2642539494093, '[Sph(20:0)-C2H4-NH2-2H]-', 'Sph-C2H4-NH2-2H',
            'Sph', 20, 0, -1]], dtype=object)



Now let’s annotate an MS2 scan with possible fragment identifications.
To do this we open an example MGF file included in the module. The
``lipyd.mgf`` module serves MS2 scans from MGF files on demand. Btw the
``lipyd.settings`` module gives easy access for and control over near
100 customizable parameters.

.. code:: python

    from lipyd import mgf
    from lipyd import settings
    mgffile = settings.get('mgf_example')
    mgfreader = mgf.MgfReader(mgffile)
    precursor = 590.45536 # this is a Cer-1P
    idx, rtdiff = mgfreader.lookup_scan_ids(precursor)

We found the following scans for precursor 590.455:

.. code:: python

    idx




.. parsed-literal::

    array([1941, 1929,  427,  423,  589,  645,  642,  308,  481,  478,  586,
            368,  696,  535,  532,  755,  752,  700,  721])



Select a scan from the ones above and annotate its fragments:

.. code:: python

    scan = mgfreader.scan_by_id(1941)
    annot = fragdb.FragmentAnnotator(
        mzs = scan[:,0],
        ionmode = 'pos',
        precursor = precursor
    )

One example of the annotations. This fragment ranks 25 by intensity.

.. code:: python

    pprint.pprint(list(annot)[24])


.. parsed-literal::

    (
      FragmentAnnotation(
        mz = 228.23219101229077,
        name = '[Sph(14:0)-H2O+H]+',
        fragtype = 'Sph-H2O+H',
        chaintype = 'Sph',
        c = 14,
        u = 0,
        charge = 1
      ),
      FragmentAnnotation(
        mz = 228.23219101229077,
        name = '[FA(14:0)+NH2-O]+',
        fragtype = 'FA+NH2-O',
        chaintype = 'FA',
        c = 14,
        u = 0,
        charge = 1
      ),
      FragmentAnnotation(
        mz = 228.23219101229077,
        name = '[Sph(14:1)-O+H]+',
        fragtype = 'Sph-O+H',
        chaintype = 'Sph',
        c = 14,
        u = 1,
        charge = 1
      )
      )


10: MS2 spectrum analysis
-------------------------

The ``lipyd.ms2.Scan`` class is able to perform the entire
identification workflow. By an alternative constructor method it can be
initialized by providing and MGF file and scan ID.

.. code:: python

    from lipyd import ms2

.. code:: python

    mgfname = settings.get('mgf_pos_examples')
    scan_id = 3626
    scan = ms2.Scan.from_mgf(mgfname, scan_id, 'pos')

If not provided the ``Scan`` instance performs the database lookup of
the precursor ion. Here are the results:

.. code:: python

    pprint.pprint(scan.ms1_records['[M+H]+'][1][0])


.. parsed-literal::

    LipidRecord(
      lab = LipidLabel(db_id=None, db='lipyd.lipid', names=('Hex2-Cer(t42:2)',)),
      hg = Headgroup(main='Cer', sub=('Hex2',)),
      chainsum = ChainSummary(
          c = 42,
          u = 2,
          typ = ('Sph', 'FA'),
          attr = (
              ChainAttr(sph='t', ether=False, oh=()),
          ChainAttr(sph='', ether=False, oh=())
              ),
          iso = None
        ),
      chains = ()
    )


The ``identify`` method attempts to confirm each of the records by
analysing the MS2 spectrum.

.. code:: python

    identity = scan.identify()

The results are grouped by lipid species and come with a score.
Hex2-Cer(t42:2) got a score of 45, which is the highest at this scan:

.. code:: python

    pprint.pprint(identity['Hex2-Cer(t42:2)'][0])


.. parsed-literal::

    MS2Identity(
      score = 45,
      hg = Headgroup(main='Cer', sub=('Hex2',)),
      chainsum = ChainSummary(
          c = 42,
          u = 2,
          typ = ('Sph', 'FA'),
          attr = (
              ChainAttr(sph='t', ether=False, oh=()),
          ChainAttr(sph='', ether=False, oh=())
              ),
          iso = None
        ),
      chains = (
          Chain(
            c = 18,
            u = 1,
            typ = 'Sph',
            attr = ChainAttr(sph='t', ether=False, oh=()),
            iso = ()
          ),
        Chain(
            c = 24,
            u = 1,
            typ = 'FA',
            attr = ChainAttr(sph='', ether=False, oh=()),
            iso = ()
          )
          ),
      details = ChainIdentificationDetails(rank = (0, None), i = (1.0, None), fragtype = ('Sph-2xH2O-H', None))
    )


At the same time there were attempts to confirm for example
Hex-Cer(d53:9-2OH) but it resulted a score of 0.

.. code:: python

    pprint.pprint(identity['Hex-Cer(d53:9-2OH)'][0])


.. parsed-literal::

    MS2Identity(
      score = 0,
      hg = Headgroup(main='Cer', sub=('Hex',)),
      chainsum = ChainSummary(
          c = 53,
          u = 9,
          typ = ('Sph', 'FAOH'),
          attr = (
              ChainAttr(sph='d', ether=False, oh=()),
          ChainAttr(sph='', ether=False, oh=('2OH',))
              ),
          iso = None
        ),
      chains = (
          Chain(
            c = 17,
            u = 1,
            typ = 'Sph',
            attr = ChainAttr(sph='d', ether=False, oh=()),
            iso = ()
          ),
        Chain(
            c = 36,
            u = 8,
            typ = 'FAOH',
            attr = ChainAttr(sph='', ether=False, oh=('2OH',)),
            iso = ()
          )
          ),
      details = ChainIdentificationDetails(
          rank = (1, None),
          i = (0.18172325900094663, None),
          fragtype = ('Sph-2xH2O+H', None)
        )
    )


Let’s see one more example.

.. code:: python

    mgfname = settings.get('mgf_neg_examples')
    scan_id = 2516
    scan = ms2.Scan.from_mgf(mgfname, scan_id, 'neg')
    identity = scan.identify()

We see that this is a PI(34:1) with score 11 and both fatty acyl
fragments are confirmed by [FA-H]- ions (see the
``ChainIdentificationDetails`` object). These fragments are the 1st and
2nd most abundant with relative intensities of 100% and 99%.

.. code:: python

    pprint.pprint(identity['PI(34:1)'][0])


.. parsed-literal::

    MS2Identity(
      score = 11.0,
      hg = Headgroup(main='PI', sub=()),
      chainsum = ChainSummary(
          c = 34,
          u = 1,
          typ = ('FA', 'FA'),
          attr = (
              ChainAttr(sph='', ether=False, oh=()),
          ChainAttr(sph='', ether=False, oh=())
              ),
          iso = None
        ),
      chains = (
          Chain(
            c = 18,
            u = 1,
            typ = 'FA',
            attr = ChainAttr(sph='', ether=False, oh=()),
            iso = ()
          ),
        Chain(
            c = 16,
            u = 0,
            typ = 'FA',
            attr = ChainAttr(sph='', ether=False, oh=()),
            iso = ()
          )
          ),
      details = ChainIdentificationDetails(rank = (0, 1), i = (1.0, 0.9897746748655405), fragtype = ('FA-H', 'FA-H'))
    )

