
import os
import imp

import ppeval

rootdir = '/'
basedir = os.path.join(rootdir, 'home', 'denes', 'documents')
msproc_basedir = os.path.join(basedir, 'lipyd', 'src')


p = ppeval.PeakPickingEvaluation(
    feature_xml_fname = os.path.join(
        msproc_basedir,
        'lipyd_ms_preproc',
        '150310_Popeye_MLH_AC_STARD10_A10_pos',
        '150310_Popeye_MLH_AC_STARD10_A10_pos__features.featureXML',
    ),
    peaks_fname = 'stard10_pos.csv',
    examples_fname = 'STARD10_pos_examples.csv',
    outfile = 'STARD10_ppeval__iterative__no-smoothing.tsv',
    progenesis_fname = 'STARD10_invivo_pos_A10_isotopes.xml',
)

p.main()


p = ppeval.PeakPickingEvaluation(
    feature_xml_fname = os.path.join(
        msproc_basedir,
        'lipyd_ms_preproc',
        '150310_Popeye_MLH_AC_STARD10_A10_pos',
        '150310_Popeye_MLH_AC_STARD10_A10_pos__features.featureXML',
    ),
    peaks_fname = 'stard10_pos.csv',
    examples_fname = 'STARD10_pos_examples.csv',
    outfile = 'STARD10_ppeval__iterative__gaussian-20.tsv',
    progenesis_fname = 'STARD10_invivo_pos_A10_isotopes.xml',
)

p.main()




p = ppeval.PeakPickingEvaluation(
    feature_xml_fname = os.path.join(
        msproc_basedir,
        'lipyd_ms_preproc',
        '150310_Popeye_MLH_AC_STARD10_A10_pos',
        '150310_Popeye_MLH_AC_STARD10_A10_pos__features.featureXML',
    ),
    peaks_fname = 'stard10_pos.csv',
    examples_fname = 'STARD10_pos_examples.csv',
    outfile = 'STARD10_ppeval__iterative__gaussian-10.tsv',
    progenesis_fname = 'STARD10_invivo_pos_A10_isotopes.xml',
)

p.main()
