#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

__revision__ = "$Id$"
import os
from setuptools import setup, find_packages
import imp

_version = imp.load_source(
    '_version',
    os.path.join('src', 'lipyd', '_version.py')
)
__version__ = _version.__version__

with open('README.rst') as f:
    readme = f.read()

with open('HISTORY.rst') as f:
    history = f.read()

setup(
    name = 'lipyd',
    version = __version__,
    maintainer = 'Dénes Türei',
    maintainer_email = 'turei.denes@gmail.com',
    author = 'Dénes Türei',
    author_email = 'turei.denes@gmail.com',
    long_description = readme + '\n\n' + history,
    keywords = [
        'lipidomics',
        'lipids',
        'mass spectrometry',
        'MS',
        'MS2',
        'LC MS/MS',
        'lipid-protein interactions',
        'lipid transfer proteins',
        'glycerophospholipids',
        'sphingolipids',
        'metabolomics',
        'metabolites',
        'SwissLipids',
        'LipidMaps',
        'molecular weight',
        'chemical formula',
        'chemistry',
    ],
    description = 'Analysis of lipidomics mass spectrometry data',
    license = 'GPLv3',
    platforms = [
        'Linux',
        'Unix',
        'MacOSX',
        'Windows',
    ],
    url = ['https://saezlab.github.io/lipyd'],
    download_url = ['https://github.com/saezlab/lipyd/archive/master.zip'],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    package_dir = {'': 'src'},
    packages = [
        'lipyd',
        'lipyd.reader',
        'lipyd.data',
        'lipyd.data.ms2_examples',
    ],
    include_package_data = True,
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'fastcluster',
        'pycurl',
        'lxml',
        'xlrd',
        'httplib2',
        'seaborn',
        'sklearn',
        'xlrd',
        'openpyxl',
        'xlsxwriter',
        'future',
        'tqdm',
        'xmltodict',
        'bs4',
        'regex',
    ],
    extras_require = {
        'tests': [
            'pytest',
        ],
        'docs': [
            'sphinx',
        ]
    }
)
