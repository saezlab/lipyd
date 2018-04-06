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
import sys
import os
from setuptools import setup, find_packages
import imp

_version = imp.load_source('_version', os.path.join('src', 'lipyd', '_version.py'))
__version__ = _version.__version__

metainfo = {
    'authors': {
        'Türei':('Dénes Türei','turei.denes@gmail.com'),
    },
    'version': __version__,
    'license': 'GPLv3',
    'download_url': ['https://github.com/saezlab/lipyd/archive/dev3.zip'],
    'url': ['https://github.com/saezlab/lipyd'],
    'description': 'Analysis of lipidomics mass spectrometry data',
    'platforms': [
        'Linux',
        'Unix',
        'MacOSX',
        'Windows'
    ],
    'keywords': [
        'lipidomics',
        'lipids',
        'mass spectrometry',
        'MS',
        'MS2',
        'lipid-protein interactions',
        'lipid transfer proteins',
        'glycerophospholipids',
        'sphingolipids',
        'metabolomics',
        'metabolites',
        'SwissLipids',
        'LipidMaps',
        'molecular weight',
        'chemical formula'
    ],
    'classifiers': [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Natural Language :: English',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry'
    ]
}

with open('README.rst') as f:
    readme = f.read()
with open('HISTORY.rst') as f:
    history = f.read()

deps = [
    'numpy',
    'scipy',
    'matplotlib',
    'fastcluster',
    'pycurl',
    'lxml',
    'xlrd',
    'httplib2',
    'seaborn',
    'rpy2',
    'sklearn',
    'xlrd',
    'openpyxl',
    'xlsxwriter',
    'future',
    'tqdm'
]

setup(
    name = 'lipyd',
    version = metainfo['version'],
    maintainer = metainfo['authors']['Türei'][0],
    maintainer_email = metainfo['authors']['Türei'][1],
    author = metainfo['authors']['Türei'][0],
    author_email = metainfo['authors']['Türei'][1],
    long_description = readme + '\n\n' + history,
    keywords = metainfo['keywords'],
    description = metainfo['description'],
    license = metainfo['license'],
    platforms = metainfo['platforms'],
    url = metainfo['url'],
    download_url = metainfo['download_url'],
    classifiers = metainfo['classifiers'],
    # package installation
    package_dir = {'':'src'},
    packages = list(set(find_packages() + ['lipyd'])),
    include_package_data = True,
    install_requires = deps
    # dependency_links = deplinks
)
