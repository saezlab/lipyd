*lipyd* – A Python module for lipidomics LC MS/MS data analysis
===============================================================

This module implements methods and workflows for MS/MS lipidomics data
analysis. Runs primarily in Python 3 but also in Python 2.7.x.

Input and preprocessing
-----------------------

At reading raw mass spec data from mzML files, peak picking and feature
detection we rely on the `OpenMS <http://openms.de/>`_ library. This ensures
computationally efficient processing by well established methods. As our
OpenMS integration is not yet complete we provide a temporary solution to
read already preprocessed features from CSV files exported by the
`PEAKS <http://www.bioinfor.com/peaks-studio/>`_ software. We are not
comfortable with the idea of building on expensive proprietary software and
in the near future we will provide complete integration with OpenMS.

Metabolite database lookup
--------------------------

The ``lipyd.modb`` module provides an unified interface to standard
databases like `SwissLipids <https://swisslipids.org/>`_ and
`SwissLipids <http://lipidmaps.org/>`_ In addition it is able to generate
custom metabolite masses.
With the default settings the database consists of more than
100 thousands of lipid species. The ``lipyd.lipid`` module
contains more than 150 predefined lipid classes and it's easy to define
new ones. The ``Sample`` and ``SampleSet`` objects in
``lipyd.sample``, which represent a series of features, support
the automatic lookup in the databases.

MS2 spectrum identification
---------------------------

The ``lipyd.ms2`` module contains generic classes to support the
analysis and identification of MS2 spectra. Based on around 50 standards
run by our group and reviewing many spectra from publications and
databases we created bult in rules for identification of more than 80
lipid classes. You can modify the methods or create new ones by writing
Python methods. However we are working on
`MFQL <https://wiki.mpi-cbg.de/lipidx/LipidXplorer_MFQL>`_ integration to
provide a more standard way of defining rules. Also we will introduce
similarity search against spectrum databases.

Feature filtering, post-processing
----------------------------------

The ``lipyd.sample`` and ``lipyd.feature`` modules provide
classes for analysis of features optionally in relation to other variables
and filter them. Analysis and filtering of the features can be done
before or after the lipid identification. Doing it before reduces the
number of MS2 spectra to be analyzed this way saving time. In the future
we will add more utilities to build arrays of features and also MS2
fragments across arbitrary number of experiments to provide opportunities
for higher level analysis.
