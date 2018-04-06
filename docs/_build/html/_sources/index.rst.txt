.. lipyd documentation master file, created by
   sphinx-quickstart on Mon Jan 30 17:10:27 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to lipyd's documentation!
=================================

.. toctree::
   :maxdepth: 5 
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Reference
=========

lipyd.main
----------

The main module contains the **Screening** class which
represents the total analysis workflow on multiple objects (proteins).
   
The Screening class
+++++++++++++++++++

.. autoclass:: lipyd.main.Screening
   :members:

lipyd.mz
--------

Submodule for calculations with m/z values.

.. automodule:: lipyd.mz
   :members:


lipyd.mass
----------

Submodule with molecular mass calculation classes.

.. automodule:: lipyd.mass
   :members:

lipyd.formula
-------------

Defines the `Formula` class which represents a chemical formula.

.. automodule:: lipyd.formula
   :members:

lipyd.metabolite
----------------

Abstract base classes for representation of metabolites with
various constant cores and substituent groups.

.. automodule:: lipyd.metabolite
   :members:

lipyd.substituent
-----------------

Classes representing certain substituents e.g. `FattyAcyl`.

.. automodule:: lipyd.substituent
   :members:

lipyd.lipid
-----------

Many classes representing generic or more specific lipids.
These are all based on `metabolite.AbstractMetabolite`.

.. automodule:: lipyd.lipid
   :members:

lipyd.moldb
-----------

Handling of molecular mass databases. Here you find `SwissLipids` and
`LipidMaps` and also the `MoleculeDatabaseAggregator` which you can
use to combine all the databases and also to build your own database
by auto-generating masses either by the module defaults or your custom
rules.

.. automodule:: lipyd.moldb
   :members:

lipyd.sdf
---------

Parsing `sdf` files.

.. automodule:: lipyd.sdf
   :members:

lipyd.mgf
---------

Parsing `mgf` files.

.. automodule:: lipyd.mgf
   :members:

lipyd.ms2
---------

Submodule with methods for MS2 spectra analysis.

The Feature class 
+++++++++++++++++

Represents one MS1 m/z detected across multiple samples (fractions).

.. autoclass:: lipyd.ms2.Feature
   :members:

The MS2Scan class
+++++++++++++++++

Represents one MS2 spectrum with all the methods for its analysis.

.. autoclass:: lipyd.ms2.MS2Scan
   :members:

lipyd.fragments
---------------

Classes representing fragment ions in MS2 spectra.

.. automodule:: lipyd.fragments
   :members:

lipyd._curl
-----------

Submodule for the built-in curl based downloader.

.. automodule:: lipyd._curl
   :members:

lipyd.common
------------

Submodule for small methods commonly used in other modules
and some generic constant data structures.

.. automodule:: lipyd.common
   :members:

lipyd.progress
--------------

In this submodule the progress bar takes place.

.. automodule:: lipyd.progress
   :members:
