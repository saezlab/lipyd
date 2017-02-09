.. emese documentation master file, created by
   sphinx-quickstart on Mon Jan 30 17:10:27 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to emese's documentation!
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

emese.main
----------

The main module contains the **Screening** class which
represents the total analysis workflow on multiple objects (proteins).
   
The Screening class
+++++++++++++++++++

.. autoclass:: emese.main.Screening
   :members:

emese.ms2
---------

Submodule with methods for MS2 spectra analysis.

The Feature class 
+++++++++++++++++

Represents one MS1 m/z detected across multiple samples (fractions).

.. autoclass:: emese.ms2.Feature
   :members:

The MS2Scan class
+++++++++++++++++

Represents one MS2 spectrum with all the methods for its analysis.

.. autoclass:: emese.ms2.MS2Scan
   :members:


emese.mz
--------

Submodule for calculations with m/z values.

.. automodule:: emese.mz
   :members:

emese.fragments
---------------

Classes representing fragment ions in MS2 spectra.

.. automodule:: emese.fragments
   :members:

emese.mass
----------

Submodule with molecular mass calculation classes.

.. automodule:: emese.mass
   :members:

emese._curl
-----------

Submodule for the built-in curl based downloader.

.. automodule:: emese._curl
   :members:

emese.common
------------

Submodule for small methods commonly used in other modules
and some generic constant data structures.

.. automodule:: emese.common
   :members:

emese.progress
--------------

In this submodule the progress bar takes place.

.. automodule:: emese.progress
   :members:
