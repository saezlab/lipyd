
Installation
============

:note: ``lipyd`` is developed and tested in Python 3.7. Although it is written
       a Python 2.7 compatible way, no thorough testing has been done.
       We recommend to use Python 3.7+ whenever it's possible. Please file
       issues if you are experiencing problems in Python 2.7.

Linux
-----

In any modern Linux distribution installation by ``pip`` should be
straightforward.

.. code:: bash
    
   pip install git+https://github.com/saezlab/lipyd.git



Or download the tarball:

.. code:: bash
    
   curl -L https://github.com/saezlab/lipyd/archive/master.zip -o lipyd.zip
   pip install lipyd.zip

Mac OS X
--------

Should simply work by ``pip``. You might need to install ``pycurl``. Always be
aware which distributions do you have on your computer and which one are you
installing to.

Windows
-------

You need to install ``pycurl`` from this unofficial repository:
https://www.lfd.uci.edu/~gohlke/pythonlibs/#pycurl

For example, for x86_64 architecture and Python 3.7 do this:

.. code:: bash
    
    pip install https://download.lfd.uci.edu/pythonlibs/h2ufg7oq/pycvodes-0.9.1-cp37-cp37m-win_amd64.whl

Then install by ``pip`` as usual.
