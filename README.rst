==========
 DISCOVER
==========

DISCOVER is a novel statistical method for analysing co-occurrence and mutual exclusivity in cancer genomics data. The details of this method are described in our paper *A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence* (`Genome Biology 2016 17:261`_).

.. _`Genome Biology 2016 17:261`: http://dx.doi.org/10.1186/s13059-016-1114-x


Installation
============

DISCOVER is available for both Python and R.


Python
------

The easiest way to install the DISCOVER Python package is by using Miniconda_ or Anaconda_. We provide precompiled DISCOVER packages for 64-bit Linux, Windows, and Mac OS X. The following steps assume that Miniconda or Anaconda has been installed.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Anaconda: https://www.continuum.io/why-anaconda

To create a new conda environment containing DISCOVER and its dependencies, execute the following command in a terminal (Linux/Mac OS X) or command prompt (Windows).

::

  conda create -n discover -c http://ccb.nki.nl/software/discover/repos/conda discover

On Linux or Mac OS X, this environment can be activated using:

::

  source activate discover

On Windows, the correct command is:

::

  activate discover

Note that this environment contains the bare minimum to use DISCOVER. It does not, for example, include IPython. Trying to run IPython anyway might start a version installed in a different environment, and hence, importing DISCOVER will not succeed. Consult the `conda documentation`_ to find out how to install additional packages (such as IPython). Alternatively, DISCOVER can be installed in an existing environment as follows.

::

  conda install -c http://ccb.nki.nl/software/discover/repos/conda discover

.. _conda documentation: http://conda.pydata.org/docs/using/pkgs.html

Check the documentation_ for instructions on how to use this package.


R
-

We provide precompiled R packages for Windows and Mac OS X, as well as a source package for installation on Linux. Installation on Linux requires gfortran version 4.6 or later. To install the DISCOVER package, execute the following in an R session.

::

  options(repos=c(getOption("repos"), "http://ccb.nki.nl/software/discover/repos/r"))
  install.packages("discover")

Check the documentation_ for instructions on how to use this package.


Documentation
=============

* Python_
* R_

.. _Python: http://ccb.nki.nl/software/discover/doc/python
.. _R: http://ccb.nki.nl/software/discover/doc/r/discover-intro.html
