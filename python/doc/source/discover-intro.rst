=============================
 An introduction to DISCOVER
=============================

.. currentmodule:: discover
.. default-domain:: py


Introduction
============

DISCOVER is a novel statistical test for detecting co-occurrence and mutual exclusivity in cancer genomics data. Unlike traditional approaches used for these tasks such as Fisher's exact test, DISCOVER is based on a statistical model that takes into account the overall tumour-specific alteration rates when deciding whether alterations co-occur more or less often than expected by chance. This improved model prevents spurious associations in co-occurrence detection, and increases the statistical power to detect mutual exclusivities.

This vignette introduces the DISCOVER package for mutual exclusivity and co-occurrence analysis. Two variants of the DISCOVER test are implemented by this package. The first variant tests pairs of genes for either co-occurrence or mutual exclusivity. The second is a mutual exclusivity test for gene sets larger than two. We will illustrate the use of both tests by performing a mutual exclusivity and co-occurrence analysis of somatic mutations in the TCGA breast cancer samples.


Set up
======

First, we load the package and the breast cancer mutation data.


.. ipython::

   In [1]: import discover

   In [2]: import discover.datasets

   In [3]: mut = discover.datasets.load_dataset("brca_mut")
    

A small fragment of the mutation matrix is shown below. It takes the form of a binary matrix, i.e. only 0 and 1 are allowed as values.


.. ipython::

   In [4]: mut.iloc[:5, :5]
   Out[4]:
           TCGA-A1-A0SB-01A-11D-A142-09  TCGA-A1-A0SD-01A-11D-A10Y-09  \
   A1CF                               0                             0
   A2M                                0                             0
   A2ML1                              0                             0
   A4GALT                             0                             0
   A4GNT                              0                             0
   
           TCGA-A1-A0SE-01A-11D-A099-09  TCGA-A1-A0SF-01A-11D-A142-09  \
   A1CF                               0                             0
   A2M                                0                             0
   A2ML1                              0                             0
   A4GALT                             0                             0
   A4GNT                              0                             0
   
           TCGA-A1-A0SG-01A-11D-A142-09
   A1CF                               0
   A2M                                0
   A2ML1                              0
   A4GALT                             0
   A4GNT                              0



An important ingredient of the DISCOVER test is the estimation of a background matrix. This is how tumour-specific alteration rates can be used by the test. To estimate this background matrix, a :class:`DiscoverMatrix` object is created and passed the mutation matrix. Depending on the size of the matrix this might take some time. This estimation step only needs to be performed once.


.. ipython::

   In [5]: events = discover.DiscoverMatrix(mut)
    



Note that for estimating this background matrix, it is important that the full mutation matrix is provided. Even if only a subset of genes will subsequently be used in the analysis, a whole-genome view of the mutations is required for this first step. Indeed, for this example we will only look for mutual exclusivity between genes that are mutated in at least 25 tumours.


.. ipython::

   In [6]: subset = mut.sum(1) > 25


Pairwise tests
==============

We can now test for mutual exclusivity between all pairs of genes using the :func:`pairwise_discover_test` function. Its only required argument is an instance of the :class:`DiscoverMatrix` class like we created above. In order to only analyse frequently mutated genes, we select the corresponding rows before passing the matrix to the function. Again, depending on the size of the matrix, this function might take some time.


.. ipython::

   In [7]: result_mutex = discover.pairwise_discover_test(events[subset])


We can get a quick overview of the results of this analysis simply by printing the result object.


.. ipython::

   In [8]: result_mutex
   Out[8]:
   Pairwise DISCOVER mutual exclusivity test
   alternative hypothesis: observed overlap is less than expected by chance
   
   number of pairs tested: 3403
   proportion of true null hypotheses: 0.928959115562
   number of significant pairs at a maximum FDR of 0.01: 24  


..
	The `print` method for the `pairwise.discover.out` class optionally takes an argument with which the false discovery rate (FDR) threshold is selected. By default, a maximum FDR of 1% is used. Below, we get a summary of the test results using a somewhat more liberal threshold of 5%.

	print(result.mutex, fdr.threshold=0.05)


To get the pairs of genes which are significantly mutually exclusive, we can use the :meth:`~discover.pairwise.PairwiseDiscoverResult.significant_pairs` method. This method, takes an optional FDR threshold as argument. Below, we use the default of 1%.



.. ipython::

   In [9]: result_mutex.significant_pairs()
   Out[9]:
        gene1   gene2        pvalue        qvalue
   0   ARID1A    TP53  1.707339e-04  3.720371e-03
   1     CDH1    FAT3  2.367758e-04  4.229179e-03
   2     CDH1   GATA3  1.013593e-05  2.693358e-04
   3     CDH1  MAP3K1  1.953768e-04  4.005492e-03
   4     CDH1    TP53  6.351581e-17  5.110803e-16
   5    CSMD1  MAP3K1  5.207218e-04  9.568275e-03
   6      DMD  MAP3K1  1.392113e-04  3.208038e-03
   7     FAT3  PIK3CA  3.152191e-04  5.572603e-03
   8    GATA3  MAP3K1  1.278231e-04  3.089085e-03
   9    GATA3   MUC16  4.622371e-04  8.737364e-03
   10   GATA3    MUC4  6.633634e-05  1.711992e-03
   11   GATA3  PIK3CA  1.149704e-06  2.738015e-05
   12   GATA3   SYNE1  9.833305e-05  2.462260e-03
   13   GATA3    TP53  5.470534e-15  3.138128e-14
   14   GATA3     TTN  1.943835e-04  4.005492e-03
   15   GATA3   XIRP2  2.429469e-04  4.229179e-03
   16  MAP2K4    TP53  4.315389e-05  1.111645e-03
   17  MAP3K1   MUC16  4.783053e-04  8.737364e-03
   18  MAP3K1   NCOR1  1.527654e-05  3.835730e-04
   19  MAP3K1     NEB  2.431831e-04  4.229179e-03
   20  MAP3K1    TP53  1.207007e-10  8.509564e-10
   21  MAP3K1     TTN  4.145708e-05  1.111645e-03
   22    MUC4   MUC5B  2.043579e-04  4.048278e-03
   23  PIK3CA    TP53  1.224574e-10  8.509564e-10


Groupwise test
==============

In the above results, we can identify a set of four genes---*TP53*, *CDH1*, *GATA3*, and *MAP3K1*---for which all pairwise combinations are found mutually exclusive. We can use the :func:`groupwise_discover_test` to asses whether this set of genes is also mutually exclusive as a group. This we can do by simply passing the subset of the :class:`DiscoverMatrix` that corresponds to the given genes.


.. ipython::

   In [10]: genes = ["TP53", "CDH1", "GATA3", "MAP3K1"]

   In [11]: discover.groupwise_discover_test(events[genes])
   Out[11]: 7.819884565135893e-22
