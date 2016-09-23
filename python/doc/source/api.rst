.. currentmodule:: discover
.. _api:

===============
 API Reference
===============

.. _api.functions:

.. autosummary::
   :toctree:

   DiscoverMatrix
   pairwise_discover_test
   ~discover.pairwise.PairwiseDiscoverResult
   row_stack


Alteration matrices
===================

.. autoclass:: DiscoverMatrix
   :members:

.. autofunction:: row_stack


Pairwise co-occurrence and mutual exclusivity tests
===================================================

.. autofunction:: pairwise_discover_test

.. autoclass:: discover.pairwise.PairwiseDiscoverResult
   :members:


Groupwise mutual exclusivity tests
==================================

.. autofunction:: groupwise_discover_test
