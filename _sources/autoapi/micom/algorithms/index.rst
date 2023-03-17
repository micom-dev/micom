:py:mod:`micom.algorithms`
==========================

.. py:module:: micom.algorithms

.. autoapi-nested-parse::

   Implements additional analysis algorithms for communities.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.algorithms.jaccard
   micom.algorithms.euclidean
   micom.algorithms.reaction_matrix
   micom.algorithms.metabolic_dist



.. py:function:: jaccard(inclusion)

   Calculate jaccard distances for a community.


.. py:function:: euclidean(inclusion)

   Calculate euclidean distances for a community.


.. py:function:: reaction_matrix(files)

   Create a matrix of reactions x models.


.. py:function:: metabolic_dist(reactions, metric=jaccard)

   Calculate the metabolic distances between all members.


