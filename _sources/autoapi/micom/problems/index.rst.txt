:py:mod:`micom.problems`
========================

.. py:module:: micom.problems

.. autoapi-nested-parse::

   Implements tradeoff optimization between community and egoistic growth.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.problems.regularize_l2_norm
   micom.problems.cooperative_tradeoff
   micom.problems.knockout_taxa



.. py:function:: regularize_l2_norm(community, min_growth)

   Add an objective to find the most "egoistic" solution.

   This adds an optimization objective finding a solution that maintains a
   (sub-)optimal community growth rate but is the closest solution to the
   community members individual maximal growth rates. So it basically finds
   the best possible tradeoff between maximizing community growth and
   individual (egoistic) growth. Here the objective is given as the sum of
   squared differences between the individuals current and maximal growth
   rate. In the linear case squares are substituted by absolute values
   (Manhattan distance).

   :param community: The community to modify.
   :type community: micom.Community
   :param min_growth: The minimal community growth rate that has to be mantained.
   :type min_growth: positive float
   :param linear: Whether to use a non-linear (sum of squares) or linear version of the
                  cooperativity cost. If set to False requires a QP-capable solver.
   :type linear: boolean
   :param max_gcs: The precomputed maximum individual growth rates.
   :type max_gcs: None or dict


.. py:function:: cooperative_tradeoff(community, min_growth, fraction, fluxes, pfba, atol, rtol)

   Find the best tradeoff between community and individual growth.


.. py:function:: knockout_taxa(community, taxa, fraction, method, progress, diag=True)

   Knockout a taxon from the community.


