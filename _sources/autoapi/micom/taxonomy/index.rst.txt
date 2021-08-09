:py:mod:`micom.taxonomy`
========================

.. py:module:: micom.taxonomy

.. autoapi-nested-parse::

   Helpers to convert external data to a MICOM taxonomy.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.taxonomy.build_from_qiime
   micom.taxonomy.qiime_to_micom



.. py:function:: build_from_qiime(abundance, taxonomy: pandas.Series, collapse_on='genus') -> pandas.DataFrame

   Build the specification for the community models.


.. py:function:: qiime_to_micom(feature_table, taxonomy, collapse_on='genus')

   Load a micom taxonomy from Qiime 2 data.

   :param feature_table: Path to a Qiime 2 FeatureTable artifact.
   :type feature_table: str
   :param taxonomy: Path to a Qiime 2 FeatureData[Taxonomy] artifact.
   :type taxonomy: str
   :param collapse_on: The taxa ranks to collapse on. This will dictate how strict the database
                       matching will be as well.
   :type collapse_on: str or List[str]

   :returns: A micom taxonomy containing abundances and taxonomy calls in long
             format.
   :rtype: pd.DataFrame


