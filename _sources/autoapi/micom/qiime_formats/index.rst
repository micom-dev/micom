:py:mod:`micom.qiime_formats`
=============================

.. py:module:: micom.qiime_formats

.. autoapi-nested-parse::

   Provides support for Qiime formats.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.qiime_formats.metadata
   micom.qiime_formats.load_qiime_model_db
   micom.qiime_formats.load_qiime_manifest
   micom.qiime_formats.load_qiime_model
   micom.qiime_formats.load_qiime_medium
   micom.qiime_formats.load_qiime_feature_table
   micom.qiime_formats.load_qiime_taxonomy



Attributes
~~~~~~~~~~

.. autoapisummary::

   micom.qiime_formats.yaml
   micom.qiime_formats._has_manifest


.. py:data:: yaml

   

.. py:data:: _has_manifest
   :value: ['CommunityModels[Pickle]', 'MetabolicModels[JSON]', 'MetabolicModels[SBML]']

   

.. py:function:: metadata(artifact)

   Read metadata from a Qiime 2 artifact.


.. py:function:: load_qiime_model_db(artifact, extract_path)

   Prepare a model database for use.


.. py:function:: load_qiime_manifest(artifact)

   Prepare community models for use.


.. py:function:: load_qiime_model(artifact, id)

   Load a model from a Qiime 2 artifact.


.. py:function:: load_qiime_medium(artifact)

   Load a growth medium/diet from a Qiime 2 artifact.


.. py:function:: load_qiime_feature_table(artifact)

   Load a feature table from a Qiime 2 artifact.


.. py:function:: load_qiime_taxonomy(artifact)

   Load taxonomy feature data from a Qiime 2 artifact.


