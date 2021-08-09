:py:mod:`micom.workflows.media`
===============================

.. py:module:: micom.workflows.media

.. autoapi-nested-parse::

   Example workflows for micom.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.workflows.media.process_medium
   micom.workflows.media._medium
   micom.workflows.media.minimal_media
   micom.workflows.media._fix_medium
   micom.workflows.media.fix_medium



.. py:function:: process_medium(medium, samples)

   Prepare a medium for simulation.


.. py:function:: _medium(args)

   Get minimal medium for a single model.


.. py:function:: minimal_media(manifest, model_folder, summarize=True, min_growth=0.1, threads=1)

   Calculate the minimal medium for a set of community models.


.. py:function:: _fix_medium(args)

   Get the fixed medium for a model.


.. py:function:: fix_medium(manifest, model_folder, medium, community_growth=0.1, min_growth=0.001, max_import=1, minimize_components=False, summarize=True, weights=None, threads=1)

   Augment a growth medium so all community members can grow in it.

   :param manifest: The manifest as returned by the `build` workflow.
   :type manifest: pandas.DataFrame
   :param model_folder: The folder in which to find the files mentioned in the manifest.
   :type model_folder: str
   :param medium: A growth medium with exchange reaction IDs as index and positive
                  import fluxes as values. If a DataFrame needs columns `flux` and
                  `reaction`.
   :type medium: pandas.Series or pandas.DataFrame
   :param community_growth: The minimum community-wide growth rate that has to be achieved on the created
                            medium.
   :type community_growth: positive float
   :param min_growth: The minimum biomass production required for growth.
   :type min_growth: positive float
   :param max_import: The maximum import rate for added imports.
   :type max_import: positive float
   :param minimize_components: Whether to minimize the number of media components rather than the
                               total flux.
   :type minimize_components: boolean
   :param summarize: Whether to summarize the medium across all samples. If False will
                     return a medium for each sample.
   :type summarize: boolean
   :param weights: Will scale the fluxes by a weight factor. Can either be "mass" which will
                   scale by molecular mass, a single element which will scale by
                   the elemental content (for instance "C" to scale by carbon content).
                   If None every metabolite will receive the same weight.
                   Will be ignored if `minimize_components` is True.
   :type weights: str
   :param threads: The number of processes to use.
   :type threads: int

   :returns: A new growth medium with the smallest amount of augmentations such
             that all members of the community can grow in it.
   :rtype: pandas.DataFrame


