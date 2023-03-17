:py:mod:`micom.workflows.core`
==============================

.. py:module:: micom.workflows.core

.. autoapi-nested-parse::

   Makes it easier to run analyses on several samples in parallel.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   micom.workflows.core.save_results
   micom.workflows.core.load_results
   micom.workflows.core.workflow



Attributes
~~~~~~~~~~

.. autoapisummary::

   micom.workflows.core.GrowthResults


.. py:data:: GrowthResults

   

.. py:function:: save_results(results, path)

   Save growth results to a file.

   This will write all tables as CSV into a single ZIP file.

   :param results: The results as returned from `grow`.
   :type results: GrowthResults
   :param path: A filepath for the generated file. Should end in `.zip`.
   :type path: str


.. py:function:: load_results(path)

   Load growth results from a file.

   :param path: Path to saved `GrowthResults`.
   :type path: str

   :returns: The saved GrowthResults.
   :rtype: GrowthResults


.. py:function:: workflow(func, args, threads=4, description=None, progress=True)

   Run analyses for several samples in parallel.

   This will analyze several samples in parallel. Includes a workaround for
   optlang memory leak.

   :param func: A function that takes a single argument (can be any object) and
                that performs your analysis for a single sample.
   :type func: function
   :param args: An array-like object (list, tuple, numpy array, pandas Series, etc.)
                that contains the arguments for each sample.
   :type args: array-like object
   :param threads: How many samples to analyze in parallel at once.
   :type threads: positive int
   :param description: The dewscription shown in front of the progress bar.
   :type description: str
   :param progress: Whether to show a progress bar.
   :type progress: bool


