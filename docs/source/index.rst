.. micom documentation master file, created by
   sphinx-quickstart on Thu Mar  9 13:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://github.com/micom-dev/micom/raw/master/docs/source/_static/micom.png
    :width: 640 px

|actions status| |coverage| |pypi status|

`MICOM` is a Python package for metabolic modeling of microbial
communities developed in the
`Gibbons Lab <https://gibbons.systemsbiology.org>`_ at the `Institute for Systems
Biology <https://systemsbiology.org>`_ and the
`Human Systems Biology Group <https://resendislab.github.io>`_ of
Prof. Osbaldo Resendis Antonio at the `National Institute of Genomic
Medicine Mexico <https://inmegen.gob.mx>`_.

`MICOM` allows you to construct a community model from a list on input
COBRA models and manages exchange fluxes between individuals and individuals
with the environment. It explicitly accounts for different abundances of
individuals in the community and can thus incorporate data from 16S amplicon
sequencing or metagenomic experiments. It allows optimization with a variety of algorithms
modeling the trade-off between egoistic growth rate maximization and
cooperative objectives.

Attribution
-----------

MICOM is published at https://msystems.asm.org/content/5/1/e00606-19.
Please cite this article when referencing MICOM.



To get an idea which assumptions and strategies MICOM uses we recommend
to start with some background on the :doc:`methods <logic>`.

Contents
--------
.. toctree::
   :maxdepth: 1

   Methods used by MICOM <logic>
   Installing MICOM <installing>
   The High Level API <high_level>
   Building communities <community>
   Growth rates and fluxes <growth_fluxes>
   Growth media <media>
   Knockouts <taxa_knockouts>
   Intervention studies <elasticities>
   Analyzing many samples in parallel <workflows>
   API <micom>


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |actions status| image:: https://github.com/micom-dev/micom/workflows/Python%20package/badge.svg
   :target: https://github.com/micom-dev/micom/actions
.. |coverage| image:: https://codecov.io/gh/micom-dev/micom/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/micom-dev/micom
.. |pypi status| image:: https://img.shields.io/pypi/v/micom.svg
   :target: https://pypi.org/project/micom/
