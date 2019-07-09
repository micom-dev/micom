.. micom documentation master file, created by
   sphinx-quickstart on Thu Mar  9 13:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to micom
================


|travis status| |appveyor status| |coverage| |pypi status|

`micom` is a Python package for metabolic modeling of microbial
communities developed in the
`Human Systems Biology Group <https://resendislab.github.io>`_ of
Prof. Osbaldo Resendis Antonio at the `National Institute of Genomic
Medicine Mexico <https://inmegen.gob.mx>`_ and the
`Gibbons Lab <https://gibbons.systemsbiology.org>`_ at the `Institute for Systems
Biology <https://systemsbiology.org>`_.

`micom` allows you to construct a community model from a list on input
COBRA models and manages exchange fluxes between individuals and individuals
with the environment. It explicitly accounts for different abundances of
individuals in the community and can thus incorporate data from 16S rRNA
sequencing experiments. It allows optimization with a variety of algorithms
modeling the trade-off between egoistic growth rate maximization and
cooperative objectives.

Attribution
-----------

Micom is described in a `recent preprint <https://doi.org/10.1101/361907>`_ that you can cite.
If you have any questions contact us at `oresendis (at) inmegen.gob.mx`. Thanks :smile:

To get an idea which assumptions and strategies `micom` uses we recommend
to start with some background on the :doc:`methods <logic>`.

Contents
--------
.. toctree::
   :maxdepth: 1

   Methods used by micom <logic>
   Installing micom <installing>
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

.. |travis status| image:: https://travis-ci.org/micom-dev/micom.svg?branch=master
   :target: https://travis-ci.org/micom-dev/micom
.. |appveyor status| image:: https://ci.appveyor.com/api/projects/status/uqcmw82uq9jtui0t?svg=true
   :target: https://ci.appveyor.com/project/cdiener/micom-uicdk
.. |coverage| image:: https://codecov.io/gh/micom-dev/micom/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/micom-dev/micom
.. |pypi status| image:: https://img.shields.io/pypi/v/micom.svg
   :target: https://pypi.org/project/micom/
