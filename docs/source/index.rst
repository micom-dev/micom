.. micom documentation master file, created by
   sphinx-quickstart on Thu Mar  9 13:00:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to micom
================

|travis status| |appveyor status| |coverage| |pypi status|

*Please note that this is still a pre-release and might still lack
functionalities.*

`micom` is a Python package for metabolic modeling of microbial
communities.

`micom` allows you to construct a community model from a list on input
COBRA models and manages exchange fluxes between individuals and individuals
with the environment. It explicitly accounts for different abundances of
individuals in the community and can thus incorporate data from 16S rRNA
sequencing experiments. It allows optimization with a variety of algorithms
modeling the trade-off between egoistic growth rate maximization and
cooperative objectives.

Contents
--------
.. toctree::
   :maxdepth: 1

   Installing micom <installing>
   Building communities <community>
   API <micom>


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |travis status| image:: https://travis-ci.org/cdiener/micom.svg?branch=master
   :target: https://travis-ci.org/cdiener/micom
.. |appveyor status| image:: https://ci.appveyor.com/api/projects/status/m2vu008h7n35ji2g/branch/master?svg=true
   :target: https://ci.appveyor.com/project/cdiener/micom/branch/master
.. |coverage| image:: https://codecov.io/gh/cdiener/micom/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/cdiener/micom
.. |pypi status| image:: https://img.shields.io/pypi/v/micom.svg
   :target: https://pypi.org/project/micom/
