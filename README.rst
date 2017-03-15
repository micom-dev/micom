.. image:: https://github.com/cdiener/micom/raw/master/micom.png

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

Installation
------------

`micom` is available on PyPi and can be installed via

.. code:: bash

    pip install micom

Note that `micom` currently requires a development version of
COBRApy (>=0.6.0a2).

Getting started
---------------

Documentation can be found at https://cdiener.github.io/micom .

.. |travis status| image:: https://travis-ci.org/cdiener/micom.svg?branch=master
   :target: https://travis-ci.org/cdiener/micom
.. |appveyor status| image:: https://ci.appveyor.com/api/projects/status/m2vu008h7n35ji2g/branch/master?svg=true
   :target: https://ci.appveyor.com/project/cdiener/micom/branch/master
.. |coverage| image:: https://codecov.io/gh/cdiener/micom/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/cdiener/micom
.. |pypi status| image:: https://img.shields.io/pypi/v/micom.svg
   :target: https://pypi.org/project/micom/
