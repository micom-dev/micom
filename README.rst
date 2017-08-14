.. image:: https://github.com/cdiener/micom/raw/master/micom.png

|travis status| |appveyor status| |coverage| |pypi status|

*Please note that this is still a pre-release and might still lack
functionalities. It comes with no guarantees.*

`micom` is a Python package for metabolic modeling of microbial
communities developed in the
[Human Systems Biology Group](https://resendislba.github.io) of
Prof. Osbaldo Resendis Antonio at the [National Institute of Genomic
Medicine Mexico](https://inmegen.gob.mx).

`micom` allows you to construct a community model from a list on input
COBRA models and manages exchange fluxes between individuals and individuals
with the environment. It explicitly accounts for different abundances of
individuals in the community and can thus incorporate data from 16S rRNA
sequencing experiments. It allows optimization with a variety of algorithms
modeling the trade-off between egoistic growth rate maximization and
cooperative objectives.

Attribution
-----------

If you want to use micom in a scientific publication the attribution clause in
the license is covered by citing our relevant publication. However, we are still
in the process in pairing micom with validations so that publication still
**does not exist**. If you want to use micom before that article is published please
contact us at `oresendis (at) inmegen.gob.mx`. Thanks :smile:

Installation
------------

`micom` is available on PyPi and can be installed via

.. code:: bash

    pip install micom

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
