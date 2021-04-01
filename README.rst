.. image:: https://github.com/micom-dev/micom/raw/master/docs/source/micom.png
    :width: 640 px

|actions status| |coverage| |pypi status|

Welcome
-------

`MICOM` is a Python package for metabolic modeling of microbial
communities currently developed in the
`Gibbons Lab <https://gibbons.systemsbiology.org>`_ at the `Institute for Systems
Biology <https://systemsbiology.org>`_ and the
`Human Systems Biology Group <https://resendislab.github.io>`_ of
Prof. Osbaldo Resendis Antonio at the `National Institute of Genomic
Medicine Mexico <https://inmegen.gob.mx>`_.

`MICOM` allows you to construct a community model from a list on input
COBRA models and manages exchange fluxes between individuals and individuals
with the environment. It explicitly accounts for different abundances of
individuals in the community and can thus incorporate data from 16S rRNA
sequencing experiments. It allows optimization with a variety of algorithms
modeling the trade-off between egoistic growth rate maximization and
cooperative objectives.

Attribution
-----------

MICOM is published in

::

      MICOM: Metagenome-Scale Modeling To Infer Metabolic Interactions in the Gut Microbiota
      Christian Diener, Sean M. Gibbons, Osbaldo Resendis-Antonio
      mSystems 5:e00606-19
      https://doi.org/10.1128/mSystems.00606-19

Please cite this publication when referencing MICOM. Thanks :smile:

Installation
------------

`MICOM` is available on PyPi and can be installed via

.. code:: bash

    pip install micom

Getting started
---------------

Documentation can be found at https://micom-dev.github.io/micom .

Getting help
------------

General questions on usage can be asked in Github Discussions
    https://github.com/micom-dev/micom/discussions

We are also available on the cobrapy Gitter channel
    https://gitter.im/opencobra/cobrapy

Questions specific to the MICOM Qiime2 plugin (q2-micom) can also be asked on the Qiime2 forum
    https://forum.qiime2.org/c/community-plugin-support/


.. |actions status| image:: https://github.com/micom-dev/micom/workflows/Python%20package/badge.svg
   :target: https://github.com/micom-dev/micom/actions
.. |coverage| image:: https://codecov.io/gh/micom-dev/micom/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/micom-dev/micom
.. |pypi status| image:: https://img.shields.io/pypi/v/micom.svg
   :target: https://pypi.org/project/micom/
