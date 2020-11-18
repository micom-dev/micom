"""setup.py for micom."""

from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="micom",
    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version="0.19.0",
    description="Microbial community modeling based on cobrapy.",
    long_description=long_description,
    # The project's main homepage.
    url="https://github.com/micom-dev/micom",
    # Author details
    author="Christian Diener",
    author_email="mail@cdiener.com",
    # Choose your license
    license="Apache License 2.0",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: Apache Software License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    # What does your project relate to?
    keywords="microbiota modeling metabolism community",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=["docs", "tests"]),
    package_data={
        "micom": [
            path.join("data", "*.csv"),
            path.join("data", "*.gz"),
            path.join("data", "templates", "*.*"),
            path.join("data", "artifacts", "*.*")
        ]
    },
    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "cobra>=0.17.1",
        "optlang>=1.4.4",
        "pandas>=0.20.0",
        "tqdm>=4.14",
        "loguru>=0.3.2",
        "jinja2>=2.10.0",
        "scikit-learn>=0.22.0",
        "scipy>=1.0.0",
        "fastcluster>=1.1.0",
        "symengine>=0.6.1"
    ],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        "dev": [""],
        "test": ["coverage", "pytest", "pytest-cov"],
    },
)
