#!/usr/bin/env python

import os
from setuptools import setup

from rpXsuite._version import __version__

setup(name='RPxSuite',
      version=__version__,
      description='RPxSuite is a dedicated, stand-alone, tool for the identification and clustering of marker genes in shotgun metagenomic assemblies.',
      url='https://github.com/alexcritschristoph/RPxSuite',
      author='Matt Olm and Alex Crits-Christoph',
      author_email='crits-christoph@berkeley.edu',
      license='MIT',
      package_data={'rpXsuite': ['helper_files/essential.hmm']},
      include_package_data=True,
      packages=['rpXsuite'],
      scripts=['rpXsuite/RPxSuite.py'],
      python_requires='>=3.4.0',
      install_requires=[
          'biopython',
      ],
      zip_safe=False)
