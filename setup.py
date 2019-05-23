#!/usr/bin/env python

import os
from setuptools import setup

from rpXsuite._version import __version__

setup(name='RPxSuite',
      version=__version__,
      description='RPxSuite is a dedicated, stand-alone, tool for the identification and clustering of marker genes in shotgun metagenomic assemblies.',
      url='https://github.com/alexcritschristoph/RPxSuite',
      author='Matt Olm, Alex Crits-Christoph, and Spencer Diamond',
      author_email='crits-christoph@berkeley.edu',
      license='MIT',
      package_data={'rpXsuite': ['helper_files/essential.hmm', 'helper_files/SupplementalTable_S2.2.tsv',
                    'helper_files/SupplementalTable_S4.2.tsv']},
      include_package_data=True,
      packages=['rpXsuite'],
      scripts=['rpXsuite/RPxSuite.py'],
      python_requires='>=3.4.0',
      install_requires=[
          'biopython',
      ],
      zip_safe=False)
