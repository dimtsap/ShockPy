#!/usr/bin/env python
import sys

version = sys.argv[1]
del sys.argv[1]
from setuptools import setup, find_packages

setup(
    name='ImpedancePy',
    version=version,
    description="Code developed to analytically model shock wave experiments with uncertainty",
    url='https://github.com/dimtsap/MSEE_MgO_UQ',
    author="Dimitris Tsapetis",
    license='MIT',
    platforms=["OSX", "Windows", "Linux"],
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "UQpy"
    ],
    extras_require={
      'dev': [
          "Sphinx",
          "sphinx-rtd-theme"
      ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
    ],
)
