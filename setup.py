#!/usr/bin/env python
import sys

version = sys.argv[1]
del sys.argv[1]
from setuptools import setup

with open('README.md') as f:
    readme = f.read()

setup(
    name='MSEE_MgO_UQ',
    version=version,
    long_description_content_type="text/markdown",
    long_description=readme,
    description="Code developed to analytically model shock wave experiments with uncertainty",
    url='https://github.com/dimtsap/MSEE_MgO_UQ',
    author="Dimitris Tsapetis",
    license='MIT',
    platforms=["OSX", "Windows", "Linux"],
    install_requires=[
        "UQpy"
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
    ],
)
