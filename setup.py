#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from setuptools import setup, find_packages

def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()


def find_version():
    with open('ssc/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Cannot find version information')
    return version

setup(
    name='ssc',
    version=find_version(),
    author='Andrey Smelter',
    author_email='andrey.smelter@gmail.com',
    description='Library for clustering peaks into spin systems',
    keywords='ssc',
    packages=find_packages(),
    package_data={'ssc': ['bin/*']},
    platforms='Linux',
    long_description=readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
