#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import glob
import os

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy',
    'scipy',
    'astropy',
    'matplotlib',
    'pandas',
    'hipercam'
]

test_requirements = [
    # TODO: put package test requirements here
]

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

setup(
    name='cam-cal',
    version='0.1.0',
    description="Python package for flux calibrating HiPERCAM and ULTRACAM light curves",
    long_description=readme + '\n\n' + history,
    author="Alex Brown",
    author_email='ajbrown2@shef.ac.uk',
    url='https://github.com/Alex-J-Brown/cam-cal',
    download_url='https://github.com/Alex-J-Brown/cam-cal/archive/v0.1.0.tar.gz',
    packages=[
        'cam-cal'
    ],
    package_dir={'cam-cal':
                 'cam-cal'},
    include_package_data=True,
    scripts=scripts,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='cam-cal',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ]
)
