# -*- coding: utf-8 -*-

from distutils.core import setup
from setuptools import find_packages

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.7
Operating System :: POSIX :: Linux
""".strip().split('\n')

setup(
    name='COSMO',
    version='0.1.0',
    packages=find_packages(exclude=['tests']),
    url='https://github.com/SANBI-SA/COSMO',
    license='GPLv3',
    author='Hocine Bendou',
    author_email='hocine@sanbi.ac.za',
    description='Tool to detect operons (relative to M. tuberculosis H37Rv) using gene expression levels.',
    keywords='Mycobacterium tuberculosis bioinformatics',
    classifiers=classifiers,
    package_dir={'COSMO': 'COSMO'},
    install_requires=[
        'pysam>=0.15.0',
    ],
    extras_require={
        'test': ['pytest>=4.3.1'],
    },
    entry_points={
        'console_scripts': [
            'COSMO = operon.user_input:main',
        ]
    }
)
