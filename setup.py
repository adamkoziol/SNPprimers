#!/usr/bin/env python
from setuptools import setup, find_packages
__author__ = 'adamkoziol'

setup(
    name="SNPprimers",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description=u'Extracts nucleotide sequence data surrounding SNP coordinates in .fasta assembly files.'
                u'Create primers to amplify extracted regions using Primer3',
    url='https://github.com/adamkoziol/metagenomeFilter',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'xlsxwriter'],
)
