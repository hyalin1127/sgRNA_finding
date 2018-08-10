#!/usr/bin/env python

'''
CRISPR_sgRNA_finding set up script
'''
from __future__ import print_function;

import os
import sys
from distutils.core import setup, Extension
from subprocess import call as subpcall
from distutils.command.install import install as DistutilsInstall

setup(
   name='sgRNA_finding',
   version='1.0',
   license="MIT",
   description='CRISPR_sgRNA_finding',
   author='Chen-Hao Chen',
   author_email='hyalin1127@gmail.com',
   packages=['sgRNA_finding'],
   scripts = ['bin/sgRNA_finding'],
   package_dir={'sgRNA_finding':'sgRNA_finding'},
)
