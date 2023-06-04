#!/usr/bin/env python
from lib.__version__ import version

from setuptools import setup, find_packages
from distutils.extension import Extension

with open('README.md') as f:
    long_description = f.read()


setup(
    name='OGAP',
    version=version,
    description='OGAP: Organelle Genome Annotation Pipeline',
    url='https://github.com/zhangrengang/OGAP/',
    author='Zhang, Ren-Gang and Wang, Zhao-Xuan',
    license='GPL-3.0',

    python_requires='==2.7:',
    packages=find_packages(),
    include_package_data=True,
    scripts=[],
    entry_points={
        'console_scripts': ['ogap = OGAP:main',
		 'ogap-makedb = makedb:main',
		'ogap-checkdb = checkdb:main', 
		'ogap-compare = compare:main',
        ],
    },
)
