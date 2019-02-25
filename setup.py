# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from glob import glob
from os.path import basename, splitext

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pymdd',
    version='0.0.3',
    description='A simple python implementation of an MDD',
    long_description=readme,
    author='Ryo Kimura',
    author_email='rkimura47@gmail.com',
    url='https://github.com/rkimura47/pymdd',
    license='MIT',
    packages=find_packages('src', exclude=('docs', 'examples', 'tests')),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')]
)
