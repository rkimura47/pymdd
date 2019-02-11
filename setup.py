# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pymdd',
    version='0.0.1',
    description='A simple python implementation of an MDD',
    long_description=readme,
    author='Ryo Kimura',
    author_email='rkimura47@gmail.com',
    url='https://github.com/rkimura47/pymdd',
    license='MIT',
    packages=find_packages(exclude=('docs', 'examples'))
)
