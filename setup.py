#!/usr/bin/env python

from setuptools import setup, find_packages

version = '1.0dev'

msg = '''------------------------------
Installing omSnpScore version {}
------------------------------
'''.format(version)
print(msg)

setup(
    name='omSnpScore',
    version=version,
    author='lx Gui',
    author_email='guilixuan@gmail.com',
    keywords=['bioinformatics', 'NGS', 'Reseq', 'SNP'],
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    scripts=[
        'scripts/snpScore', 'scripts/table2hd5', 'scripts/table2pkl',
        'scripts/snpScore-cli', 'scripts/merge_vcf_feature', 'scripts/snpInf',
        'scripts/snpInf-cli'
    ],
    install_requires=[
        'click',
        'pandas',
        'loguru',
        'delegator.py',
        'pybedtools',
        'numpy',
        'tables',
        'attrs',
    ],
)

msg = '''------------------------------
omSnpScore installation complete!
------------------------------
'''
print(msg)
