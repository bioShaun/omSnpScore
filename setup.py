#!/usr/bin/env python

from setuptools import setup, find_packages

version = "1.0dev"

msg = """------------------------------
Installing omSnpScore version {}
------------------------------
""".format(
    version
)
print(msg)

setup(
    name="omSnpScore",
    version=version,
    author="lx Gui",
    author_email="guilixuan@gmail.com",
    keywords=["bioinformatics", "NGS", "Reseq", "SNP"],
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    scripts=[
        "scripts/publicSample",
        "scripts/table2csv-mp",
        "scripts/splitVcf-mp",
        "scripts/snpScore-bychr",
        "scripts/snpScore-mp2",
        "scripts/varDensityCompare",
        "scripts/varDensityCompare-cli",
        "scripts/snpFilter-cli",
        "scripts/snpFilter-bychr",
        "scripts/snpFilter-mp",
    ],
    install_requires=[
        "fire",
        "click",
        "pandas",
        "loguru",
        "delegator.py",
        "pybedtools",
        "numpy",
        "tables",
        "attrs",
        "jinja2",
    ],
)

msg = """------------------------------
omSnpScore installation complete!
------------------------------
"""
print(msg)
