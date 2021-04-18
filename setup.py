#!/usr/bin/env python

from setuptools import setup, find_packages

version = "1.0.3"

msg = """------------------------------
Installing varBscore version {}
------------------------------
""".format(
    version
)
print(msg)

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="varBscore",
    version=version,
    description="A python package to calculate varBscore using vcf input",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioShaun/omSnpScore",
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
        "scripts/vcfValidator",
    ],
    install_requires=[
        "fire",
        "click",
        "pandas",
        "loguru",
        "delegator.py",
        "pybedtools",
        "numpy",
        "attrs",
        "jinja2",
        "typer",
        "vcfpy",
    ],
)

msg = """------------------------------
varBscore installation complete!
------------------------------
"""
print(msg)
