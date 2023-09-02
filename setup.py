import setuptools
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(
    name = "proteinTools",
    version = "1.9.1",
    author = "Christian de Frondeville",
    description = "Lightweight, object-oriented backend bioinformatics package which simplifies interacting with proteins.",
    long_description = long_description,
    packages = ["proteinTools"],
    install_requires=['pandas','urllib3','mygene', 'pubchempy', 'stringdb'] #'chembl-webresource-client',
)
