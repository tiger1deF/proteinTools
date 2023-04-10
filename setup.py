import setuptools
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(
    name = "proteinTools",
    version = "0.1.8",
    author = "Christian de Frondeville",
    description = "Lightweight package which simplifies interacting with proteins.",
    long_description = long_description,
    packages = ["proteinTools"],
    

)