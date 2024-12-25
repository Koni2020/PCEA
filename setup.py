from setuptools import setup

__author__ = "Hanyu Jin"
__version__ = '0.1'
__email__ = "henrychin2006@gmail.com"
__license__ = "MIT"
__copyright__ = "Copyright Hanyu Jin"

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = "pyCEA",
    version = __version__,
    author = __author__,
    author_email = __email__,
    description = ("A python package for compound event analysis."),
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/Koni2020/pyCEA",
    packages = ["pyLiang"],
    license = __license__,
    install_requires = ["numpy"],
    classifiers = [
		"Programming Language :: Python :: 3.10",
                "License :: OSI Approved :: MIT License",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Topic :: Scientific/Engineering",
    ]
)