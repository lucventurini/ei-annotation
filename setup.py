# from distutils.core import setup
from setuptools import setup, find_packages
from distutils.extension import Extension
import os
from os import path
from glob import glob

USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.c'

here = path.abspath(path.dirname("__file__"))

if os.path.exists(path.join(here, "DESCRIPTION.md")):
    with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
        long_description = description.read()
else:
    long_description = ""

version = None
with open(path.join(here, "eiannot", "__init__.py"), "rt") as fp:
    # exec(fp.read(), globals(), version)
    for line in fp:
        if line.startswith("__version__"):
            version = line.rstrip().split(" = ")[1].strip('"')
            break


extensions = [Extension("eiannot.util.bam2hints", sources=[os.path.join("eiannot", "util", "bam2hints.pyx")],
                        language="c++")]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='eiannot',
    version=version,
    description="A Python3 complete annotation pipeline",
    long_description=long_description,
    ext_modules=extensions,
    scripts=glob(os.path.join("eiannot", "util", "*py")) + glob(os.path.join("eiannot", "util", "*pl")),
    entry_points={"console_scripts": ["eiannot = eiannot.cli:main"]},
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    packages=find_packages(".")
)
