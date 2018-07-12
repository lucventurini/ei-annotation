# from distutils.core import setup
from setuptools import setup, find_packages
from distutils.extension import Extension
import os
from glob import glob

USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("eiannot.util.bam2hints", sources=[os.path.join("eiannot", "util", "bam2hints.pyx")],
                        language="c++")]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='eiannot',
    version="0.0.2",
    ext_modules=extensions,
    scripts=glob(os.path.join("eiannot", "util", "*py")) + glob(os.path.join("eiannot", "util", "*pl")),
    entry_points={"console_scripts": ["eiannot = eiannot.cli:main"]},
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    packages=find_packages(".")
)
