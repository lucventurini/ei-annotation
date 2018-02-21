from distutils.core import setup
from distutils.extension import Extension
import os

USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("eiannot.util.bam2hints", sources=[os.path.join("eiannot", "util", "bam2hints.pyx")], language="c++")]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    ext_modules = extensions
)
