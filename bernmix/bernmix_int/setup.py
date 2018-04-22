from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


ext = Extension("bernmix_int",
                sources=["bernmix_int.pyx", "bernmix_fourier.c"])

setup(
      ext_modules=cythonize([ext]))
