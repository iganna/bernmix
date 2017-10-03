from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize



ext = Extension("bernvect",
                sources=["bernvect.pyx", "bernvect_evol.cpp"],
                language="c++",
                extra_compile_args=["-std=c++11"],
                extra_link_args=["-std=c++11"]
                )

setup(
      ext_modules = cythonize([ext]))
