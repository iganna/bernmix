# setup.py
from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import build_ext


PACKAGE = "bernmix"
NAME = "bernmix"
DESCRIPTION = "Methods to compute PMF and CDF values of a weighted sum of " \
              "i.ni.d. BRVs"
AUTHOR = "Anna Igolkina"
AUTHOR_EMAIL = "igolkinaanna11@gmail.com"
URL = "https://github.com/iganna/bernmix"
VERSION = "1.11.44"

#with open('requirements.txt') as f:
#    reqs = f.read().splitlines()
reqs = ['cvxopt', 'Cython', 'numpy', 'scipy']

with open("README.md", "r") as fh:
    long_description = fh.read()


src_dir = 'bernmix/bernmix_int'
ext = Extension(src_dir + '.bernmix_int',
                [src_dir + '/bernmix_int.pyx',
                 src_dir + '/bernmix_fourier.c'],
                libraries=['fftw3'],
                extra_compile_args=["-lfftw3"])
#
# src_dir2 = 'bernmix/bernmix_fancy'
# ext2 = Extension(src_dir2 + '.bm_dynam_fast',
#                  [src_dir2 + '/bm_dynam_fast.pyx'])

package_data = {'bernmix_int': ['*.pxd', '*.pyx', '*.c', '*.h'],
                'bm_dynam_fast': ['*.pyx', '*.c']}


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license="MIT",
    url=URL,
    packages=find_packages(),
    install_requires=reqs,
    # ext_modules=[ext, ext2],
    ext_modules=[ext],
    cmdclass={'build_ext': build_ext},
    package_data = package_data,
    include_package_data=True,
    long_description=long_description,
    long_description_content_type="text/markdown", 
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    zip_safe=False,
)
