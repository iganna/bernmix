# setup.py
from setuptools import setup


PACKAGE = "bernmix"
NAME = "bernmix"
DESCRIPTION = "Methods to compute PMF and CDF values of a weighted sum of " \
              "i.ni.d. BRVs"
AUTHOR = "Anna Igolkina"
AUTHOR_EMAIL = "igolkinaanna11@gmail.com"
URL = "https://github.com/iganna/bernmix"
VERSION = "0.1"

with open('requirements.txt') as f:
    reqs = f.read().splitlines()


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license="MIT",
    url=URL,
    packages=[PACKAGE],
    install_requires=reqs,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    zip_safe=False,
)
