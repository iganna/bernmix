# BernMix

Computation of PMF and CDF for a weighted sum of dependent and independent ni.d. Bernoulli random variables

**Abbreviations**

BRV - Bernoulli Random variable  
PMF -  Probability mass function  
CDF - Cumulative distribution function  
ni.d. - independent non-identically distributed 


## Description

The BernMix package includes two efficient algorithms to calculate the exact distribution of a weighted sum of ni.d. BRV – the first is for integer weights and the second is for non-integer weights. The discussed distribution includes, as particular cases, Binomial and Poisson Binomial distributions together with their linear combinations. 

For integer weights, we present two algorithms to calculate a probability mass function (PMF) and a cumulative distribution function (CDF) of a weighted sum of BRVs: utilising the Discrete Fourier transform and convolution method. For non-integer weights we suggest the heuristic approach to compute a pointwise CDF using rounding and integer linear programming. We also propose the transformation algorithm to estimate the joint distribution of weighted sums of BRVs and suggest the heuristics to estimate PMF and CDF, when BRVs are non-independent. Together with numerical studies we discuss possible application in bioinformatics analysis.
  
The BernMix package provides Python implementations of all developed algorithms; C++ library for using Fast Fourier transform is wrapped with Cython. Code is available on GitHub and via PyPi. 

## Implemented methods

* `pmf_int` - computation of the PMF for a **integer-weighted** sum of BRVs by the FFT-based} method
* `pmf_int_conv` - computation of the PMF for a **integer-weighted** sum of BRVs by the **convolution** method
* `pmf_int_bf` - computation of the PMF for a integer-weighted sum of BRVs by the brute-force search
* `pmf_int_dep` - computation of the PMF for a integer-weighted sum of **dependent BRVs** by the convolution method
* `pmf_int_joint` - computation of the **joint** PMF for integer-weighted sums of BRVs
* `pmf_int_prod` - computation of the PMF for the **product** of two integer-weighted sums of BRVs
* `cdf_int` - computation of the CDF for a integer-weighted sum of BRVs by the FFT-based method
* `cdf_double` - computation of the pointwise corrected CDF for a weighted sum of BRVs with **real weights** by the FFT-based method
* `cdf_permut` - computation of the CDF for a weighted sum of BRVs by the permutation

There are two required parameters in each function: a list of success probabilities for BRVs and a vector of weights. Other parameters have default values but can be set by the user.



## Requirements

To run BernMix methods you need Python 3.4 or later. A list of required Python packages that the BernMix depends on, are in `requirements.txt`.  
The BernMix also required the [FFTW3](http://www.fftw.org/download.html) library (a C library for computing the discrete Fourier transform) and Cython.

## Installation

BernMix can be installed via PyPi:
```
pip install bernmix
```

or by the following commands:
```
git clone https://github.com/iganna/bernmix.git
cd bernmix
python setup.py sdist bdist_wheel
cd dist
pip install *.whl
```

## Running the tests

To demonstrate the use of methods we created a Python notebook `tests/bernmix_demo.ipynb`.  
All tests that were used in the below article, are presented in a Python notebook `tests/bernmix_test.ipynb` and in a R notebook `tests/gpb_test.ipynb`

## References

The mathematical inference of the algorithm implemented in the BernMix package is described in A.A.Igolkina, *On distributions of weighted sums of binary random variables*

## Authors

**Anna Igolkina** developed the BernMix package, [e-mail](mailto:igolkinaanna11@gmail.com).    

**Max Kovalev**  contributed in `bernmix_int/bernmix_fourier.c`.


## License information

The BernMix package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).
