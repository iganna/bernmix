# BernMix

Computation of PMF and CDF for a weighted sum of i.ni.d. Bernoulli random variables

**Abbreviations**

BRV - Bernoulli Random variable  
PMF -  Probability mass function  
CDF - Cumulative distribution function  
i.ni.d. - independent non-identically distributed 


## Description

The BernMix package includes two efficient algorithms to calculate the exact distribution of a weighted sum of i.ni.d. BRV – the first is for integer weights and the second is for non-integer weights. The discussed distribution includes, as particular cases, Binomial and Poisson Binomial distributions together with their linear combinations. For integer weights we present the algorithm to calculate a PMF and a CDF of a weighted sum of BRVs utilising the Discrete Fourier transform of the characteristic function. For non-integer weights we suggest the heuristic approach to compute pointwise CDF using rounding and integer linear programming.  
  
The BernMix package provides a Python implementation of the algorithms to calculate PMFs and CDFs for both cases (integer and non-integer weights); C++ library for using Fast Fourier transform is wrapped with Cython. We analyse the time complexity of the algorithms and demonstrate their performance and accuracy.  

## Implemented methods

* `bernmix_pmf_int` - computation of the PMF for a integer-weighted sum of BRVs by the developed method
* `bernmix_cdf_int` - computation of the CDF for a integer-weighted sum of BRVs by the developed method
* `bernmix_cdf_double` - computation of the CDF for a weighted sum of BRVs with real weights by the developed method
* `conv_pmf_int` - computation of the PMF for a integer-weighted sum of BRVs by the convolution
* `permut_cdf` - computation of the CDF for a weighted sum of BRVs by the permutation


## Requirements

To run BernMix methods you need Python 3.4 or later. A list of required Python packages that the BernMix depends on, are in `requirements.txt`.  
The BernMix also required the [FFTW3](http://www.fftw.org/download.html) library (a C library for computing the discrete Fourier transform) and Cython.

## Installation


To install the BernMix package, run the following commands:
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

The mathematical inference of the algorithm implemented in the BernMix package is described in A.A.Igolkina et al., *BernMix: the distribution of a weighted sum of independent Bernoulli random variables with different success probabilities*

## Authors

**Anna Igolkina** developed the BernMix package, [e-mail](mailto:igolkinaanna11@gmail.com).    
**Max Kovalev**  contributed in `bernmix_int/bernmix_fourier.c`.


## License information

The BernMix package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).
