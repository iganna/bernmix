# BernMix
## Computation of PMF and CDF for a weighted sum of i.ni.d. Bernoulli random variables

PMF --  Probability mass function  
CDF -- Cumulative distribution function  
i.ni.d. -- independent non-identically distributed  


## Description
We consider two efficient algorithms to calculate the exact distribution of a weighted sum of i.ni.d. Bernoulli random variables (BRVs) – the first is for integer weights and the second is for non-integer weights. The discussed distribution includes, as particular cases, Binomial and Poisson Binomial distributions together with their linear combinations. For integer weights we present the algorithm to calculate a PMF and a CDF of a weighted sum of BRVs utilising the Discrete Fourier transform of the characteristic function. For non-integer weights we suggest the heuristic approach to compute pointwise CDF using rounding and integer linear programming.  
The bernmix package provides a Python implementation of the algorithms to calculate PMFs and CDFs for both cases (integer and non-integer weights); C++ library for using Fast Fourier transform is wrapped with Cython. We analyse the time complexity of the algorithms and demonstrate their performance and accuracy.

## Requirements
Separate file


## Installation
To install, run the following commands:
```
git clone https://github.com/iganna/bernmix.git
cd bernmix
python setup.py install
```


## Running the tests
Python notebook


## References
- The mathematical inference of the algorithm implemented in BernMix is described in XXX


## License
BernMix package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).