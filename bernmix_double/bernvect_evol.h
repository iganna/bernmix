//
//  gan_alg.h
//  
//
//  Created by Anna Igolkina on 02/10/2017.
//
//

#ifndef bernvect_evol_h
#define bernvect_evol_h

#include <stdio.h>
// unique algorithm example
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <random>
#include <algorithm>
#include <iterator>
#include <numeric>


void evolving (std::vector<std::vector<int> >& bernvectors,
               std::vector<double>& fitness,
               std::vector<double> weights,
               double mut_prob,
               int n_epoch);


#endif /* bernvect_evol_h */
