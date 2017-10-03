//
//  gan_alg.c
//  
//
//  Created by Anna Igolkina on 02/10/2017.
//
//

#include "bernvect_evol.h"




/*!
 * This function returns 1 if vector of Bernoulli trails are euqal or not
 *
 */
static bool bernvect_cmp (std::vector<int>& x, std::vector<int>& y)
{

    // ANNA asserts
    int n = x.size();
    for(int i = 0; i < n; ++i)
    {
        if ((x.at(i) & 1) ^ (y.at(i) & 1)) //if bits on the 0-th porituin are different
            return 0; //then bernvectors are different
    }
    return 1;
}


static double vect_sum(std::vector<double> v)
{
    double s = 0;
    for(int i = 0; i < (int)v.size(); ++i)
        s += std::abs(v.at(i));
    return s;
}

/*!
 * This function generates new population of bernvectors
 * fitnef is signed value around zero
 */
static void bernvect_generate (std::vector<std::vector<int> >& bernvectors,
                        std::vector<double>& fitness,
                        std::vector<double> weights,
                        double& mut_prob)
{

    //set the random generator + bernoulli distribition
    std::default_random_engine generator;
    std::bernoulli_distribution distribution(mut_prob);
    
    // ANNA asserts
    int n_vect = fitness.size();
    int n_bern = weights.size();
    int counts;
    
    
    for(int i = 0; i < n_vect; ++i) // mutate each vector of Bernoulli trails
    {
        // coly the current vector to mutate
        bernvectors.push_back(bernvectors.at(i));
        fitness.push_back(fitness.at(i));
        counts = 0;
        for (int j = 0; j < n_bern; ++j)
        {
            
            if (distribution(generator) == 0) // no changes for j-th Bernoulli trail
                continue;

            ++counts; // flag for changes
                
            if (bernvectors.back().at(j) & 1) // if the Bernoulli value was true - remove it
            {
                fitness.back() -= weights.at(j);
                bernvectors.back().at(j) = 0;
            }
            else  // if the Bernoulli value was false - addd it
            {
                fitness.back() += weights.at(j);
                bernvectors.back().at(j)  = 1;
            }
            
        }
        if (counts == 0)
        {
            bernvectors.pop_back();
            fitness.pop_back();
            continue;
        }
        // chech for unique
        // if not unique - remove
        for (int k = 0; k < ((int)fitness.size() - 1); ++k)
        {
            if (k == i)
                continue;
            if (bernvect_cmp(bernvectors.back(), bernvectors.at(k)))
            {
                bernvectors.pop_back();
                fitness.pop_back();
                break;
            }
        }
    }


    // Remain only top n_vect vectors by the lowest p_value
    // 1. get indexes of sorted fintess
    // 2. slice n_vec indexes

    std::vector<int> index(fitness.size(), 0);
    for (int i = 0 ; i < (int)index.size() ; i++) index[i] = i;
    
    std::sort(std::begin(index), std::end(index),
              [&](const int& a, const int& b)
              {
                  return (std::abs(fitness[a]) < std::abs(fitness[b]));
              });
    

    std::vector<std::vector<int> > bernvectors_new;
    std::vector<double> fitness_new;
    for(int i = 0; i < n_vect; ++i)
    {
        bernvectors_new.push_back(bernvectors.at(index.at(i)));
        fitness_new.push_back(fitness.at(index.at(i)));
    }
    bernvectors = bernvectors_new;
    fitness = fitness_new;
    
    
}



void evolving (std::vector<std::vector<int> >& bernvectors,
                             std::vector<double>& fitness,
                             std::vector<double> weights,
                             double mut_prob,
                             int n_epoch)
{

    for(int i = 0; i < n_epoch; ++i)
    {
        bernvect_generate (bernvectors,
                           fitness,
                           weights,
                           mut_prob);
        mut_prob = mut_prob + (0.5-mut_prob)*0.1;

        if ( i % 100 == 0 )
        {
            
            std::cout << i << " " << vect_sum(fitness) << " " << mut_prob <<"\n";
        }
        
    }
        
}
