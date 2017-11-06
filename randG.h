//
// Created by molin on 17-11-6.
//

#ifndef MOEA_D_RANDG_H
#define MOEA_D_RANDG_H

#include <ctime>
#include <gsl/gsl_rng.h>
//Debug scope

double randG(){
    clock_t t;
    t = clock();
    //std::cout<<t<<std::endl;
    static int counter = 3;
    gsl_rng * engine = gsl_rng_alloc(gsl_rng_taus);
    int seed = counter+t;
    gsl_rng_set( engine, seed);
    double lambda = gsl_rng_uniform(engine);
    gsl_rng_free(engine);
    counter ++;

    //@Todo: Consider parallism computing
    return lambda;
}

#endif //MOEA_D_RANDG_H
