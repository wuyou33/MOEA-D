//
//  instance.cpp
//  MOEA/D
//
//  Created by Molin on 04/08/2017.
//  Copyright © 2017 merlin. All rights reserved.
//
//
#include "loader.h"
#include "randG.h"
#include "moead.h"
#include <iostream>
using namespace std;

int main(){
    std::string FILE_PATH = "/home/molin/CLionProjects/MOEA-D/DataSet/portreb1.txt";
    //std::string FILE_PATH = "/Users/mirror/ClionProjects/MOEA-D/DataSet/portreb1.txt";
    struct Constraint port1_constraint;
    double port1_correlations[31][31] = {};
    vector <struct asset> assetArray;
    bool s[4] = {0,0,0,0};
    for(int i =0; i<4; i++){
        char buffer;
        //buffer = getchar();
        //s[i] = (buffer =='1'?true:false);
    }

    //
    //Check input
    //

    if(loadItem(  FILE_PATH,
                assetArray,
                port1_constraint,
                port1_correlations)) {
        util_preprocess(assetArray);
        /*
        for (int i = 0; i < port1_constraint.num_assets; i++) {
            cout << assetArray[i].holding<<"\t"<<
                 assetArray[i].min_buy<<"\t"<<
                 assetArray[i].max_buy<< "\t"<<
                 assetArray[i].max_sell<<"\t"<<
                 assetArray[i].mean_income<<endl;
        }
         */

    }

    setting(s, assetArray);
    //
    //Check random number genrator
    //
    /*
    for(int i = 0; i<10; i++){
        cout<<randG()<<endl;
    }
     */
    //
    //Check init_lamb
    //
    vector <lamb> lamblist;
    if(init_lamb(2, 99, lamblist)){
        for(auto item:lamblist){
            //cout<<item.v[0]<<"\t"<<item.v[1]<<endl;
        }
    }
    init_distance(lamblist);
    //display_nn_id(lamblist);
    //display_nn(lamblist);

    //1.3 Initialization
    population x;
    init_population(x, assetArray, port1_constraint, port1_correlations);
    population ep;
    ep = x;
    //1.4 Initialize Z solution
    double z_population[2];
    double max_income = 0;
    double min_risk = x.xi[0].fitness[1];
    for(int i = 0; i<x.xi.size(); i++){
        if(x.xi[i].fitness[0]>max_income){
            max_income = x.xi[i].fitness[0];
        }
        if(x.xi[i].fitness[1]<min_risk){
            min_risk = x.xi[i].fitness[1];
        }
    }
    z_population[0] = max_income;
    z_population[1] = min_risk;
    cerr<<"[1.4]:\tInitialize z pupulation:\t"
        <<"\n\t\tMax income:\t"<<z_population[0]<<"\t"<<"\tMin risk:\t"
        <<z_population[1]<<endl;
    int generation = 1000;

    for(int i = 0; i<generation; i++){
        //[2] Update population
        //cerr<<"\nGENRATION\t"<<i<<endl;
        population new_x;
        process_updateP(lamblist, z_population, x, new_x, assetArray);
        process_updateZ(new_x, z_population);

        /*
        int test_set[32] = {0,0,191,82,176,0,0,64,0,0,0,69,70,0,226, 0, 76, 0, 0, 58, 0, 0, 0, 0,
                          0, 227, 0, 0, 0, 83, 65};
        solution test_solution;
        for(int i = 0; i<32; i++){
            test_solution.gene.push_back(test_set[i]);
        }
        //util_print_gene(test_solution);
        util_repair_gene(test_solution, assetArray);
         */
        process_updateN(new_x, x, assetArray, lamblist, z_population);
        process_updateEP(x, ep);
        if(i%10==1){
            for(int j = 0; j<ep.xi.size(); j++){
                cerr<<ep.xi[j].fitness[0]<<"\t"<<ep.xi[j].fitness[1]<<"\t"<<i<<"\n";
            }
        }

    }
}