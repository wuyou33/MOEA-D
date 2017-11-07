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

    struct Constraint port1_constraint;
    double port1_correlations[31][31] = {};
    vector <struct asset> assetArray;
    //
    //Check input
    //

    if(loadItem(  FILE_PATH,
                assetArray,
                port1_constraint,
                port1_correlations)) {
        for (int i = 0; i < port1_constraint.num_assets; i++) {
            //cout << assetArray[i].holding << endl;
        }
    }

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
    solution x;
    init_solutions(x, assetArray, port1_constraint, port1_correlations);
}