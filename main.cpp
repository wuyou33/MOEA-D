//
//  instance.cpp
//  MOEA/D
//
//  Created by Molin on 04/08/2017.
//  Copyright © 2017 merlin. All rights reserved.
//
//
#include "loader.hpp"
#include <iostream>
using namespace std;

int main(){
    std::string FILE_PATH = "/home/molin/CLionProjects/MOEA-D/DataSet/portreb1.txt";

    struct Constraint port1_constraint;
    double port1_correlations[31][31] = {};
    vector <struct asset> assetArray;

    if(loadItem(  FILE_PATH,
                assetArray,
                port1_constraint,
                port1_correlations)) {
        for (int i = 0; i < port1_constraint.num_assets; i++) {
            cout << assetArray[i].holding << endl;
        }
    }
}