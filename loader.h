//
//  loader.hpp
//  MOEA
//
//  Created by Molin on 04/11/2017.
//  Copyright © 2017 merlin. All rights reserved.
//
#ifndef MOEA_D_LOADER_H
#define MOEA_D_LOADER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstring>

//
//Define structs here
//
struct Constraint
{
    unsigned int num_assets = 0;
    unsigned int max_assets = 1;
    size_t transaction_limit = 0;
    size_t cash_change = 0;
};

struct asset
{
    double current_price = 0;
    double holding = 0;
    double cost_buy = 0;
    double cost_sell = 0;
    double vcost_buy = 0;
    double vcost_sell = 0;
    double min_buy = 0;
    double min_sell = 0;
    double max_buy = 0;
    double max_sell = 0;

    double mean_income = 0;
    double diviation_r = 0;

    bool sellable = false;

    struct asset &operator=(struct asset x)
    {
        this->current_price = x.current_price;
        this->holding = x.holding;
        this->cost_buy = x.cost_buy;
        this->cost_sell = x.cost_sell;
        this->vcost_buy = x.vcost_buy;
        this->vcost_sell = x.vcost_sell;
        this->min_buy = x.min_buy;
        this->min_sell = x.min_sell;
        this->max_buy = x.max_buy;
        this->max_sell = x.max_sell;
        return *this;
    }
};

//
//Loader function
//
bool loadItem( std::string path,
                std::vector<struct asset> &assetArray,
                struct Constraint &constraint,
                double (&correlation)[31][31])
{
    //
    //Initialize openning file
    //
    //Line number starts from 1
    //
    std::fstream data(path, std::ios::in);
    std::string input_coach;
    unsigned int line = 1;
    if(!data.is_open()){
        std::cerr<<"Error: "<<strerror(errno)<<std::endl;
        return false;
    }

    //
    //Read from file
    //
    //Part 2
    //Line format:
    //${current_price} ${current holding} ${buying cost} ${selling cost} ${variable buying cost} ${variable selling cost} ${min_buy} ${min_sell} ${max_buy} ${max_sell}
    //
    while( std::getline(data, input_coach)){
        std::istringstream line_coach(input_coach);
        if(line == 1){
            line_coach>>
                      constraint.num_assets>>
                      constraint.max_assets>>
                      constraint.transaction_limit>>
                      constraint.cash_change;

            line++;
            continue;
        }
        if(line<=constraint.num_assets+1){
            struct asset assetInput;
            line_coach>>assetInput.current_price>>
                      assetInput.holding>>
                      assetInput.cost_buy>>
                      assetInput.cost_sell>>
                      assetInput.vcost_buy>>
                      assetInput.vcost_sell>>
                      assetInput.min_buy>>
                      assetInput.min_sell>>
                      assetInput.max_buy>>
                      assetInput.max_sell;
            if(assetInput.max_sell != 0)
                assetInput.sellable = true;
            assetArray.push_back(assetInput);
        }
        //
        //Part 3
        //Asset properties
        //${mean_return} ${standard deviation of return}
        if(line>constraint.num_assets+1&&line<=(2*constraint.num_assets+1)){
            int tran = line - constraint.num_assets - 2;
            line_coach>>assetArray[tran].mean_income>>
                      assetArray[tran].diviation_r;

        }
        //
        //Part 4
        //Correlation
        if(line>(2*constraint.num_assets +1)){
            int i,j = 0;
            double corr_buffer = 0;

            line_coach>>i>>j>>corr_buffer;
            correlation[i-1][j-1] = corr_buffer;

        }
        line++;
    }
    return true;
}

#endif
