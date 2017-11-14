//
//  MOEA-D.cpp
//  MOEA/D
//
//  Created by Molin on 04/08/2017.
//  Copyright © 2017 merlin. All rights reserved.
//
//
#ifndef MOEA_D_MOEAD_H
#define MOEA_D_MOEAD_H
#include <queue>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <random>
#include <exception>
#include "randG.h"
#include "loader.h"
//Population's sizetrue
#define N 100
//K nearest neighbors
#define K 11//Cannot assign 10. I don't know why.
bool mainprocess = false;
bool p_1_3_detail = false;
bool show_gene = false;
bool initial_detail = false;
bool genetic_state = false;
std::vector<asset>raw_asset;
void setting(bool *s, const std::vector<asset>assetList){
    mainprocess = s[0];
    p_1_3_detail = s[1];
    show_gene = s[2];
    initial_detail = s[3];
    for(auto asset_item:assetList){
        raw_asset.push_back(asset_item);
    }
}
struct solution{//i.e. x_i
    std::vector<int> gene;
    std::vector<asset> raw_data;
    double fitness[2]={};
    void init(){
        for(auto asset_item:raw_asset){
            raw_data.push_back(asset_item);
        }
    }
    solution(){
        init();
    }

    struct solution &operator = (const struct solution x){
        this->gene.clear();
        for(int i = 0; i<x.gene.size(); i++){
            this->gene.push_back(x.gene[i]);
        }
        this->fitness[0] = x.fitness[0];
        this->fitness[1] = x.fitness[1];
        return *this;
        for(int i = 0; i<x.raw_data.size(); i++){
            this->raw_data[i] = x.raw_data[i];
        }
    }
};
struct population{
    std::vector<solution> xi;
    struct population &operator = (const struct population x){
        for(int i = 0; i<x.xi.size(); i++){
            solution temp_push;
            temp_push = x.xi[i];
            this->xi.push_back(temp_push);
        }
        return *this;
    }
};
struct lamb{
    //
    //Size of weight vector is determined by the number of objective of MOP.
    //To make it easy to implement, we only consider profit and risk, where k = 2.
    //According to the contribution of [2], the corresponding
    // H = 99, and N = 100;
    //
    double v[2] = {};
    int id = 0;
    std::vector<lamb>k_nearest;
    lamb &operator = (const lamb&y){
        for(int i = 0; i<2; i++){
            this->v[i] = y.v[i];
            i++;
        }
        this->v[0] = y.v[0];
        this->v[1] = y.v[1];
        for(int i = 0; i<this->k_nearest.size(); i++){
            this->k_nearest[i] = y.k_nearest[i];
        }
        this->id = y.id;
        return *this;
    }

};
bool operator<(lamb a, lamb b){
    return a.v[0] < b.v[0];
}

void util_evalue( solution &xi, std::vector<asset>assetList){
    double income = 0;
    double risk = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            income += xi.gene[i] * assetList[i].mean_income;
            risk += xi.gene[i] * assetList[i].diviation_r;
        }
    }
    xi.fitness[0] = income;
    xi.fitness[1] = risk;
}
void util_print_gene( const solution&xi){
    int counter = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i] != 0) counter++;
        std::cerr<<xi.gene[i]<<"\t";
    }
    std::cerr<<"Number: "<<counter<<"\tFitness: "<<xi.fitness[0]<<" "<<xi.fitness[1]<<std::endl;
    /*
    for(int i = 0; i<xi.raw_data.size(); i++){
        std::cerr<<xi.raw_data[i].id<<"\t";
    }
    std::cerr<<std::endl;
     */
    for(int i = 0; i<xi.raw_data.size(); i++){
        std::cerr<<xi.raw_data[i].max_buy<<"\t";
    }
    std::cerr<<std::endl;
}
bool util_dominate(const struct solution &x, const struct solution&y){
    if(x.fitness[0]>y.fitness[0]&&x.fitness[1]<y.fitness[1])
        return true;
    return false;
}
void util_repair_gene( solution &xi, const std::vector<asset>&assetList ){
    //1. Check Cardinality Constraint
    //2. Check Fund Constraint
    //3. Check min&max Buy Limit Constraint
    //
    //Check cardinality constraint
    //
    int M = 10;
    bool cac = false;
    int counter = 0;
    std::vector<int>port;

    //Caculate Fund Constraint
    int fund = 0;
    for(int i = 0; i<assetList.size(); i++){
        fund += assetList[i].current_price*assetList[i].holding;
    }

    //Caculate Solution's Fund
    int local_fund = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            counter++;
            port.push_back(i);
        }
    }
    solution local_best;
    if(counter>M){
        int differ = counter - M;
        for(int i = 0; i<differ; i++){
            solution local_best_2;
            for(int j = 0; j<port.size(); j++){
                solution xi_buffer;
                if(i == 0)
                    xi_buffer= xi;
                else
                    xi_buffer = local_best;
                if(xi_buffer.gene[port[j]]!=0)
                    xi_buffer.gene[port[j]] = 0;
                else
                    continue;
                //std::cerr<<"\nx"<<j<<" buffer:\n";
                util_evalue(xi_buffer, assetList);
                //util_print_gene(xi_buffer);
                if(j==0){
                    local_best_2 = xi_buffer;
                    //util_print_gene(local_best_2);
                }else{
                    if(util_dominate(xi_buffer, local_best_2)){
                        local_best_2 = xi_buffer;
                    }
                }
            }
            local_best = local_best_2;
            /*
            std::cerr<<"\nLocal best 2\n";
            util_print_gene(local_best_2);
            std::cerr<<"\nLocal best:\n";
            util_print_gene(local_best);
            std::cerr<<std::endl;
             */
        }
        xi = local_best;
    }

    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            if(xi.gene[i]>assetList[i].max_buy){
                xi.gene[i] = assetList[i].max_buy;
            }
            if(xi.gene[i]<assetList[i].min_buy){
                xi.gene[i] = assetList[i].min_buy;
            }
        }
    }

    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            local_fund+=xi.gene[i]*assetList[i].current_price;
        }
    }
    if(local_fund>fund){
        int differ = local_fund-fund;
        for(int i = 0; i<xi.gene.size(); i++){
            int minus = xi.gene[i]/local_fund*differ;
            if(xi.gene[i] - minus >=assetList[i].min_buy)
                xi.gene[i]-=minus;
        }
    }
    util_evalue(xi, assetList);
}
void util_findK( const int i, int k, int upper, int*h){
    int u = 0;
    int d[2] = {};
    bool up = false;
    if(k%2!=0)
        up = true;
    if(i+k/2+(up?1:0)>=upper){
        d[1] = upper-(i+k/2)+(up?1:0);
        d[0] = k-d[1];
    }
    else{
        d[1] = k/2+(up?1:0);
        if(i-k/2-(up?1:0)>0)
            d[0] = k/2+(up?1:0);
        else{
            d[0] = i;
            d[1] = k-i;
        }
    }
    for (int m = 1; m<=d[1]; m++) {
        int buffer = i+m;
        h[u] = buffer;
        u++;
    }
    for(int m = 1; m<=d[0]; m++){
        if(i == 0 )
            std::cerr<<"error"<<std::endl;
        int buffer = i-m;
        h[u] = buffer;
        u++;
    }
}
bool util_isfeasible( const solution &xi, const std::vector<asset>&assetList){
    int fund = 0, local_fund = 0, counter = 0;
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0)
            counter++;
        fund+=assetList[i].current_price*assetList[i].holding;
        local_fund +=xi.gene[i]*assetList[i].current_price;
        if(xi.gene[i]>assetList[i].max_buy||xi.gene[i]<assetList[i].min_buy)
            return false;
    }
    if(counter>10)  return false;
    if(local_fund>fund) return false;
    return true;
}
void util_display_nn(const std::vector<lamb> &lamblist){
    for(int i = 0; i<lamblist.size(); i++) {
        std::cout << lamblist[i].v[0] << ", " << lamblist[i].v[1] << ":\t\t";
        for (int j = 0; j < lamblist[i].k_nearest.size(); j++) {
            std::cout << lamblist[i].k_nearest[j].v[0] << "," <<
                      lamblist[i].k_nearest[j].v[1] << "\t";
        }
        std::cout << std::endl;
    }
}
void util_display_nn_id( const std::vector<lamb> &lamblist){
    for(int i = 0; i<lamblist.size(); i++){
        std::cout<<lamblist[i].id<<":\t\t";
        for(int j = 0; j<lamblist[i].k_nearest.size(); j++){
            std::cout<<lamblist[i].k_nearest[j].id<<",\t";
        }
        std::cout<<std::endl;
    }
}
void util_fundDistribute( std::vector<struct asset>&p, double&fund, int n, int mode){

    if(p[n].max_buy<p[n].min_buy){
        check:
        bool buyable = false;
        for (int i = 0; i < p.size(); i++) {
            if (p[i].max_buy > p[i].min_buy&&fund/p[i].current_price>p[i].min_buy) {
                buyable = true;
            }
        }
        if(buyable){
            for (int i = 0; i < p.size(); i++) {
                if (p[i].max_buy > p[i].min_buy&&fund/p[i].current_price>p[i].min_buy) {
                    util_fundDistribute(p, fund, i, 0);
                }
            }
        }
        if(buyable) goto check;
        return;
    }

    if(fund/p[n].current_price<p[n].min_buy){
        for(int i = 0; i<p.size(); i++){
            if(fund/p[i].current_price>=p[i].min_buy){
                util_fundDistribute(p, fund, i, 0);
            }
        }
        return;
    }

    double dice = randG();
    double raw_buy = (dice * (p[n].max_buy-p[n].min_buy));
    if(raw_buy<0){
        std::cerr<<"!!!!Error\n";
    }
    int buy_number = static_cast<int>(dice * (p[n].max_buy-p[n].min_buy) + p[n].min_buy);

    double local_fund = buy_number*p[n].current_price;
    if(local_fund>fund){
        local_fund = fund;
        buy_number = local_fund/p[n].current_price;
        if(buy_number == 0){
            //std::cerr<<local_fund<<std::endl;
        }
        return;
    }
    p[n].fundpool+=local_fund;
    p[n].buy_asset_number += buy_number;
    p[n].max_buy -=buy_number;
    p[n].history.push_back(-buy_number);
    fund -=local_fund;
    if(n<=0||mode == 1){
        rn:
        double s = randG();
        util_fundDistribute(p, fund, int(s*p.size()),1);
    }
    else {
        util_fundDistribute(p, fund, n - 1, 0);
    }
}
void util_preprocess(std::vector<struct asset> &assetArray){
    for(int i = 0; i<assetArray.size(); i++){
        assetArray[i].max_buy = assetArray[i].max_buy+assetArray[i].holding;
        //assetArray[i].holding = 0;
    }
}
void util_ffunction(const solution &xi,
                    const std::vector<asset> &assetArray,
                    double &income,
                    double &risk){
    for(int i = 0; i<xi.gene.size(); i++){
        if(xi.gene[i]!=0){
            income+=xi.gene[i]*assetArray[i].mean_income;
            risk+=xi.gene[i]*assetArray[i].diviation_r;
        }
    }
    if(p_1_3_detail) {
        std::cerr << "[1.3]: Income:\t" << income
                  << "\tRisk:\t" << risk << std::endl;
    }
}
double max(const double a, const double b){
    if(a>b)
        return a;
    return b;
}
double util_tchebycheff(const struct solution &x,
                        const struct lamb lb,
                        const double *z_star){
    double result;
    double temp[2];
    for(int i = 0; i<2; i++){
        temp[i] = std::abs((x.fitness[i]-z_star[i])*lb.v[i]);
    }
    result = max(temp[0], temp[1]);
    return result;
}

void process_updateZ( const population &x, double z[]){
    for(int i = 0; i<x.xi.size(); i++){
        if(x.xi[i].fitness[0]>z[0]){
            z[0] = x.xi[i].fitness[0];
        }
        if(x.xi[i].fitness[1]<z[1]){
            z[1] = x.xi[i].fitness[1];
        }
    }
    if(mainprocess){
        std::cerr<<"[2.3]:\tUpdate z:\t"
            <<"\n\t\tMax income:\t"<<z[0]<<"\t"<<"\tMin risk:\t"
            <<z[1]<<std::endl;
    }
}
void process_genetic( solution &m, solution &n, const std::vector<asset>&assetList, solution&trail){
    //
    //Crossover - one point
    //
    int M = 10;
    solution trail_y;
    int dot = randG()*m.gene.size();
    for(int i = 0; i<m.gene.size(); i++) {
        int trail_y_buffer = 0;
        if (i < dot) {
            trail_y_buffer = m.gene[i];
            trail_y.raw_data[i] = m.raw_data[i];
        }
        else{
            trail_y_buffer = n.gene[i];
            trail_y.raw_data[i] = n.raw_data[i];
        }
        trail_y.gene.push_back(trail_y_buffer);
    }
    if(genetic_state) {
        std::cerr << "M:\n";
        util_print_gene(m);
        std::cerr << "N:\n";
        util_print_gene(n);
        std::cerr << "Trail before repair\n";
        util_print_gene(trail_y);
        util_repair_gene(trail_y, assetList);
        std::cerr << "Trail after repair\n";
        util_print_gene(trail_y);
    }
    double roll_dice = randG();
    double mu = 1;
    if(roll_dice<mu){
        //
        //Mutation
        //
        int counter = 0;
        std::vector<asset>trail_asset;
        std::vector<int>random_list;
        double trail_fund = 0;
        double fund = 0;
        for(int i = 0; i<trail_y.gene.size(); i++){
            fund += assetList[i].current_price*assetList[i].holding;
            if(trail_y.gene[i]!=0){
                counter++;
                trail_fund += trail_y.gene[i] * assetList[i].current_price;
            }
            else{
                random_list.push_back(i);
            }
        }
        //std::cerr<<"Counter:\t"<<counter<<"\tBuffer fund:\t"<<trail_fund<<std::endl;
        std::vector<int>index_random;
        for(int i = 0; i<M-counter; i++){
            int temp_index;
            roll:
            temp_index = randG()*random_list.size();
            for(auto item:index_random){
                if(temp_index == item)
                    goto roll;
            }
            index_random.push_back(temp_index);
            asset temp_trail_meta = trail_y.raw_data[random_list[temp_index]];
            //std::cerr<<"Empty index:\t"<<random_list[temp_index]<<"\t";
            trail_asset.push_back(temp_trail_meta);
        }
        if(counter<M&&trail_fund<fund){
            int differ = counter-M;
            double buffer_fund = fund - trail_fund;
            //std::cerr<<"Size of trail asset\t"<<trail_asset.size()<<std::endl;
            util_fundDistribute(trail_asset, buffer_fund, trail_asset.size()-1, 0);
        }
        /*
        for(int i = 0; i<trail_asset.size(); i++){
            std::cerr<<trail_asset[i].buy_asset_number<<"\t";
        }
         */
        for(int i = 0; i<trail_y.gene.size(); i++){
            for(int j = 0; j<trail_asset.size(); j++){
                if(i == trail_asset[j].id){
                    if(trail_y.gene[i] != 0)
                        std::cerr<<"Error\n";
                    trail_y.gene[i] = trail_asset[j].buy_asset_number;
                    trail_y.raw_data[i] = trail_asset[j];
                }
            }
        }
    }
    util_evalue(trail_y, assetList);
    //util_print_gene(trail_y);
    trail = trail_y;
}
void process_updateP(const std::vector<lamb>&lamblist,
             double z[2],
             const population &x,
             population &new_x,
             const std::vector<asset>&assetList ){
    if(mainprocess){
        std::cerr<<"[2]:\tUpdate\n";
    }
    //
    //  [2.1] Reproduction
    //
    if(mainprocess){
        std::cerr<<"[2.1]:\tStart\tReporduction\n";
    }
    for(int i = 0; i<N; i++){
        solution trail_buffer;
        int m = lamblist[i].k_nearest[static_cast<int>(randG() * lamblist[i].k_nearest.size())].id;
        int n = lamblist[i].k_nearest[static_cast<int>(randG() * lamblist[i].k_nearest.size())].id;

        solution m_buffer = x.xi[m];
        solution n_buffer = x.xi[n];
        //if(mainprocess) std::cerr<<"[2.1]:\tm:"<<m<<"\tn:"<<n<<std::endl;
        process_genetic(m_buffer, n_buffer, assetList, trail_buffer);
        new_x.xi.push_back(trail_buffer);
    }
}
void process_updateN( population&new_x,
                      population&x,
                      const std::vector<asset>&assetList,
                      const std::vector<lamb>&lambList,
                      const double z[]){
    if(mainprocess){
        std::cerr<<"[2.4]: Start\tUpdate Neighborhood Population"<<std::endl;
    }
    for(int i = 0; i<N; i++){
        for(int j = 0; j<lambList[i].k_nearest.size(); j++){
            double y = util_tchebycheff(new_x.xi[i], lambList[i].k_nearest[j], z);
            double x_j = util_tchebycheff(x.xi[lambList[i].k_nearest[j].id], lambList[i].k_nearest[j], z);
            if (y<=x_j){
                x.xi[lambList[i].k_nearest[j].id] = new_x.xi[i];
            }
        }
    }
    if(mainprocess)   std::cerr<<"[1.3]: *End*\tUpdate Neighborhood Population"<<std::endl;
}

void process_updateEP( const population&x, population&ep){
    if(mainprocess){
        std::cerr<<"[2.5]: Start\tUpdate EP\n";
    }
    for(int i = 0; i<x.xi.size(); i++){
        for(int j = 0; j<x.xi.size(); j++){
            if(util_dominate(x.xi[i], ep.xi[j])){
                ep.xi[j] = x.xi[i];
            }
        }
    }
    if(mainprocess){
        std::cerr<<"[2.5]: *End*\tUpdate EP\n";
    }
}
bool init_lamb( const int &k, const int &H, std::vector<lamb>&lamb_list ){
    for(int i = 0; i<N; i++){
        double source_number = randG();
        struct lamb buffer;
        buffer.id = i;
        buffer.v[0] = static_cast<double>(static_cast<int>(source_number*H))/ static_cast<double>(H);
        //std::cerr<<i<<":"<<buffer.v[0]<<std::endl;
        buffer.v[1] = 1 - buffer.v[0];
        lamb_list.push_back(buffer);
    }
    return true;
}
void init_distance( std::vector<lamb> &lamb_list ){
    if(mainprocess)   std::cerr<<"[1.2]: Start\tInit Distance"<<std::endl;
    std::priority_queue<lamb>nearest_list;
    int length = lamb_list.size();
    for(int i = 0; i<length; i++){
        lamb buffer = lamb_list.back();
        //std::cout<<buffer.v[0]<<"\t"<<buffer.v[1]<<std::endl;
        lamb_list.pop_back();
        nearest_list.push(buffer);
    }
    std::vector<lamb> toHandle;
    while(!nearest_list.empty()){
        lamb buffer = nearest_list.top();
        toHandle.push_back(buffer);
        nearest_list.pop();
    }

    for(int i = 0; i<toHandle.size(); i++){
        lamb buffer = toHandle[i];
        int h[K] = {};
        util_findK(i, K, toHandle.size(), h);
        for(int j = 0; j<K; j++){
            buffer.k_nearest.push_back(toHandle[h[j]]);
        }
        lamb_list.push_back(buffer);
    }
    if(mainprocess)   std::cerr<<"[1.2]: *End*\tInit Distance"<<std::endl;
}


void init_solution( std::vector<struct asset> &tobuy,
                    double&fund){
    util_fundDistribute(tobuy, fund, tobuy.size()-1, 0);
    if(initial_detail){
        for(int i = 0; i<tobuy.size(); i++){
            std::cerr<<"ID: "<<tobuy[i].id<<"\tMax: "<<tobuy[i].max_buy<<"\tMin: "
                     <<tobuy[i].min_buy<<"\tBuy: "<<tobuy[i].buy_asset_number<<"\t";
            for(auto item:tobuy[i].history){
                std::cerr<<item<<"\t";
            }
            std::cerr<<std::endl;
        }
    }
}
void init_population(struct population &x,
                     const std::vector <struct asset> &asset,
                     const struct Constraint constraint,
                     const double (&correlation)[31][31]){
    if(mainprocess)   std::cerr<<"[1.3]: Start\tInit Population"<<std::endl;
    //Cardinality Constraint M:
    int M = constraint.max_assets;
    int num_assets = constraint.num_assets;
    //std::cerr<<"Number of Assets:\t"<<num_assets<<std::endl;
    //
    double fundpool = 0;
    size_t fundpool_l = 0;
    for(int i = 0; i<asset.size(); i++) {
        fundpool += asset[i].current_price * asset[i].holding;
        fundpool_l += asset[i].current_price * asset[i].holding;
    }
    for(int i = 0; i<N; i++){
        solution xi_buffer;
        struct Set candidate;
        //Handling cardinality constraint by using pointers.
        //std::cerr<<"[1.3]: \t\tStart\tassign assets\n";
        std::vector<struct asset>tobuy;
        for(int i = 0; i<M; i++){
            int buffer = static_cast<int>(randG()*(num_assets+1));
            if(buffer>=num_assets)
                buffer = -1;//Hold cash
            //std::cerr<<buffer<<"\t";
            if(buffer>=0) {
                if(!candidate.isin(buffer)){
                    candidate.data.push_back(buffer);
                }
            }
        }
        for(int i = 0; i<candidate.data.size(); i++){
            int buffer_id = candidate.data[i];
            struct asset asset_buffer = asset[buffer_id];
            asset_buffer.id = buffer_id;
            tobuy.push_back(asset_buffer);
        }
        std::priority_queue<int, std::vector<int>, std::greater<int>>tobuy_ase;
        for(int i = 0; i<tobuy.size(); i++){
            //std::cerr<<tobuy[i].id<<"\t";
            int buffer = tobuy[i].id;
            tobuy_ase.push(buffer);
        }
        if(p_1_3_detail) {
            std::cerr << "[1.3]:\tPort " << i << ":\t";
            while (!tobuy_ase.empty()) {
                std::cerr << tobuy_ase.top() << "\t";
                tobuy_ase.pop();
            }
            std::cerr << std::endl;
        }
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tEnd\t\tassign assets\n";
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tStart\tInitialize Solution\n";
        double fund_buffer = fundpool;
        init_solution(tobuy, fund_buffer);
        if(p_1_3_detail)   std::cerr<<"[1.3]: Fund remain:"<<fund_buffer<<std::endl;
        //if(mainprocess)   std::cerr<<"[1.3]: \t\tEnd\t\tInitialize Solution\n";
        for(int j = 0; j<asset.size(); j++){
            xi_buffer.gene.push_back(0);
        }
        for(auto item:tobuy){
            xi_buffer.gene[item.id] = item.buy_asset_number;
            xi_buffer.raw_data[item.id] = item;
        }
        double income_buffer = 0;
        double risk_buffer;
        util_ffunction(xi_buffer, asset, income_buffer, risk_buffer);
        xi_buffer.fitness[0] = income_buffer;
        xi_buffer.fitness[1] = risk_buffer;
        x.xi.push_back(xi_buffer);
        if(show_gene) {
            for (int i = 0; i < xi_buffer.gene.size(); i++) {
                std::cout << xi_buffer.gene[i] << " ";
            }
            std::cout << std::endl;
        }
    }
    if(mainprocess)   std::cerr<<"[1.3]: *End*\tInit Population"<<std::endl;
}


void test_generator( std::vector<asset>&assetList ){

}



void de_process(std::vector<int>h, std::vector<int>m, std::vector<int>n){
    solution trail;
    //
    //Mutation
    //
    double F = 0.5;
    double CR = 0.4;
    for(int i = 0; i<h.size(); i++){

        double dice = randG();
        if(dice<CR||i==h.size()-1) {
            int trail_buffer = static_cast<int>(h[i] + F * (m[i] - n[i]));
            trail.gene.push_back(trail_buffer);
        }
        else{

        }
    }
    //
    //Crossover
    //
    for(int i = 0; i<h.size(); i++){

    }
}
size_t covariance(const std::vector<int>&x,
                  const std::vector<struct asset> &asset,
                  const double (&correlation)[31][31]){
    size_t cov = 0;
    for(int i = 0; i<31; i++){
        for(int j = 0; j<31; j++){
            cov +=asset[i].current_price*(x[i]+asset[i].holding)*asset[j].current_price*(x[j]+asset[j].holding)*correlation[i][j];
        }
    }
    return cov;
}




void init_all(struct population &candidate,
                const std::vector<struct asset> &asset,
                const double (&correlation)[31][31]){
    //N is the number of subproblems.
    //m is the dimension of solution space.
    for(int i = 0; i<N; i++){
        //init_solutions(candidate.x[i], asset, correlation);
    }
}
int calNumber(const std::vector<asset>&assetArray, const std::vector<int>&gene){
    int counter = 0;
    for(int i = 0; i< assetArray.size();i++){
        if(assetArray[i].holding + gene[i] !=0)
            counter++;
    }
    return counter;
}

int calNumber2(const std::vector<asset>&assetArray){
    int counter = 0;
    for(int i = 0; i < assetArray.size(); i++){
        if(assetArray[i].holding !=0)
            counter++;
    }
    return counter;
}
void randomN2(const std::vector<asset>&assetArray,
              const Constraint&constraints,
              std::vector<int>&output){
    std::cout<<"Origin:"<<calNumber2(assetArray)<<std::endl;
    re:
    for(int i = 0; i<assetArray.size(); i++){
        int max_sell = assetArray[i].max_sell;
        int max_buy = assetArray[i].max_buy;
        int min_sell = assetArray[i].min_sell;
        int min_buy = assetArray[i].min_buy;
        int holding = assetArray[i].holding;

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(0, 10);
        int isbuy = static_cast<int>((dist(mt)))%2;
        if(isbuy == 0){
            int len = max_sell - min_sell + 1;
            std::uniform_real_distribution<double> dist(0, 100000);
            int divide = 1;
            if(holding > max_sell){
                divide = len;
            }
            else{
                if(holding != 0)
                    divide = holding - min_sell + 1;
            }
            int result = static_cast<int>(dist(mt))%divide + 1;
            if(holding <=0){
                output.push_back(0);
            }
            else{
                output.push_back(-result);
            }
        }

        if(isbuy == 1){
            int len = max_buy - min_buy+1;
            std::uniform_real_distribution<double> dist(0, 100000);
            int result = static_cast<int>(dist(mt))%len+1;

            //Log
            /*
             std::cout<<"max sell:"<<max_sell<<std::endl<<
             "min sell:"<<min_sell<<std::endl<<
             "max buy:"<<max_buy<<std::endl<<
             "min buy:"<<min_buy<<std::endl;

             std::cout<<"length:"<<len<<std::endl<<
             "result:"<<result<<std::endl;
             */
            output.push_back(result);
        }

    }
    std::cout<<"Number"<<calNumber(assetArray, output)<<std::endl;
    /*
    if(calNumber(assetArray, output)>constraints.max_assets){
        output.clear();
        goto re;
    }
     */
}
int randomN(const std::vector <asset>&assetArray, int index){
    /*
    int len;
    len = max_sell - min_sell + 1 + max_buy - min_buy + 1;
    int l_len = max_sell - min_sell + 1;
    int r_len = max_buy - min_buy + 1;
     */
    int max_sell = assetArray[index].max_sell;
    int max_buy = assetArray[index].max_buy;
    int min_sell = assetArray[index].min_sell;
    int min_buy = assetArray[index].min_buy;
    int holding = assetArray[index].holding;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0, 10);
    int isbuy = static_cast<int>((dist(mt)))%2;
    if(isbuy == 0){
        int len = max_sell - min_sell + 1;
        if(len>0){
            if(holding > max_sell){
                std::uniform_real_distribution<double> dist(0, 100000);
                int result = static_cast<int>(dist(mt))%len+1;

                //Log
                /*
                std::cout<<"max sell:"<<max_sell<<std::endl<<
                "min sell:"<<min_sell<<std::endl<<
                "max buy:"<<max_buy<<std::endl<<
                "min buy:"<<min_buy<<std::endl;

                std::cout<<"length:"<<len<<std::endl<<
                "result:"<<result<<std::endl;
                 */
                return -result;
            }
            else{
                if(holding!=0){
                    std::uniform_real_distribution<double> dist(0, 100000);
                    int result = static_cast<int>(dist(mt))%holding+1;

                    //Log
                    /*
                    std::cout<<"max sell:"<<max_sell<<std::endl<<
                    "min sell:"<<min_sell<<std::endl<<
                    "max buy:"<<max_buy<<std::endl<<
                    "min buy:"<<min_buy<<std::endl;

                    std::cout<<"length:"<<len<<std::endl<<
                    "result:"<<result<<std::endl;
                     */
                    return -result;
                }
                else
                    return 0;
            }
        }
        return 0;
    }

    if(isbuy == 1){
        int len = max_buy - min_buy+1;
        std::uniform_real_distribution<double> dist(0, 100000);
        int result = static_cast<int>(dist(mt))%len+1;

        //Log
        /*
        std::cout<<"max sell:"<<max_sell<<std::endl<<
        "min sell:"<<min_sell<<std::endl<<
        "max buy:"<<max_buy<<std::endl<<
        "min buy:"<<min_buy<<std::endl;

        std::cout<<"length:"<<len<<std::endl<<
        "result:"<<result<<std::endl;
         */
        return result;
    }

    return 0;
}

int randomIndex(int x){
    int result;
    srand((unsigned)time(NULL));

    result = (rand()%x);
    return result;
}

double randomDouble(){
    double x = ((double) rand() / (RAND_MAX));
    return x;
}


struct closet{
    double distance = 0;
    int index = 0;

    friend bool operator < (closet x, closet y){
        return (x.distance < y.distance);
    }
};
struct solution crossover(const std::vector<int>&x, const std::vector<int> &y){
    struct solution offspring;
    int point = randomIndex(x.size());
    for(int i = 0; i<point; i++){
        offspring.gene.push_back(x[i]);
    }
    for(int i = point; i<y.size(); i++){
        offspring.gene.push_back(y[i]);
    }
    return offspring;
}
void mutation(struct solution&origin,
              const std::vector<asset> &asset){
    int i = randomIndex(31);
    origin.gene[i] = randomN(asset, i);
}




struct solution geneticOperation(const struct solution&x,
                                 const struct solution&y,
                                 double rate,
                                 double mutate_rate,
                                 const std::vector<asset> &asset,
                                 const double (&correlation)[31][31]){
    struct solution offspring1, offspring2;
    double p = randomDouble();
    if(p<rate){
        offspring1 = crossover(x.gene, y.gene);
        offspring2 = crossover(y.gene, x.gene);
    }
    else{
        offspring1 = x;
        offspring2 = y;
    }
    double p2 = randomDouble();
    if(p2<mutate_rate){
        mutation(offspring1, asset);
        mutation(offspring2, asset);
    }
    //init_solutions(offspring1, asset, correlation);
    //init_solutions(offspring2, asset, correlation);
    if(util_dominate(offspring1, offspring2))
        return offspring1;
    return offspring2;
}

void deleteV(std::vector<struct solution> x, int index){
    if(index == x.size()-1){
        x.pop_back();
        return;
    }
    std::vector<struct solution>temp;

    for(int i = index+1; i<x.size(); i++){
        temp.push_back(x[i]);
    }
    for(int i = 0; i<x.size()-index; i++){
        x.pop_back();
    }
    for(auto item:temp){
        x.push_back(item);
    }
}

//create a fund pool as constraint

void mainProcess(const std::vector<struct asset>&asset,
                 const double (&correlation)[31][31],
                 const int T,
                 double rate,
                 double mutate_rate,
                 const Constraint&constraints){
    //Read data from file
    struct population EP;
    struct population S;//Initial population

    for(int i = 0; i<N; i++){
        struct solution solu_buffer;
        int buffer = 0;
        for(int j = 0; j<asset.size();j++){
            int cache_buffer = randomN(asset, j);
            buffer = cache_buffer;
            solu_buffer.gene.push_back(buffer);
        }
        //init_solutions(solu_buffer, asset, correlation);
        S.xi.push_back(solu_buffer);
    }
    //Generate lamb
    std::vector<struct lamb> lambList;
    for(int i = 1; i<=N; i++){
        struct lamb lamb_cache;
        lamb_cache.v[0] = static_cast<double>(i)/static_cast<double>(100);
        lamb_cache.v[1] = 1-lamb_cache.v[0];
        lambList.push_back(lamb_cache);
    }

    //Euclidean distances
    int B[N][T];
    memset(B, 0, sizeof(B));
    int Ed[N][N] = {};
    for(int i = 0; i<N; i++){
        std::priority_queue<closet> B_buffer;
        for(int j = 0; j<N; j++){
            if(Ed[i][j] == 0){
                struct closet c_buffer;
                //Ed[j][i] = Ed[i][j] = EDistance(lambList[i], lambList[j]);
                //c_buffer.index = j;
                //c_buffer.distance = Ed[i][j];
                B_buffer.push(c_buffer);
            }
        }
        for(int k = 0; k < T; k++){
            B[i][k] = B_buffer.top().index;
            B_buffer.pop();
        }
        B_buffer.empty();
    }

    double z_star[2] = {};
    for(int i = 0; i<S.xi.size(); i++){
        if(S.xi[i].fitness[0]>z_star[0]){
            z_star[0] = S.xi[i].fitness[0];
        }
        z_star[1] = S.xi[i].fitness[1];
    }

    int t = 0;
    while(t<20){
        //Update
        t++;
        for(int i = 0; i<N; i++){
            int k_b = randomIndex(T);
            int l_b = randomIndex(T);
            int k = B[i][k_b];
            int l = B[i][l_b];

            struct solution candidate = geneticOperation(S.xi[k],
                                                         S.xi[l],
                                                         rate,
                                                         mutate_rate,
                                                         asset,
                                                         correlation);
            //Update z_star
            for(int j = 0; j<2; j++){
                if(candidate.fitness[j]>z_star[j]){
                    z_star[j] = candidate.fitness[j];
                }
            }
            for(int j = 0; j<T; j++){
                double Ta = util_tchebycheff(S.xi[B[i][j]], lambList[B[i][j]], z_star);
                double Tb = util_tchebycheff(candidate, lambList[B[i][j]], z_star);
                if(Tb>Ta){
                    S.xi[B[i][j]] = candidate;
                }
            }
            if(EP.xi.size() == 0){
                EP.xi.push_back(candidate);
            }
            else{
                bool dominateY = false;
                std::vector<int> rmlist;
                for(int j = 0; j<EP.xi.size(); j++){
                    if(util_dominate(candidate, EP.xi[j])){
                        rmlist.push_back(j);
                    }
                    else
                    if(util_dominate(EP.xi[j], candidate)){
                        dominateY = true;
                    }
                }
                if(dominateY == false){
                    EP.xi.push_back(candidate);
                    for(auto item:rmlist){
                        deleteV(EP.xi, item);
                    }
                }
            }
        }
    }
    for(int j = 0; j<N; j++){
        for(int i = 0; i<2; i++){
            if(i == 0)
                std::cout<<EP.xi[j].fitness[i]<<"\t";
            else
                std::cout<<EP.xi[j].fitness[i]/1000000000000000000<<std::endl;
        }
    }
}


#endif //MOEA_D_MOEAD_H
