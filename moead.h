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
#include <gsl>
#include "randG.h"
#include "loader.hpp"
#define N 200

struct lamb{
    //
    //Size of weight vector is determined by the number of objective of MOP.
    //To make it easy to implement, we only consider profit and risk, where k = 2.
    //According to the contribution of [2], the corresponding H = 99, and N = 100;
    //
    double v[2];
};

bool init_lamb( const int &k, const int &H ){

}
double evalue(const std::vector<struct asset> &assetArray){
    double value = 0;
    double income = 0;
    for(int i = 0; i<assetArray.size(); i++){
        value+=assetArray[i].holding*assetArray[i].current_price;
        income += assetArray[i].holding*assetArray[i].mean_income;
    }
    return income/value;
}

void current2target(const std::vector<struct asset> &assetArray,
                    const std::vector<struct asset> &targetArray,
                    std::vector<int> &offset,
                    const Constraint&constraint){
    int buffer;
    for(int i = 0; i<assetArray.size(); i++){
        buffer = targetArray[i].holding - assetArray[i].holding;
        offset.push_back(buffer);
        //For change
    }
}

double surgeOut(const std::vector<struct asset> &assetArray){
    double fundPool = 0;
    for(int i = 0; i < assetArray.size(); i++){
        fundPool += assetArray[i].holding * assetArray[i].current_price;
    }
    return fundPool;
}

void surgeIn(const std::vector<struct asset> &assetArray, const Constraint&constraints, double fundPool){
    if(fundPool<=0)
        return;
    std::vector<int> skeleton;

}



struct solution{
    std::vector<int> gene;
    double fitness[2]={};

    struct solution &operator = (struct solution x){
        this->gene.clear();
        for(int i = 0; i<x.gene.size(); i++){
            this->gene.push_back(x.gene[i]);
        }
        this->fitness[0] = x.fitness[0];
        this->fitness[1] = x.fitness[1];
        return *this;
    }
};
struct Individual{
    std::vector<solution> x;
    int m = 0;//Number of asset & dimansion of solution space

};
size_t expectReturn(const std::vector<int>&x,
                    const std::vector<struct asset> &asset,
                    const size_t solution_size){
    size_t sum = 0;
    for(int i = 0; i<solution_size; i++){
        //sum += (asset[i].current_price*asset[i].holding)+(x[i]*asset[i].current_price);
        sum += (asset[i].mean_income*asset[i].holding);
    }
    return sum;
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

bool dominate(const struct solution &x, const struct solution&y){
    for(int i = 0; i<2; i++){
        if(x.fitness[i]<y.fitness[i])
            return false;
    }
    return true;
}

void solutionInitial(struct solution &x,
                     const std::vector <struct asset> &asset,
                     const double (&correlation)[31][31]){
    x.fitness[0] = expectReturn(x.gene, asset, 31);
    x.fitness[1] = -covariance(x.gene, asset, correlation);
}

void init_all(struct Individual &candidate,
                const std::vector<struct asset> &asset,
                const double (&correlation)[31][31]){
    //N is the number of subproblems.
    //m is the dimension of solution space.
    for(int i = 0; i<N; i++){
        solutionInitial(candidate.x[i], asset, correlation);
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

double EDistance(struct lamb x, struct lamb y){
    double sum = 0;
    sum += (y.v[0] - x.v[0])*(y.v[0] - x.v[0]) + (y.v[1] - x.v[1])*(y.v[1] - x.v[1]);
    return sqrt(sum);
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


double max(const double a, const double b){
    if(a>b)
        return a;
    return b;
}
double Tchebycheff(const struct solution&x,
                   const struct lamb lb,
                   const double z_star[2]){
    double result;
    double temp[2];
    for(int i = 0; i<2; i++){
        temp[i] = std::abs((x.fitness[i]-z_star[i])*lb.v[i]);
    }
    result = max(temp[0], temp[1]);
    return result;
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
    solutionInitial(offspring1, asset, correlation);
    solutionInitial(offspring2, asset, correlation);
    if(dominate(offspring1, offspring2))
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
    struct Individual EP;
    struct Individual S;//Initial population

    for(int i = 0; i<N; i++){
        struct solution solu_buffer;
        int buffer = 0;
        for(int j = 0; j<asset.size();j++){
            int cache_buffer = randomN(asset, j);
            buffer = cache_buffer;
            solu_buffer.gene.push_back(buffer);
        }
        solutionInitial(solu_buffer, asset, correlation);
        S.x.push_back(solu_buffer);
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
                Ed[j][i] = Ed[i][j] = EDistance(lambList[i], lambList[j]);
                c_buffer.index = j;
                c_buffer.distance = Ed[i][j];
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
    for(int i = 0; i<S.x.size(); i++){
        if(S.x[i].fitness[0]>z_star[0]){
            z_star[0] = S.x[i].fitness[0];
        }
        z_star[1] = S.x[i].fitness[1];
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

            struct solution candidate = geneticOperation(S.x[k],
                                                         S.x[l],
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
                double Ta = Tchebycheff(S.x[B[i][j]], lambList[B[i][j]], z_star);
                double Tb = Tchebycheff(candidate, lambList[B[i][j]], z_star);
                if(Tb>Ta){
                    S.x[B[i][j]] = candidate;
                }
            }
            if(EP.x.size() == 0){
                EP.x.push_back(candidate);
            }
            else{
                bool dominateY = false;
                std::vector<int> rmlist;
                for(int j = 0; j<EP.x.size(); j++){
                    if(dominate(candidate, EP.x[j])){
                        rmlist.push_back(j);
                    }
                    else
                    if(dominate(EP.x[j], candidate)){
                        dominateY = true;
                    }
                }
                if(dominateY == false){
                    EP.x.push_back(candidate);
                    for(auto item:rmlist){
                        deleteV(EP.x, item);
                    }
                }
            }
        }
    }
    for(int j = 0; j<N; j++){
        for(int i = 0; i<2; i++){
            if(i == 0)
                std::cout<<EP.x[j].fitness[i]<<"\t";
            else
                std::cout<<EP.x[j].fitness[i]/1000000000000000000<<std::endl;
        }
    }
}


#endif //MOEA_D_MOEAD_H
