//
//  main.cpp
//  AonoEnc
//
//  Created by Mayank Rathee on 21/08/17.
//  Copyright © 2017 Mayank Rathee. All rights reserved.
//

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <pari/pari.h>
#include <time.h>
#include <string.h>
#define PARI_OLD_NAMES
#include <iostream>
#include "knuthYaoSampler.h"

// TODO : Use NTL to discrete gaussian sampling

using namespace std;

// Some functions for random sample generation
double Uniform(void) {
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double Normal(void) {
    return sqrt( -log(Uniform())*2.0 ) * sin( 2.0*M_PI*Uniform() );
}

double Gauss(double mu, double sigma) {
    double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
    return mu + sigma*z;
}

// Following the steps written under algorithm D in the paper - https://www.sav.sk/journals/uploads/0212094402follat.pdf
void SampleAlgorithmD(double mu, double sigma){
    
}

GEN Sample(int n, double sigma)
{
    GEN ret	= cgetg(n + 1, t_VEC);
    double z;
    int i;
    
    for (i = 1; i <= n; i++) {
        z = Gauss(0, sigma); z = abs(round(z));
        ret[i] = (long) stoi((long) z);
    }
    
    return ret;
}

GEN randomElement(int n){
    GEN ret;
    ret = cgetg(n + 1, t_VEC);
    for(int i=0; i<n; i++){
        gel(ret, i+1) = lift(gmodulo(stoi(rand()), stoi(300000)));
    }
    return ret;
}

struct pp{
    GEN q;
    GEN l;
    GEN p;
};

int main(){
    pari_init(2000000000,2);
    
    getProbabilityMatrix(4, "3.455555554534535353253425234543534535345245235312345678901234567890", 6, 4);
    cout<<endl;
    SampleKnuthYao(4, 6, 3, 4);
    
    GEN l, p, n, s, q;
    l = stoi(64); // l is the message length
    n = stoi(3530);
    // More demanding parameters
    //l = stoi(16128);
    //n = stoi(3530);
    
    int lambda = 212; // lambda is the security parameter for the homomorphic encryption
    
    // employing 128 bit security by taking n as 3530
    s = stoi(8);
    q = nextprime(gpowgs(stoi(2), lambda));
    p = gadd(gpowgs(stoi(2), 30), stoi(1));
    
    pp *pp1 = new pp;
    pp1->q = q;
    pp1->l = l;
    //ppq->p = p;
    
    
    
    
    cout<<"Parameter generation has been done"<<endl;
    
    cout<<"Cleaning up the Pari stack. Ending program.";
    pari_close();
    return 0;
}
