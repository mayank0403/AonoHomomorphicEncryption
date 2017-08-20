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



int main(){
    pari_init(2000000000,2);
    GEN l, p, n, s, q;
    l = stoi(64); // l is the message length
    n = stoi(3530);
    // More demanding parameters
    //l = stoi(16128);
    //n = stoi(3530);
    
    int lambda = 128; // lambda is the security parameter for the homomorphic encryption
    q = nextprime(gpowgs(stoi(2), lambda));
    p = stoi(2);
    
    cout<<"Parameter generation has been done"<<endl;
    
    cout<<"Cleaning up the Pari stack. Ending program.";
    pari_close();
    return 0;
}
