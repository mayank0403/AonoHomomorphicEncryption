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

// To know the amount of error in the code
void printErrorTerm(GEN p, GEN e1, GEN e2, GEN e3, GEN R, GEN S){
    cout<<"The error term should be small compared to q and it (the error term) is - "<<GENtostr(gmul(p, gadd(gadd(gmul(e1, R), gmul(e2, S)), e3)))<<endl<<"---------------------------"<<endl;
}

int main(){
    pari_init(2000000000,2);
    
    
    //cout<<endl;
    //SampleKnuthYao(4, 6, 3, 4);
    
    GEN l, p, n, s, q;
    l = stoi(64); // l is the message length
    n = stoi(530);
    // More demanding parameters
    //l = stoi(16128);
    //n = stoi(3530);
    
    int lambda = 212; // lambda is the security parameter for the homomorphic encryption
    
    // employing 128 bit security by taking n as 3530
    s = stoi(8);
    q = nextprime(gpowgs(stoi(2), lambda));
    p = gadd(gpowgs(stoi(2), 3), stoi(1));
    
    pp *pp1 = new pp;
    pp1->q = q;
    pp1->l = l;
    pp1->p = p;
    

    cout<<"Parameter generation has been done"<<endl;
    
    GEN R, S, A;
    // R and S are nxl
    R = zeromatcopy(itos(n), itos(l));
    S = zeromatcopy(itos(n), itos(l));
    
    //cout<<GENtostr(gmodulo(R, s))<<endl;
    //cout<<GENtostr(R)<<endl;
    getProbabilityMatrix(4, "0.00", 6, itos(s));
    cout<<pPack->startPos.size()<<endl;
    
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=itos(n); j++){
            //cout<<i<<" "<<j<<endl;
            gel(gel(R, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
            gel(gel(S, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
        }
    }
    
    // A is nxn
    A = zeromatcopy(itos(n), itos(n));
    for(int i = 1; i <= itos(n); i++){
        for(int j=1; j<=itos(n); j++){
            gel(gel(A, i), j) = gmodulo(stoi(rand()%20), q);
        }
    }
    
    GEN P, temp;
    
    temp = RgM_mul(A, S);
    P = gsub(gmul(p, R), temp);
    cout<<"Matrix is "<<lg(gel(P, 1))-1<<"x"<<lg(P)-1<<endl;
    cout<<"Key generation has been done\n";
    
    GEN m = zeromatcopy(1, itos(l));
    
    GEN e1 = zeromatcopy(1, itos(n));
    GEN e2 = zeromatcopy(1, itos(n));
    GEN e3 = zeromatcopy(1, itos(l));
    
    for(int i = 1; i <= itos(n); i++){
        for(int j=1; j<=1; j++){
            //cout<<i<<" "<<j<<endl;
            gel(gel(e1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
            gel(gel(e2, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
        }
    }
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=1; j++){
            gel(gel(e3, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
        }
    }
    cout<<"Errors generated\n";
    GEN c1, c2;
    c1 = gadd(RgM_mul(e1, A), gmul(p, e2));
    c2 = gadd(gadd(RgM_mul(e1, P), gmul(p, e3)), m);
    //cout<<GENtostr(c2)<<endl;
    cout<<"Message matrix is "<<lg(gel(c2, 1))-1<<"x"<<lg(c2)-1<<endl;
    cout<<"Decrypting the encrypted message now"<<endl;
    
    GEN decryptedmessage;
    decryptedmessage = lift(gmodulo(lift(gadd(gmul(c1, S), c2)), p));
    cout<<"The actual message is "<<GENtostr(m)<<endl<<"---------------------------"<<endl;
    cout<<"The decrypted message is "<<GENtostr(decryptedmessage)<<endl<<"---------------------------"<<endl;
    
    printErrorTerm(p, e1, e2, e3, R, S);
    
    cout<<"Cleaning up the Pari stack. Ending program.";
    pari_close();
    return 0;
}
