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

GEN power2(GEN x, int n, int kappa, int l, GEN q){
    GEN power2mat = zeromatcopy(n*kappa, l);
    long long int nkappa = n*kappa;
    GEN pow2 = stoi(1);
    for(int i = 1; i <= l; i++){
        //cout<<i<<endl;
        pow2 = stoi(1);
        for(int j=1; j <= kappa; j++){
            for(int k=1; k<=n; k++){
                gel(gel(power2mat, i), (j-1)*n+k) = gmodulo(gmul(gel(gel(x, i), k), pow2), q);
            }
            pow2 = gmul(pow2, stoi(2));
        }
    }
    
    return power2mat;
}

GEN appendmat(GEN m1, GEN m2, int col1, int col2, int row){
    GEN mat = zeromatcopy(row, col1+col2);
    for(int i =1; i<= col1+col2; i++){
        for(int j=1; j<=row; j++){
            if(i<=col1){
                gel(gel(mat, i), j) = gel(gel(m1, i), j);
            }
            else{
                gel(gel(mat, i), j) = gel(gel(m2, i-col1), j);
            }
        }
    }
    
    return mat;
}

GEN bits(GEN m, int kappa, int n){
    GEN mat = zeromatcopy(1, n*kappa);
    long long int nkappa = n*kappa;
    for(int i =1; i<= n; i++){
        GEN bintemp = binary_zv(gel(gel(m, i), 1));
        GEN binx = vecreverse(bintemp); // This is now LSB to MSB
        int size = lg(binx)-1;
        for(int j=1; j<=kappa; j++){
            if(j>size){
                gel(gel(mat, (i-1)*kappa+j), 1) = stoi(0);
            }
            else{
                gel(gel(mat, (i-1)*kappa+j), 1) = gel(gel(m, i), 1);
            }
        }
    }
    
    return mat;
}

int main(){
    pari_init(2500000000,2);
    
    
    //cout<<endl;
    //SampleKnuthYao(4, 6, 3, 4);
    
    GEN l, p, n, s, q;
    l = stoi(64); // l is the message length
    n = stoi(100);
    // More demanding parameters
    //l = stoi(16128);
    //n = stoi(3530);
    
    int lambda = 100; // lambda is the security parameter for the homomorphic encryption
    
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
    
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=1; j++){
            gel(gel(m, i), j) = stoi(2);
        }
    }
    
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
    
    // Additive Homomorphism
    
    /*GEN m1 = zeromatcopy(1, itos(l));
    
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=1; j++){
            gel(gel(m1, i), j) = stoi(3);
        }
    }
    
    GEN e1_1 = zeromatcopy(1, itos(n));
    GEN e2_1 = zeromatcopy(1, itos(n));
    GEN e3_1 = zeromatcopy(1, itos(l));
    
    for(int i = 1; i <= itos(n); i++){
        for(int j=1; j<=1; j++){
            //cout<<i<<" "<<j<<endl;
            gel(gel(e1_1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
            gel(gel(e2_1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
        }
    }
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=1; j++){
            gel(gel(e3_1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s), 4, 0, 6)), s));
        }
    }
    
    GEN c1_1, c2_1;
    c1_1 = gadd(RgM_mul(e1_1, A), gmul(p, e2_1));
    c2_1 = gadd(gadd(RgM_mul(e1_1, P), gmul(p, e3_1)), m1);
    
    GEN c1add, c2add;
    
    c1add = gadd(c1, c1_1);
    c2add = gadd(c2, c2_1);
    
    decryptedmessage = lift(gmodulo(lift(gadd(gmul(c1add, S), c2add)), p));
    cout<<"The decrypted message after additive homomorphism is "<<GENtostr(decryptedmessage)<<endl<<"---------------------------"<<endl;
    
    GEN cmul;
    
    // Append c2 after c1 to make just one ciphertext cmul
    GEN cbeforemul = zeromatcopy(1, itos(n)+itos(l));
    GEN c_1beforemul = zeromatcopy(1, itos(n)+itos(l));
    for(int i = 1; i <= itos(l)+itos(n); i++){
        for(int j=1; j<=1; j++){
            if(i<=itos(n)){
                gel(gel(cbeforemul, i), j) = gel(gel(c1, i), j);
                gel(gel(c_1beforemul, i), j) = gel(gel(c1_1, i), j);
                
            }
            else{
                gel(gel(cbeforemul, i), j) = gel(gel(c2, i-itos(n)), j);
                gel(gel(c_1beforemul, i), j) = gel(gel(c2_1, i-itos(n)), j);
            }
        }
    }
    
    cmul = RgM_transmul(cbeforemul, c_1beforemul);
    
    // This matrix can be precomputed
    GEN SIMatrix = zeromatcopy(itos(n)+itos(l), itos(l));
    GEN I = matid(itos(l));
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=itos(l)+itos(n); j++){
            if(j<=itos(n)){
                gel(gel(SIMatrix, i), j) = gel(gel(S, i), j);
                
            }
            else{
                gel(gel(SIMatrix, i), j) = gel(gel(I, i), j-itos(n));
            }
        }
    }
    
    decryptedmessage = lift(gmodulo(lift(gmul(RgM_transmul(SIMatrix, cmul), SIMatrix)), p));
    cout<<"The decrypted message after multiplicative homomorphism is "<<GENtostr(decryptedmessage)<<endl<<"---------------------------"<<endl;
    cout<<"Message matrix is "<<lg(gel(decryptedmessage, 1))-1<<"x"<<lg(decryptedmessage)-1<<endl;
    */
    
    // Performing key switching
    GEN n1 = stoi(90);
    GEN s1 = stoi(8);
    
    GEN R1, S1, A1;
    // R and S are nxl
    R1 = zeromatcopy(itos(n1), itos(l));
    S1 = zeromatcopy(itos(n1), itos(l));
    
    //cout<<GENtostr(gmodulo(R, s))<<endl;
    //cout<<GENtostr(R)<<endl;
    
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=itos(n1); j++){
            //cout<<i<<" "<<j<<endl;
            gel(gel(R1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
            gel(gel(S1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
        }
    }
    
    // A is nxn
    A1 = zeromatcopy(itos(n1), itos(n1));
    for(int i = 1; i <= itos(n1); i++){
        for(int j=1; j<=itos(n1); j++){
            // TODO remove this hardcoded modulo
            gel(gel(A1, i), j) = gmodulo(stoi(rand()%20), q);
        }
    }
    
    GEN P1, temp1;
    
    temp1 = RgM_mul(A1, S1);
    P1 = gsub(gmul(p, R1), temp1);
    // New set of keys have been generated
    //cout<<GENtostr(gdiv(glog(stoi(17), 4), mplog2(4)))<<endl;
    // NOTE: Add sufficient precision here if you get incorrect results.
    GEN kappa = gceil(gdiv(glog(q, 10), mplog2(10)));
    if(itos(kappa) != lambda+1){
        cout<<"log incorrect\n";
    }
    GEN X, E, Y;
    X = zeromatcopy(itos(gmul(n,kappa)), itos(n1));
    long long int nkappa = itos(gmul(n,kappa));
    for(int i = 1; i <= itos(n1); i++){
        for(int j=1; j<=nkappa; j++){
            // TODO remove this hardcoded modulo
            gel(gel(X, i), j) = gmodulo(stoi(rand()%30), q);
        }
    }
    E = zeromatcopy(nkappa, itos(l));
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=nkappa; j++){
            gel(gel(E, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
        }
    }
    
    Y = gsub(gadd(gmul(p, E), power2(S, itos(n), itos(kappa), itos(l), q)), gmul(X, S1));
    
    cout<<"Y matrix is "<<lg(gel(Y, 1))-1<<"x"<<lg(Y)-1<<endl;
    cout<<"Rotation keys have been generated"<<endl;
    
    GEN f1, f2, f3, E0, F, cdash;
    f1 = zeromatcopy(1, itos(n1));
    for(int i = 1; i <= itos(n1); i++){
        for(int j=1; j<=1; j++){
            gel(gel(f1, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
        }
    }
    f2 = zeromatcopy(1, itos(n1));
    for(int i = 1; i <= itos(n1); i++){
        for(int j=1; j<=1; j++){
            gel(gel(f2, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
        }
    }
    f3 = zeromatcopy(1, itos(l));
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=1; j++){
            gel(gel(f3, i), j) = lift(gmodulo(stoi(SampleKnuthYao(itos(s1), 4, 0, 6)), s1));
        }
    }
    E0 = gadd(gmul(f1, appendmat(A1, P1, itos(n1), itos(l), itos(n1))), gmul(p, appendmat(f2, f3, itos(n1), itos(l), 1)));
    
    F = appendmat(gmul(bits(c1, itos(kappa), itos(n)), X), gadd(gmul(bits(c1, itos(kappa), itos(n)), Y), c2), itos(n1), itos(l), 1);
    
    cdash = gadd(E0, F);
    //cout<<GENtostr(lift(cdash))<<endl;
    
    GEN SIMatrix = zeromatcopy(itos(n1)+itos(l), itos(l));
    GEN I = matid(itos(l));
    for(int i = 1; i <= itos(l); i++){
        for(int j=1; j<=itos(l)+itos(n1); j++){
            if(j<=itos(n1)){
                gel(gel(SIMatrix, i), j) = gel(gel(S1, i), j);
                
            }
            else{
                gel(gel(SIMatrix, i), j) = gel(gel(I, i), j-itos(n1));
            }
        }
    }
    decryptedmessage = lift(gmodulo(lift(gmul(F, SIMatrix)), p));
    cout<<GENtostr(decryptedmessage)<<endl;
    
    cout<<"Cleaning up the Pari stack. Ending program.";
    pari_close();
    return 0;
}
