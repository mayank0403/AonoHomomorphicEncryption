//
//  knuthYaoSampler.h
//  AonoEnc
//
//  Created by Mayank Rathee on 25/08/17.
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

// Function to generate the probability matrix which is a substitute of the DDT Tree. Some really fine optimizations are not considered here which help prune this matrix

void getProbabilityMatrix(int sigma, char* center, long long int precision, int tailprune){
    GEN temp;
    temp = strtor(center, precision);
    
    
}
