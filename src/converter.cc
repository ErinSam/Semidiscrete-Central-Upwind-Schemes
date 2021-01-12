// converter.cc

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor



#include "converter.h"



void conservedToPrimitive(const std::vector<double>& U, std::vector<double>& W) {
    /**
     * Function that converts a vector of conserved variable to vector of primitive variables
     * Meant for 2D compressible Euler equations 
     * 
     * Args:
     *  @param1 : conserved variables vector
     *  @param2 : primitive variables vector
     *
     */

    W[0] = U[0];
    W[1] = U[1] / U[0];
    W[2] = U[2] / U[0];
    W[3] = 0.4 * (U[3] - 0.5*U[0]*(pow(U[1],2) + pow(U[2],2)));
}


void primitiveToConserved(const std::vector<double>& W, std::vector<double>& U) {
    /**
     * Function that converts a vector of primitive variables to a vector of conserved variables. 
     * Meant for 2D compressible Euler equations
     * 
     * Args:
     *  @param1 : primitive variables vector
     *  @param2 : conserved variables vector 
     *
     */

    U[0] = W[0];
    U[1] = W[0]*W[1];
    U[2] = W[0]*W[2];
    U[3] = 2.5*W[3] + 0.5*W[0]*(pow(W[1], 2) + pow(W[2], 2));
}

