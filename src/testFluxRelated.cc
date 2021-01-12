#include "fluxRelated.h"
#include "converter.h"
#include <iostream>



int main() {

    std::vector<double> U{0, 0, 0, 0};
    std::vector<double> W{0.5197, -0.7259, 0.0, 0.4};

    for ( auto i: W) {
        std::cout << i << std::endl;
    }

    primitiveToConserved(W, U);
    double maxspeed = eigvalMaxMinF(U);
    std::cout << "The max eigval: " << maxspeed << std::endl;
    std::cout << "The min eigval: " << eigvalMaxMinF(U, true)<< std::endl;
    std::cout << "The max eigval (Y): " << eigvalMaxMinG(U) << std::endl;
    std::cout << "The max eigval (Y): " << eigvalMaxMinG(U, true) << std::endl;
    
    std::vector<double> flux1{0, 0, 0, 0};
    std::vector<double> flux2{0, 0, 0, 0};
    fluxF(U, flux1);
    fluxG(U, flux2);
    std::cout << "x-direction Flux: " << std::endl;
    for ( auto i: flux1) { 
        std::cout << i << std::endl;
    }
    std::cout << "y-direction Flux: " << std::endl;
    for ( auto i: flux2) { 
        std::cout << i << std::endl;
    }
    
    




return 0;
}
