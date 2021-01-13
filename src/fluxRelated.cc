// fluxRelated.cc

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


#include "fluxRelated.h"


double eigvalMaxMinF(const std::vector<double>& U, bool minim) {
    /**
     * Function that finds the maximum or the minimum eigenvalue (based on choice) of 
     * the Jacobian matrix of the x-direction flux term for Compressible 2D Euler 
     * Equations. Eigenvalues of \partial F/ \partial U. 
     * 
     * Args:
     *  @param1 : conserved flow field variables 
     *  @param2 : whether to provide max. or min. eigenvalue of the Jacobian
     *
     * Returns: 
     *  eval (double) : max. or min. eigenvalue of the Jacobian
     */     

    
    // For my ease
    const double u1 = U[0];
    const double u2 = U[1];
    const double u3 = U[2];
    const double u4 = U[3];

    // Creating the Jacobian matrix
    Eigen::Matrix<double, 4, 4> jacb;

    jacb(0,0) = 0; 
    jacb(0,1) = 1; 
    jacb(0,2) = 0; 
    jacb(0,3) = 0;
    
    jacb(1,0) = -pow(u2,2)/pow(u1,2) + 0.2*(pow(u2,2)/pow(u1,2) + pow(u3,2)/pow(u1,2)); 
    jacb(1,1) = 2*u2/u1 - 0.4*u2/u1; 
    jacb(1,2) = -0.4*u3/u1; 
    jacb(1,3) = 0.4;   

    jacb(2,0) = -u2*u3/pow(u1,2); 
    jacb(2,1) = u3/u1; 
    jacb(2,2) = u2/u1; 
    jacb(2,3) = 0.0;   

    jacb(3,0) = (-1.4*u4*u2/pow(u1,2) + 0.4*(pow(u2,3)/pow(u1,3) + pow(u3,3)/pow(u1,3))); 
    jacb(3,1) = (1.4*u4/u1 - 1.5*0.4*pow(u2,2)/pow(u1,2));
    jacb(3,2) = -1.5*0.4*pow(u3,2)/pow(u1,2); 
    jacb(3,3) = 1.4*u2/u1;   

    // Obtaining the max and the min of the eigenvalue
    Eigen::EigenSolver<Eigen::Matrix<double,4,4>> es(jacb, false);
    std::vector<double> eigvals;
    eigvals.reserve(4);
    for ( int i = 0; i<4; i++ ) {
        eigvals.push_back(es.eigenvalues()[i].real());
    }         

    if ( minim ) {
        return *min_element(eigvals.begin(), eigvals.end());
    }
    else {
        return *max_element(eigvals.begin(), eigvals.end());
    } 
}


double eigvalMaxMinG(const std::vector<double>& U, bool minim) {
    /**
     * Function that finds the maximum or the minimum eigenvalue (based on choice) of 
     * the Jacobian matrix of the y-direction flux term for Compressible 2D Euler 
     * Equations. Eigenvalues of \partial G/ \partial U. 
     *
     * Args:
     *  @param1 : conserved flow field variables 
     *  @param2 : whether to provide max. or min. eigenvalue of the Jacobian
     *
     * Returns: 
     *  eval (float) : max. or min. eigenvalue of the Jacobian
     */

    // For my ease
    const double u1 = U[0];
    const double u2 = U[1];
    const double u3 = U[2];
    const double u4 = U[3];

    // Creating the Jacobian matrix
    Eigen::Matrix<double, 4, 4> jacb;

    jacb(0,0) = 0; 
    jacb(0,1) = 0; 
    jacb(0,2) = 1; 
    jacb(0,3) = 0;
    
    jacb(1,0) = -u2*u3/pow(u1,2);
    jacb(1,1) = u3/u1;
    jacb(1,2) = u2/u1;
    jacb(1,3) = 0.0;

    jacb(2,0) = -pow(u3,2)/pow(u1,2) + 0.2*(pow(u2,2)/pow(u1,2) + pow(u3,2)/pow(u1,2));
    jacb(2,1) = -0.4*u3/u1;
    jacb(2,2) = 2*u3/u1 - 0.4*u3/u1;
    jacb(2,3) = 0.4;

    jacb(3,0) = -1.4*u4*u3/pow(u1,2) + 0.4*(pow(u2,3)/pow(u1,3) + pow(u3,3)/pow(u1,3));
    jacb(3,1) = -1.5*0.4*pow(u2,2)/pow(u1,2);
    jacb(3,2) = 1.4*u4/u1 - 1.5*0.4*pow(u3,2)/pow(u1,2);
    jacb(3,3) = 1.4*u3/u1;

    // Obtaining the max and the min of the eigenvalue
    Eigen::EigenSolver<Eigen::Matrix<double,4,4>> es(jacb, false);
    std::vector<double> eigvals;
    eigvals.reserve(4);
    for ( int i = 0; i<4; i++ ) {
        eigvals.push_back(es.eigenvalues()[i].real());
    }         

    if ( minim ) {
        return *min_element(eigvals.begin(), eigvals.end());
    }
    else {
        return *max_element(eigvals.begin(), eigvals.end());
    } 
}


void fluxF(const std::vector<double>& U, std::vector<double>& flux) {
    /**
     * Function that calculated the flux from a given conserved flowfield
     * variable vector in the x-direction for Compressible 2D Euler Equations.
     * F(U) from 
     *          U_t + F(U)_x + G(U)_y = 0
     *
     * Args:
     *  @param1 : conserved flow field variables 
     *  @param2 : flux vector
     *
     * Returns:
     *  flux (ndarray(4,)) : flux in the x-direction
     */

    flux[0] = U[1];
    flux[1] = pow(U[1],2)/U[0] + 0.4 * (U[3] - 0.5*pow(U[1],2)/U[0] - 0.5*pow(U[2],2)/U[0]);
    flux[2] = U[1]*U[2]/U[0];
    flux[3] = U[1]/U[0] * (1.4*U[3] - 0.4 * (0.5*pow(U[1],2)/U[0] - 0.5*pow(U[2],2)/U[0]));
}


void fluxG(const std::vector<double>& U, std::vector<double>& flux) {
    /**
     * Function that calculated the flux from a given conserved flowfield
     * variable vector in the x-direction for Compressible 2D Euler Equations.
     * G(U) from 
     *          U_t + F(U)_x + G(U)_y = 0
     *
     * Args:
     *  @param1 : conserved flow field variables 
     *  @param2 : flux vector
     *
     * Returns:
     *  flux  : flux in the y-direction
     */

    flux[0] = U[2];
    flux[1] = U[1]*U[2]/U[0]; 
    flux[2] = pow(U[2],2)/U[0] + 0.4 * (U[3] - 0.5*pow(U[1],2)/U[0] - 0.5*pow(U[2],2)/U[0]);
    flux[3] = U[2]/U[0] * (1.4*U[3] - 0.4 * (0.5*pow(U[1],2)/U[0] -0.5*pow(U[2],2)/U[0]));
}
