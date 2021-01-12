// minmodForVLLim.cc

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


#include "minmodForVLLim.h"


double minmodScalar(const double& var1, const double& var2, const double& var3) {
    /**
     * Function that calculates the minmod of 3 variables 
     */
    if ((var1<0) && (var2<0) && (var3<0)){
        if ( (var1>var2) && (var1>var3) )
            return var1;
        else if ( (var2>var1) && (var2>var3) )
            return var2;
        else
            return var3;
    }
    else if ((var1>0) && (var2>0) && (var3>0)){
        if ( (var1<var2) && (var1<var3) )
            return var1;
        else if ( (var2<var1) && (var2<var3) )
            return var2;
        else 
            return var3;
    }
    else
        return 0.0;
}



void minmod(const std::vector<double>& var1,
            const std::vector<double>& var2, 
            const std::vector<double>& var3,
            std::vector<double>& mmd) {
    /** 
     * Function that calculates the minmod for 2 4-element vectors
     */

    for ( int i = 0; i < 4; i++ ) {
        mmd[i] = minmodScalar(var1[i], var2[i], var3[i]);
    }
}
