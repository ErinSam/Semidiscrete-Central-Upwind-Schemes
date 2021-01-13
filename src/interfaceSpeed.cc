// interfaceSpeed.cc

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


#include "interfaceSpeed.h"


double speedx(const std::vector<double>& fieldEast, const std::vector<double>& fieldWest, 
                bool minim) { 
    /**
     * Function that calculates the interface speed of the right interface of the cell.
     * 
     * Args:
     *  @param1 : the east flowField of the owner cell
     *  @param2 : the west flowField of the neighbor cell
     *  @param3 : the 
     * 
     * Returns:
     *  intfSpeed : the required interface speed
     */
    
    double eigval1, eigval2, intfSpeed;

    if ( minim ) {
        // Case that we want the minimum interface speed
        eigval1 = eigvalMaxMinF(fieldWest, true);
        eigval2 = eigvalMaxMinF(fieldEast, true);

        if ( (eigval1 < eigval2) && (eigval1 < 0.0) ) 
            intfSpeed = eigval1;
        else if ( (eigval2 < eigval1) && (eigval2 < 0.0) )
            intfSpeed = eigval2;
        else
            intfSpeed = 0.0;

        return intfSpeed;
    }
    else {
        // Case that we want the maximum interface speed
        eigval1 = eigvalMaxMinF(fieldWest, false);
        eigval2 = eigvalMaxMinF(fieldEast, false);

        if ( (eigval1 > eigval2) && (eigval1 > 0.0) ) 
            intfSpeed = eigval1;
        else if ( (eigval2 > eigval1) && (eigval2 > 0.0) )
            intfSpeed = eigval2;
        else
            intfSpeed = 0.0;

        return intfSpeed;
    } 
}


double speedy(const std::vector<double>& fieldNorth, const std::vector<double>& fieldSouth, 
                bool minim) {
    /**
     * Function that calculates the interface speed of the right interface of the cell.
     * 
     * Args:
     *  @param1 : the east flowField of the owner cell
     *  @param2 : the west flowField of the neighbor cell
     *  @param3 : the 
     * 
     * Returns:
     *  intfSpeed : the required interface speed
     */
   
    double eigval1, eigval2, intfSpeed;

    if ( minim ) {
        // Case that we want the minimum interface speed
        eigval1 = eigvalMaxMinG(fieldSouth, true);
        eigval2 = eigvalMaxMinG(fieldNorth, true);

        if ( (eigval1 < eigval2) && (eigval1 < 0.0) ) 
            intfSpeed = eigval1;
        else if ( (eigval2 < eigval1) && (eigval2 < 0.0) )
            intfSpeed = eigval2;
        else
            intfSpeed = 0.0;

        return intfSpeed;
    }
    else {
        // Case that we want the maximum interface speed
        eigval1 = eigvalMaxMinG(fieldSouth, false);
        eigval2 = eigvalMaxMinG(fieldNorth, false);

        if ( (eigval1 > eigval2) && (eigval1 > 0.0) ) 
            intfSpeed = eigval1;
        else if ( (eigval2 > eigval1) && (eigval2 > 0.0) )
            intfSpeed = eigval2;
        else
            intfSpeed = 0.0;

        return intfSpeed;
    }
}
