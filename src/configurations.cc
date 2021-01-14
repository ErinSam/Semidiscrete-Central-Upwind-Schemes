// configurations.cc

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


#include "configurations.h"



void initialConfig(int config, std::vector<double>& I, std::vector<double>& II, 
                   std::vector<double>& III, std::vector<double>& IV) {
    /**
     * Differnt initial configuration that have been taken from the numerical experiments 
     * that were conducted in
     * Alexander K. and Eitan T. 'Solution of two-dimensional reimann problems for gas dynamics
     * without reimann problem solvers'. (2012).
     */

    if ( config == 1 ) {
        std::vector<double> _I{1.0, 0.0, 0.0, 1.0};
        std::vector<double> _II{0.5197, -0.7259, 0.0, 0.4};
        std::vector<double> _III{0.1072, -0.7259, -1.4045, 0.0439};
        std::vector<double> _IV{0.2579, 0.0, -1.4045, 0.15};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 2 ) {
        std::vector<double> _I{1.0, 0.0, 0.0, 1.0};
        std::vector<double> _II{0.5197, -0.7259, 0.0, 0.4};
        std::vector<double> _III{1.0, -0.7259, -0.7259, 1.0};
        std::vector<double> _IV{0.5197, 0.0, -0.7259, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 3 ) {
        std::vector<double> _I{1.5, 0.0, 0.0, 1.5};
        std::vector<double> _II{0.5323, 1.206, 0.0, 0.3};
        std::vector<double> _III{0.138, 1.206, 1.206, 0.029};
        std::vector<double> _IV{0.5323, 0.0, 1.206, 0.3};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 4 ) {
        std::vector<double> _I{1.1, 0.0, 0.0, 1.1};
        std::vector<double> _II{0.5065, 0.8939, 0.0, 0.35};
        std::vector<double> _III{1.1, 0.8939, 0.8939, 1.1};
        std::vector<double> _IV{0.5065, 0.0, 0.8939, 0.35};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 5 ) {
        std::vector<double> _I{1.0, -0.75, -0.75, 1.0};
        std::vector<double> _II{2.0, -0.75, 0.5, 1.0};
        std::vector<double> _III{1.0, 0.75, 0.5, 1.0};
        std::vector<double> _IV{3.0, 0.75, -0.5, 1.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return; 
    }
    else if ( config == 6 ) {
        std::vector<double> _I{1.0, 0.75, -0.5, 1.0};
        std::vector<double> _II{2.0, 0.75, 0.5, 1.0};
        std::vector<double> _III{1.0, -0.75, 0.5, 1.0};
        std::vector<double> _IV{3.0, -0.75, -0.5, 1.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 7 ) {
        std::vector<double> _I{1.0, 0.1, 0.1, 1.0};
        std::vector<double> _II{0.5197, -0.6259, 0.1, 0.4};
        std::vector<double> _III{0.8, 0.1, 0.1, 0.4};
        std::vector<double> _IV{0.5197, 0.1, -0.6259, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 8 ) {
        std::vector<double> _I{0.5197, 0.1, 0.1, 0.4};
        std::vector<double> _II{1.0, -0.6259, 0.1, 1.0};
        std::vector<double> _III{0.8, 0.1, 0.1, 1.0};
        std::vector<double> _IV{1.0, 0.1, -0.6259, 1.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 9 ) {
        std::vector<double> _I{1.0, 0.0, 0.3, 1.0};
        std::vector<double> _II{2.0, 0.0, -0.3, 1.0};
        std::vector<double> _III{1.039, 0.0, -0.8133, 0.4};
        std::vector<double> _IV{0.5197, 0.0, -0.4259, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 10 ) {
        std::vector<double> _I{1.0, 0.0, 0.4297, 1.0};
        std::vector<double> _II{0.5, 0.0, 0.6076, 1.0};
        std::vector<double> _III{0.2281, 0.0, -0.6076, 0.3333};
        std::vector<double> _IV{0.4562, 0.0, -0.4297, 0.3333};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 11 ) {
        std::vector<double> _I{1.0, 0.1, 0.0, 1.0};
        std::vector<double> _II{0.5313, 0.8276, 0.0, 0.4};
        std::vector<double> _III{0.8, 0.1, 0.0, 0.4};
        std::vector<double> _IV{0.5313, 0.1, 0.7276, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 12 ) {
        std::vector<double> _I{0.5313, 0.0, 0.0, 0.4};
        std::vector<double> _II{1.0, 0.7276, 0.0, 1.0};
        std::vector<double> _III{0.8, 0.0, 0.0, 1.0};
        std::vector<double> _IV{1.0, 0.0, 0.7276, 1.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 13 ) {
        std::vector<double> _I{1.0, 0.0, -0.3, 1.0};
        std::vector<double> _II{2.0, 0.0, 0.3, 1.0};
        std::vector<double> _III{1.0625, 0.0, 0.8145, 0.4};
        std::vector<double> _IV{0.5313, 0.0, 0.4276, 0.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 14 ) {
        std::vector<double> _I{2.0, 0.0, -0.5606, 8.0};
        std::vector<double> _II{1.0, 0.0, -1.2172, 8.0};
        std::vector<double> _III{0.4736, 0.0, 1.2172, 2.6667};
        std::vector<double> _IV{0.9474, 0.0, 1.1606, 2.6667};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 15 ) {
        std::vector<double> _I{1.0, 0.1, -0.3, 1.0};
        std::vector<double> _II{0.5197, -0.6259, -0.3, 0.4};
        std::vector<double> _III{0.8, 0.1, -0.3, 0.4};
        std::vector<double> _IV{0.5313, 0.1, 0.4276, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return; 
    }
    else if ( config == 16 ) {
        std::vector<double> _I{0.5313, 0.1, 0.1, 0.4};
        std::vector<double> _II{1.0222, -0.6179, 0.1, 1.0};
        std::vector<double> _III{0.8, 0.1, 0.1, 1.0};
        std::vector<double> _IV{1.0, 0.1, 0.8276, 1.0};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 17 ) {
        std::vector<double> _I{1.0, 0.0, -0.4, 1.0};
        std::vector<double> _II{2.0, 0.0, -0.3, 1.0};
        std::vector<double> _III{1.0625, 0.0, 0.2145, 0.4};
        std::vector<double> _IV{0.5197, 0.0, -1.1259, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else if ( config == 18 ) {
        std::vector<double> _I{1.0, 0.0, 1.0, 1.0};
        std::vector<double> _II{2.0, 0.0, -0.3, 1.0};
        std::vector<double> _III{1.0625, 0.0, 0.2145, 0.4};
        std::vector<double> _IV{0.5197, 0.0, 0.2741, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
    else {
        std::vector<double> _I{1.0, 0.0, 0.3, 1.0};
        std::vector<double> _II{2.0, 0.0, -0.3, 1.0};
        std::vector<double> _III{1.0625, 0.0, 0.2145, 0.4};
        std::vector<double> _IV{0.5197, 0.0, -0.4259, 0.4};

        primitiveToConserved(_I, I);
        primitiveToConserved(_II, II);
        primitiveToConserved(_III, III);
        primitiveToConserved(_IV, IV);

        return;
    }
}
