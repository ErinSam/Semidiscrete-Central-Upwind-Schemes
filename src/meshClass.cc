// meshClass.cc

/**
 * File contains the definition of the mesh class and the cell class from which the main 
 * application of semi-discrete central upwind schemes are used 
 */

// Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


#include "meshClass.h"


struct cell {
    /** 
     * Structure contains all relevant data for a mesh cell to be used in semi-discrete 
     * central upwind scheme
     *
     * Attributes : 
     *  center (std::vector<double>) : location of the center of the cell
     *  flowField (std::vector<double>) : flow field of conserved variables of the cell
     *  flowFieldTemp (std::vector<double>) : temp. storage of flowField to update the cell
     *  boundaryTag (bool) : whether the cell is on the outermost boundary layer 
     *  innerBoundaryTag (bool) : whether the cell is on the outer by one boundary layer
     */

    std::vector<double> center;
    std::vector<double> flowField{0.0, 0.0, 0.0, 0.0};
    std::vector<double> flowFieldTemp{0.0, 0.0, 0.0, 0.0}; 
    bool boundaryTag = false;
    bool innerBoundaryTag = false;

    cell(double x, double y) {
        center.push_back(x);
        center.push_back(y);
    }
};



class Mesh {
    /** 
     * This is the main mesh class. 
     * 
     * Private Attributes:
     *  cells (std::vector<struct cell>) : vector contains all the cells of the mesh
     *  dx (double) : mesh resolution in x-direction
     *  dy (double) : mesh resolution in y-direction 
     *  time (double) : the time at which the mesh stores values for 
     *
     * Private Functions: 
     *  vanLeerLimx(i,j)
     *  vanLeerLimy(i,j)
     *  cellNorth(i,j) 
     *  cellSouth(i,j)
     *  cellWest(i,j)
     *  cellEast(i,j)
     *  numFlux_x(i,j)
     *  numFlux_y(i,j)
     * 
     * Public Functions:
     *  updateCells(dt)
     *  applyBoundaryConditions()
     *  save()
     *  initialise()
     * 
     */

private:
    // Attribute
    std::vector<struct cell> cells;
    double dx, dy;
    double time;
    int cellCount;

    // Private Functions
    void 


public:

    Mesh(int CELLCOUNT, double DX, double DY) {
        cellCount = CELLCOUNT;
        dx = DX;
        dy = DY;
        time = 0.0;

        // Initialising the cell data 
        for ( int i=0; i<cellCount; i++ ) {
            for ( int j=0; j<cellCount; j++ ) {
                cells.push_back(cell(dx*(i+0.5), dy*(j+0.5)));
            }
        }        
    }

    void updateCells(double dt);
    void applyBoundaryConditions();
    void save();
    void initialise(std::vector<double> I, std::vector<double> II, std::vector<double> III,
                    std::vector<double> IV);
};


















































int main() {

    Mesh mesh(10, 0.1, 0.1);



return 0;
}

























