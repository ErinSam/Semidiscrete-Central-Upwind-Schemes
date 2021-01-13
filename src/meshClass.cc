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
    void vanLeerLimx(int i, int j, std::vector<double>& mmd);
    void vanLeerLimy(int i, int j, std::vector<double>& mmd); 
    void cellNorth(int i, int j, std::vector<double>& field);
    void cellSouth(int i, int j, std::vector<double>& field);
    void cellWest(int i, int j, std::vector<double>& field);
    void cellEast(int i, int j, std::vector<double>& field);
    void numFlux_x(int i, int j, std::vector<double>& flux);
    void numFlux_y(int i, int j, std::vector<double>& flux);

    void _updateCells(double dt);
    void _applyBoundaryConditions();
    void _save();
    void _initialise(int config);

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

        // Setting boundary tags 
        for ( int k = 0; k<cellCount; k++ ) {
            cells[cellCount*k].boundaryTag = true;
            cells[cellCount*k + cellCount-1].boundaryTag = true;
            cells[k].boundaryTag = true;
            cells[cellCount*(cellCount-1) + k].boundaryTag = true;

            if ( (k>0) && (k<cellCount-1) ) {
                cells[cellCount*k + 1].innerBoundaryTag = true;
                cells[cellCount*k + cellCount-2].innerBoundaryTag = true;
                cells[cellCount + k].innerBoundaryTag = true;
                cells[cellCount*(cellCount-2) + k].innerBoundaryTag = true;
            }
        }
    }

    void updateCells(double dt);
    void applyBoundaryConditions();
    void save();
    void initialise(int config);
};



void Mesh::vanLeerLimx(int i, int j, std::vector<double>& mmd) {
    /**
     * Function that implements van Leer's one-parameter family of the minmod limiters 
     * to limit (prevent) numerical oscillations for the x-direction
     * 
     * Args :
     *  @param1 : index 
     *  @param2 : index
     *  @param3 : the "limited" partial derivative of flowField in the x-direction
     */

    // Arbitrarily setting the value of theta to 1.25
    double theta = 1.25;

    const std::vector<double>& field_ = cells[i*cellCount + j].flowField;
    const std::vector<double>& fieldN = cells[i*cellCount + (j+1)].flowField;
    const std::vector<double>& fieldS = cells[i*cellCount + (j-1)].flowField;

    std::vector<double> var1, var2, var3;

    var1.push_back(theta * (fieldN[0] - field_[0])/dx);
    var1.push_back(theta * (fieldN[1] - field_[1])/dx);
    var1.push_back(theta * (fieldN[2] - field_[2])/dx);

    var2.push_back((fieldN[0] - fieldS[0])/(2*dx));
    var2.push_back((fieldN[1] - fieldS[1])/(2*dx));
    var2.push_back((fieldN[2] - fieldS[2])/(2*dx));

    var3.push_back(theta * (field_[0] - fieldS[0])/dx);
    var3.push_back(theta * (field_[1] - fieldS[1])/dx);
    var3.push_back(theta * (field_[2] - fieldS[2])/dx);
    
    minmod(var1, var2, var3, mmd); 
}


void Mesh::vanLeerLimy(int i, int j, std::vector<double>& mmd) {
    /**
     * Function that implements van Leer's one-parameter family of the minmod limiters 
     * to limit (prevent) numerical oscillations for the x-direction
     * 
     * Args :
     *  @param1 : index 
     *  @param2 : index
     *  @param3 : the "limited" partial derivative of flowField in the x-direction
     */

    // Arbitrarily setting the value of theta to 1.25
    double theta = 1.25;
    
    const std::vector<double>& field_ = cells[i*cellCount + j].flowField;
    const std::vector<double>& fieldW = cells[(i-1)*cellCount + j].flowField;
    const std::vector<double>& fieldE = cells[(i+1)*cellCount + j].flowField;

    std::vector<double> var1, var2, var3;

    var1.push_back(theta * (fieldW[0] - field_[0])/dy);
    var1.push_back(theta * (fieldW[1] - field_[1])/dy);
    var1.push_back(theta * (fieldW[2] - field_[2])/dy);

    var2.push_back((fieldW[0] - fieldE[0])/(2*dy));
    var2.push_back((fieldW[1] - fieldE[1])/(2*dy));
    var2.push_back((fieldW[2] - fieldE[2])/(2*dy));

    var3.push_back(theta * (field_[0] - fieldE[0])/dy);
    var3.push_back(theta * (field_[1] - fieldE[1])/dy);
    var3.push_back(theta * (field_[2] - fieldE[2])/dy);

    minmod(var1, var2, var3, mmd);
}


void Mesh::cellNorth(int i, int j, std::vector<double>& field) {
    /**
     * Function that returns the flowField of the cell above the cell of the cellIndex
     * 
     * Args :
     *  @param1 : index
     *  @param2 : index
     *  @param3 : flowField of the north cell
     */
    
    std::vector<double> mmd{0.0, 0.0, 0.0, 0.0};

    // Derivative 
    vanLeerLimy(i, j, mmd);

    // Calculating the final
    field[0] = cells[cellCount*i + j].flowField[0] + 0.5*dy*mmd[0]; 
    field[1] = cells[cellCount*i + j].flowField[1] + 0.5*dy*mmd[1]; 
    field[2] = cells[cellCount*i + j].flowField[2] + 0.5*dy*mmd[2]; 
    field[3] = cells[cellCount*i + j].flowField[3] + 0.5*dy*mmd[3]; 
}


void Mesh::cellSouth(int i, int j, std::vector<double>& field) { 
    /**
     * Function that returns the flowField of the cell above the cell of the cellIndex
     * 
     * Args :
     *  @param1 : index
     *  @param2 : index
     *  @param3 : flowField of the north cell
     */
    
    std::vector<double> mmd{0.0, 0.0, 0.0, 0.0};

    // Derivative 
    vanLeerLimy(i, j, mmd);
    
    // Calculating the final 
    field[0] = cells[cellCount*i + j].flowField[0] - 0.5*dy*mmd[0]; 
    field[1] = cells[cellCount*i + j].flowField[1] - 0.5*dy*mmd[1]; 
    field[2] = cells[cellCount*i + j].flowField[2] - 0.5*dy*mmd[2]; 
    field[3] = cells[cellCount*i + j].flowField[3] - 0.5*dy*mmd[3]; 
}


void Mesh::cellWest(int i, int j, std::vector<double>& field) { 
    /**
     * Function that returns the flowField of the cell to the left of the cell of the cellIndex
     * 
     * Args :
     *  @param1 : index
     *  @param2 : index
     *  @param3 : flowField of the west cell
     */
    
    std::vector<double> mmd{0.0, 0.0, 0.0, 0.0};

    // Derivative
    vanLeerLimx(i, j, mmd);
    
    // Calculating the final 
    field[0] = cells[cellCount*i + j].flowField[0] - 0.5*dx*mmd[0]; 
    field[1] = cells[cellCount*i + j].flowField[1] - 0.5*dx*mmd[1]; 
    field[2] = cells[cellCount*i + j].flowField[2] - 0.5*dx*mmd[2]; 
    field[3] = cells[cellCount*i + j].flowField[3] - 0.5*dx*mmd[3]; 
}    

    
void Mesh::cellEast(int i, int j, std::vector<double>& field) { 
    /**
     * Function that returns the flowField of the cell to the right of the cell of the cellIndex
     * 
     * Args :
     *  @param1 : index
     *  @param2 : index
     *  @param3 : flowField of the east cell
     */
    
    std::vector<double> mmd{0.0, 0.0, 0.0, 0.0};

    // Derivative
    vanLeerLimx(i, j, mmd);
    
    // Calculating the final 
    field[0] = cells[cellCount*i + j].flowField[0] + 0.5*dx*mmd[0]; 
    field[1] = cells[cellCount*i + j].flowField[1] + 0.5*dx*mmd[1]; 
    field[2] = cells[cellCount*i + j].flowField[2] + 0.5*dx*mmd[2]; 
    field[3] = cells[cellCount*i + j].flowField[3] + 0.5*dx*mmd[3]; 
}    


void Mesh::numFlux_x(int i, int j, std::vector<double>& flux) {
    /**
     * Function that calculated the intercell numerical flux along the right face of a 
     * quadrilateral 3D mesh element for the 2D compressible Euler equations
     * 
     * Args:
     *  @param1 : index
     *  @param2 : index
     *  @param3 : intercell numerical flux 
     */

    std::vector<double> fieldEast{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fieldWest{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxFEast{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxFWest{0.0, 0.0, 0.0, 0.0};

    // Obtaining the flow fields of the required cells
    cellEast(i, j, fieldEast);
    cellWest(i+1, j, fieldWest);

    // Obtaining the maximum and the minimum interface speed of the interface 
    double maxIntf = speedx(fieldEast, fieldWest, false);
    double minIntf = speedx(fieldEast, fieldWest, true);

    // Obtaining the fluxF for fieldEast and fieldWest
    fluxF(fieldEast, fluxFEast);
    fluxF(fieldWest, fluxFWest);

    // Calculating numerical flux 
    if ((maxIntf == 0) && (minIntf == 0)) {
        flux[0] = 0.0;
        flux[1] = 0.0;
        flux[2] = 0.0;
        flux[3] = 0.0;
    }
    else { 
        flux[0] = (maxIntf * fluxFEast[0] - minIntf * fluxFWest[0]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldWest[0] - fieldEast[0]);
        flux[1] = (maxIntf * fluxFEast[1] - minIntf * fluxFWest[1]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldWest[1] - fieldEast[1]);
        flux[2] = (maxIntf * fluxFEast[2] - minIntf * fluxFWest[2]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldWest[2] - fieldEast[2]);
        flux[3] = (maxIntf * fluxFEast[3] - minIntf * fluxFWest[3]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldWest[3] - fieldEast[3]);
    }
}


void Mesh::numFlux_y(int i, int j, std::vector<double>& flux) { 
    /**
     * Function that calculated the intercell numerical flux along the right face of a 
     * quadrilateral 3D mesh element for the 2D compressible Euler equations
     * 
     * Args:
     *  @param1 : index
     *  @param2 : index
     *  @param3 : intercell numerical flux 
     */
    
    std::vector<double> fieldNorth{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fieldSouth{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxGNorth{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxGSouth{0.0, 0.0, 0.0, 0.0};

    // Obtaining the flow fields of the required cells
    cellNorth(i, j, fieldNorth);
    cellSouth(i, j+1, fieldSouth);

    // Obtaining the maximum and minimum interace speed of the interface 
    double maxIntf = speedy(fieldNorth, fieldSouth, false);
    double minIntf = speedy(fieldNorth, fieldSouth, true);

    // Obtaining the fluxG for fieldNorth and fieldSouth
    fluxG(fieldNorth, fluxGNorth);
    fluxG(fieldSouth, fluxGSouth);

    // Calculating the numerical flux 
    if ((maxIntf == 0) && (minIntf == 0)) {
        flux[0] = 0.0;
        flux[1] = 0.0;
        flux[2] = 0.0;
        flux[3] = 0.0;
    }
    else {
        flux[0] = (maxIntf * fluxGNorth[0] - minIntf * fluxGSouth[0]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldSouth[0] - fieldNorth[0]);
        flux[1] = (maxIntf * fluxGNorth[1] - minIntf * fluxGSouth[1]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldSouth[1] - fieldNorth[1]);
        flux[2] = (maxIntf * fluxGNorth[2] - minIntf * fluxGSouth[2]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldSouth[2] - fieldNorth[2]);
        flux[3] = (maxIntf * fluxGNorth[3] - minIntf * fluxGSouth[3]) / (maxIntf - minIntf)
                    + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldSouth[3] - fieldNorth[3]);
    }
}


void Mesh::_updateCells(double dt) {
    /**
     * Function that updates the cells according to the semi-discrete scheme and to the best 
     * of my interpretation of it
     * 
     * Args : 
     *  @param1 : the size of the time step
     */
    
    // Updating the mesh time 
    time += dt;

    // Looping over all the cells and updating if it is not boundary tagged 
    std::vector<double> fluxLeft{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxRight{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxUp{0.0, 0.0, 0.0, 0.0};
    std::vector<double> fluxDown{0.0, 0.0, 0.0, 0.0};

    for ( int i = 0; i<cellCount; i++ ) {
        for ( int j = 0; j<cellCount; j++ ) {
            if ( !cells[cellCount*i + j].boundaryTag && !cells[cellCount*i + j].innerBoundaryTag ) {
                // Calculating the fluxes
                numFlux_x(i-1, j, fluxLeft);
                numFlux_x(i, j, fluxRight);
                numFlux_y(i, j, fluxUp);
                numFlux_y(i, j-1, fluxDown);
            
                // Saving temporary flow 
                cells[cellCount*i + j].flowFieldTemp[0] = cells[cellCount*i + j].flowField[0]
                                                          - dt*(fluxRight[0] - fluxLeft[0])/dx
                                                          - dt*(fluxUp[0] - fluxDown[0])/dy;
                cells[cellCount*i + j].flowFieldTemp[1] = cells[cellCount*i + j].flowField[1]
                                                          - dt*(fluxRight[1] - fluxLeft[1])/dx
                                                          - dt*(fluxUp[1] - fluxDown[1])/dy;
                cells[cellCount*i + j].flowFieldTemp[2] = cells[cellCount*i + j].flowField[2]
                                                          - dt*(fluxRight[2] - fluxLeft[2])/dx
                                                          - dt*(fluxUp[2] - fluxDown[2])/dy;
                cells[cellCount*i + j].flowFieldTemp[0] = cells[cellCount*i + j].flowField[3]
                                                          - dt*(fluxRight[3] - fluxLeft[3])/dx
                                                          - dt*(fluxUp[3] - fluxDown[3])/dy;
            }
        }
    }

  # pragma omp parallel for 
    for ( int i = 0; i<cellCount*cellCount; i++ ) { 
        if ( !cells[i].boundaryTag && !cells[i].innerBoundaryTag ) {
            cells[i].flowField[0] = cells[i].flowFieldTemp[0];
            cells[i].flowField[1] = cells[i].flowFieldTemp[1];
            cells[i].flowField[2] = cells[i].flowFieldTemp[2];
            cells[i].flowField[3] = cells[i].flowFieldTemp[3];
        }
    }
}
    

void Mesh::_applyBoundaryConditions() {
    /**
     * Function that applies the zro gradient boundar conditions on the cells 
     * CURRENT IMPLMENTATION IS UNABLE TO HANDLE CORNER BOUNDARY TAGGED ELEMENTS 
     * AND CORNER INNER BOUNDARY TAGGED ELEMENTS
     * Note that the current implementation of the mesh has 2 layers of boundary
     * 
     */

    // Boundary conditions on the inner boundary elements
  # pragma omp parallel for 
    for ( int k = 2; k < (cellCount-2); k++ ) {
        // Left Boundary 
        cells[cellCount + k].flowField[0] = cells[cellCount*2 + k].flowField[0];    
        cells[cellCount + k].flowField[1] = cells[cellCount*2 + k].flowField[1];    
        cells[cellCount + k].flowField[2] = cells[cellCount*2 + k].flowField[2];    
        cells[cellCount + k].flowField[3] = cells[cellCount*2 + k].flowField[3];    

        // Right Boundary
        cells[cellCount*(cellCount-2)+k].flowField[0]=cells[cellCount*(cellCount-3)+k].flowField[0];
        cells[cellCount*(cellCount-2)+k].flowField[1]=cells[cellCount*(cellCount-3)+k].flowField[1];
        cells[cellCount*(cellCount-2)+k].flowField[2]=cells[cellCount*(cellCount-3)+k].flowField[2];
        cells[cellCount*(cellCount-2)+k].flowField[3]=cells[cellCount*(cellCount-3)+k].flowField[3];

        // Bottom Boundary
        cells[cellCount*k + 1].flowField[0] = cells[cellCount*k + 2].flowField[0]; 
        cells[cellCount*k + 1].flowField[1] = cells[cellCount*k + 2].flowField[1]; 
        cells[cellCount*k + 1].flowField[2] = cells[cellCount*k + 2].flowField[2]; 
        cells[cellCount*k + 1].flowField[3] = cells[cellCount*k + 2].flowField[3]; 

        // Top Boundary
        cells[cellCount*k + cellCount-2].flowField[0]=cells[cellCount*k + cellCount-3].flowField[0];
        cells[cellCount*k + cellCount-2].flowField[1]=cells[cellCount*k + cellCount-3].flowField[1];
        cells[cellCount*k + cellCount-2].flowField[2]=cells[cellCount*k + cellCount-3].flowField[2];
        cells[cellCount*k + cellCount-2].flowField[3]=cells[cellCount*k + cellCount-3].flowField[3];
    }
    
    // Boudary conditions on the outer boundary elements
  # pragma omp parallel for 
    for ( int k = 1; k < (cellCount-1); k++ ) {
        // Left Boundary 
        cells[k].flowField[0] = cells[cellCount + k].flowField[0];    
        cells[k].flowField[1] = cells[cellCount + k].flowField[1];    
        cells[k].flowField[2] = cells[cellCount + k].flowField[2];    
        cells[k].flowField[3] = cells[cellCount + k].flowField[3];    

        // Right Boundary
        cells[cellCount*(cellCount-1)+k].flowField[0]=cells[cellCount*(cellCount-2)+k].flowField[0];
        cells[cellCount*(cellCount-1)+k].flowField[1]=cells[cellCount*(cellCount-2)+k].flowField[1];
        cells[cellCount*(cellCount-1)+k].flowField[2]=cells[cellCount*(cellCount-2)+k].flowField[2];
        cells[cellCount*(cellCount-1)+k].flowField[3]=cells[cellCount*(cellCount-2)+k].flowField[3];

        // Bottom Boundary
        cells[cellCount*k].flowField[0] = cells[cellCount*k + 1].flowField[0]; 
        cells[cellCount*k].flowField[1] = cells[cellCount*k + 1].flowField[1]; 
        cells[cellCount*k].flowField[2] = cells[cellCount*k + 1].flowField[2]; 
        cells[cellCount*k].flowField[3] = cells[cellCount*k + 1].flowField[3]; 

        // Top Boundary
        cells[cellCount*k + cellCount-1].flowField[0]=cells[cellCount*k + cellCount-2].flowField[0];
        cells[cellCount*k + cellCount-1].flowField[1]=cells[cellCount*k + cellCount-2].flowField[1];
        cells[cellCount*k + cellCount-1].flowField[2]=cells[cellCount*k + cellCount-2].flowField[2];
        cells[cellCount*k + cellCount-1].flowField[3]=cells[cellCount*k + cellCount-2].flowField[3];
    }
}        


void Mesh::_initialise(int config) { 
    /**
     * Function that initialises the cells according to provided data 
     * 
     * Args : 
     *  @param1 : 1st quadrant of the mesh
     *  @param2 : 2nd quadrant of the mesh
     *  @param3 : 3rd quadrant of the mesh
     *  @param4 : 4th quadrant of the mesh
     */


    std::vector<double> I{0.0, 0.0, 0.0, 0.0};
    std::vector<double> II{0.0, 0.0, 0.0, 0.0};
    std::vector<double> III{0.0, 0.0, 0.0, 0.0};
    std::vector<double> IV{0.0, 0.0, 0.0, 0.0};

    initialConfig(config, I, II, III, IV);


    // Looping over all the cells to initialise them 
    // Ist Quadrant 
  # pragma omp parallel for collapse(2)        
    for ( int i = cellCount/2; i < cellCount; i++ ) { 
        for ( int j = cellCount/2; j < cellCount; j++ ) {
            cells[cellCount*i + j].flowField[0] = I[0];
            cells[cellCount*i + j].flowField[1] = I[1];
            cells[cellCount*i + j].flowField[2] = I[2];
            cells[cellCount*i + j].flowField[3] = I[3];
        }
    } 

    // IInd Quadrant 
  # pragma omp parallel for collapse(2)        
    for ( int i = 0; i < cellCount/2; i++ ) { 
        for ( int j = cellCount/2; j < cellCount; j++ ) {
            cells[cellCount*i + j].flowField[0] = II[0];
            cells[cellCount*i + j].flowField[1] = II[1];
            cells[cellCount*i + j].flowField[2] = II[2];
            cells[cellCount*i + j].flowField[3] = II[3];
        }
    } 

    // IIIrd Quadrant 
  # pragma omp parallel for collapse(2)        
    for ( int i = 0; i < cellCount/2; i++ ) { 
        for ( int j = 0; j < cellCount/2; j++ ) {
            cells[cellCount*i + j].flowField[0] = III[0];
            cells[cellCount*i + j].flowField[1] = III[1];
            cells[cellCount*i + j].flowField[2] = III[2];
            cells[cellCount*i + j].flowField[3] = III[3];
        }
    } 

    // IVth Quadrant 
  # pragma omp parallel for collapse(2)        
    for ( int i = cellCount/2; i < cellCount; i++ ) { 
        for ( int j = 0; j < cellCount/2; j++ ) {
            cells[cellCount*i + j].flowField[0] = IV[0];
            cells[cellCount*i + j].flowField[1] = IV[1];
            cells[cellCount*i + j].flowField[2] = IV[2];
            cells[cellCount*i + j].flowField[3] = IV[3];
        }
    } 
}


void Mesh::_save() {
    /** 
     * Function that save all the required mesh details
     * Primitive flow field variables are sotred in .csv files where the file name is the 
     * present simulation time and 
     */

    // Creating the file and opening it 
    std::ofstream saveFile;
    std::string filename = std::to_string(time) + ".csv";
    saveFile.open(filename);

    // Writing the file headers
    saveFile << "index,x,y,rho,u,v,p" << std::endl;

    // Creating the temp variable for primitve variable vectors
    std::vector<double> prim{0.0, 0.0, 0.0, 0.0};

    // Writing all the data into the file
    for ( int i = 0; i<cellCount*cellCount; i++ ) { 
        conservedToPrimitive(cells[i].flowField, prim);
        saveFile << i+1 << "," << cells[i].center[0] << ","
                               << cells[i].center[1] << ","
                               << prim[0] << ","
                               << prim[1] << ","
                               << prim[2] << ","
                               << prim[3] << std::endl;
    }
        
    // Closing the outpute file
    saveFile.close();
} 

void Mesh::updateCells(double dt) {
    _updateCells(dt);
}

void Mesh::applyBoundaryConditions() {
    _applyBoundaryConditions();
}

void Mesh::save() {
    _save();
}

void Mesh::initialise(int config) {
    _initialise(config);
}
    


// pybind11 to make python module
PYBIND11_MODULE(meshClass, m) {
    namespace py = pybind11;

    m.doc() = "Contains class definition for a unfiorm, quadrilateral mesh that uses \
                Semi-Discrete Central Upwind Scheme to solve the 2D compressible Euler Equations";

    py::class_<Mesh>(m, "Mesh")
        .def(py::init<int, double, double>())
        .def("update", &Mesh::updateCells)
        .def("applyBC", &Mesh::applyBoundaryConditions)
        .def("save", &Mesh::save)
        .def("initialise", &Mesh::initialise);
}

