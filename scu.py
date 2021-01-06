"""
    Python script meant for prototyping the Mesh class for Semidiscrete Central Upwind
    Schemes
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor

from simple_colors import *
from IPython import embed

import sys
import numpy as np
import math as m
import pandas as pd
import gmshparser
import meshio

sys.path.append('./../Mesh-Handlers-and-Pre-Processing/Topology')
import cartesian_geom as cgm

sys.path.append('./dependencies/fluxRelated_lib')
import fluxRelated as frtd

sys.path.append('./dependencies/minmodForVLLim_lib')
import minmodForVLLim as vll

sys.path.append('./dependencies/interfaceSpeed_lib')
import interfaceSpeed as speed

sys.path.append('./dependencies/converter_lib')
import converter as convert



class vertex:
    """
        Class that stores relevant information about the vertex

        Data:
            loc: ndarray(2,); location of the vertex in the global coordinate system of the mesh
    """

    def __init__(self, index, locx, locy):
        self.index = index
        self.loc = np.array([locx, locy])



class surface:
    """
        Class that stores relevant information about the surface 

        Attributes: 
            index: int; global index of the surface
            vertices: ndarray(2,); indices of the vertices that make create the surface
            length: float; length of surface
            center: float; location of the center of the surface in the global coordinate
                           system of the mesh
            cell_indices: ndarray(2,); the global index of the cells that straddle the surface
            boundaryTag: (bool); False if the surface is not a boundary element and True otherwise
    """
    
    def __init__(self, index, vertices, vert1_loc, vert2_loc):
        self.index = index
        self.vertices = vertices
        self.length = cgm.dist(vert1_loc, vert2_loc)
        self.center = (vert1_loc + vert2_loc)/2
        self.straddling_cells = []
        self.boundaryTag = False

    def add_straddling_cell(self, cell_index):
        """ Method that adds one of the cells that straddle the surface
            Args:
                cell_index: int; index of the cell that straddles it
        """
        self.straddling_cells.append(cell_index)



class boundary():
    """
        Class that stores all data related to the boundary of the mesh. Main purpose is to serve
        as a means of applying boundary conditions

        Attributes:
            element: ndarray(boundaryLength,); indices of the surfaces that comprise the boundary
            type: string; the type that the boundary is (eg: Reflective, Neumann, Dirichlet)
            openEnds: int; the number of incomplete ends that the boundary has. essential in 
                           constructing boundaries
        
        Methods:
            isOpen(): 
            checkAppend():
            plot():
            assignType():
            merge():
    """

    def __init__(self, index):
        self.index = index
        self.element = np.array([])
        self.openEnds = 2
        self.type = "Unassigned"

    
    def isOpen(self):
        """ Method that returns True if the boundary construction is incomplete
        """
        if ( self.openEnds == 0 ):
            return False 
        elif ( self.openEnds > 0 ):
            return True 
        else:
            print(red("ERROR. self.openEnds for boundary %d has reached a negative value." 
                        %self.index))

    def checkAppend(self, mesh, surfaceIndex):
        """ Pivotal method that checks if the provided surface index neighbors the ends of the 
            current construction of the boundary. 
            If not, return false.
            If it does, return true and append the element at the front or back of element ndarray
            according to whether it is adjacent to the first or last boundary element. Check must
            be performed to see is these are concave boundary elements.
        
            Args:
                self
                mesh: (Mesh Object); required. Varied usage shall require usage for a variety of 
                                     classes that inherit from the main mesh class
                surfaceIndex: int; the index of the surface that we are checking it with

            Returns:
                Boolean True (or) False
        """
        
        # Checking all elements
        if ( (mesh.surfaces[surfaceIndex].vertices[0] not in mesh.surfaces[self.element[0]].vertices) & \
             (mesh.surfaces[surfaceIndex].vertices[1] not in mesh.surfaces[self.element[0]].vertices) & \
             (mesh.surfaces[surfaceIndex].vertices[0] not in mesh.surfaces[self.element[-1]].vertices) & \
             (mesh.surfaces[surfaceIndex].vertices[1] not in mesh.surfaces[self.element[-1]].vertices)):
            # Returning False since not neighboring
            return False

        else:
            # Checking for concavity (means they are of different boundaries)
            if ((mesh.surfaces[surfaceIndex].straddling_cells[0] in mesh.surfaces[self.element[0]].straddling_cells) | \
                (mesh.surfaces[surfaceIndex].straddling_cells[0] in mesh.surfaces[self.element[-1]].straddling_cells)):
                # surfaceIndex corresponds to element of different boundary
                self.openEnds -= 1
                return False

            else:
                if ( (mesh.surfaces[surfaceIndex].vertices[0] in mesh.surfaces[self.element[0]].vertices) | \
                    (mesh.surfaces[surfaceIndex].vertices[1] in mesh.surfaces[self.element[0]].vertices) ):
                    # Appending to the front of self.element ndarray
                    self.element = np.insert(self.element, 0, surfaceIndex, axis=0) 
                    return True

                elif ( (mesh.surfaces[surfaceIndex].vertices[0] in mesh.surfaces[self.element[-1]].vertices) | \
                    (mesh.surfaces[surfaceIndex].vertices[1] in mesh.surfaces[self.element[-1]].vertices) ):
                    # Appending to the back of self.element ndarray
                    self.element = np.append(self.element, surfaceIndex, axis=0) 
                    return True

                else:
                    print(red("ERROR. Huge mistake in the implementation of checkAppend() \
                            in boundary class"))
                    return None

    
    def plot(self, mesh):
        """ Method that plots the surface using a series of points in matplotlib. Helps visualize
            location and serves as a check
        """

        x = []
        y = []
        for i, val in enumerate(self.element):
            x.append(mesh.surfaces[val].center[0])
            y.append(mesh.surfaces[val].center[1])
        
        plt.scatter(x, y)
        plt.xlim(-1, 4)
        plt.ylim(-1, 4)
        plt.show()


    def imposeBC(self, mesh):
        """ Salient method of the boundary class. Most likely to be removed from this base class
            and implemented in problem specific situations
        """
        
        raise NotImplementedError



class cell:
    """
        Class that stores relevant information about the cell

        Data:
            index: int; index of the current cell
            vert_indices: ndarray(n,); indices of the vertices that make it up
            surface_indices: ndarray(4,); indices of the surfaces that enclose the cell
            area: float; area of the cell
            center: float; location of the cell center in the global coordinate system of the mesh
            vertex_neighbours: list; indices of the cells that share the same vertices as self
            surface_neighbours: list; indices of the cells that share the same surfaces as self
    """

    def __init__(self, index, vert_indices, vert1_loc, vert2_loc, vert3_loc, vert4_loc):
        self.index = index
        self.vert_indices = vert_indices 
        self.area = np.linalg.norm((vert2_loc - vert1_loc) * (vert3_loc - vert1_loc))
        self.center = 0.25 * (vert1_loc + vert2_loc + vert3_loc + vert4_loc)
        self.surface_indices = []
        self.vertex_neighbours = []
        self.surface_neighbours = []
        self.boundaryTag = False
        self.innerBoundaryTag = False
        self.flowField = np.zeros(4)
        self.flowFieldTemp = np.zeros(4)
    
    def add_surfaces(self, surface_index):
        """ Method that adds a surface for the cell
        """
        self.surface_indices.append(surface_index)

    def add_surface_neighbours(self, neb_index):
        """ Method that adds index of surface-sharing-cell to the neighbor
        """
        self.surface_neighbours.append(neb_index)

    def add_vertex_neighbours(self, neb_index):
        """ Method that adds index of a vertex-sharing-cell to the neighbor
        """
        self.vertex_neighbours.append(neb_index)



class Mesh: 
    """
        Class that stores relevant information of the whole mesh
        
        The mesh defined must be a structured, uniform (no grading!) 2D mesh
    """


    def __init__(self, mesh_file, dx, dy):
        # MAIN ATTRIBUTES
        self.time = 0.0
        self.dx = dx
        self.dy = dy
        self.vertices = [] 
        self.surfaces = []
        self.cells = []
        self.boundaries = []
        
        # Reading .msh file
        print(yellow("Reading in the mesh file"))
        print(mesh_file)
        mesh = meshio.read(mesh_file)        
        print(green("Completed reading in mesh file"))
        
        # INITIALISING VERTEX DATA 
        print(yellow("Initialising vertex topology data"))
        for i, val in enumerate(mesh.points):
            vertex_object = vertex(i, val[0], val[1])
            self.vertices.append(vertex_object)
        print(green("Completed initialising vertex topology data"))
        
        
        # INITIALISING CELL DATA
        print(yellow("Initialising cell topology data"))
        cell_count_index = 0
        for i, cell_block in enumerate(mesh.cells):
            if ( cell_block.type == "quad" ):
                for j, verti in enumerate(cell_block[1]):
                    new_mesh_cell = cell(cell_count_index, verti, self.vertices[verti[0]].loc,
                                                                  self.vertices[verti[1]].loc,
                                                                  self.vertices[verti[2]].loc,
                                                                  self.vertices[verti[3]].loc)
                    self.cells.append(new_mesh_cell)
                    # Incrementing the mesh cell index count
                    cell_count_index += 1
        print(green("Completed initialising cell topology data"))
                    

        # INITIALISING SURFACE DATA
        print(yellow("Initialising surface topology data"))
        surface_count = 0

        for mesh_cell in self.cells:
            # Obtaing surfaces that could be potentially added to the surface list
            trial_surfaces = np.array([[mesh_cell.vert_indices[0], mesh_cell.vert_indices[1]],
                                       [mesh_cell.vert_indices[1], mesh_cell.vert_indices[2]],
                                       [mesh_cell.vert_indices[2], mesh_cell.vert_indices[3]],
                                       [mesh_cell.vert_indices[3], mesh_cell.vert_indices[0]]])

            for i, trial_surface in enumerate(trial_surfaces):
                ver1 = trial_surface[0]
                ver2 = trial_surface[1]
                
                # Required to check in the case that we are checking the first element
                if ( len(self.surfaces) != 0):
                    surface_index = self.surface_check(ver1, ver2)
                else:
                    surface_index = -1

                if ( surface_index >= 0 ):
                    # Surface is already present

                    # Adding the current cell as a cell that straddles the surface 
                    self.surfaces[surface_index].add_straddling_cell(mesh_cell.index)

                    # Adding the surface to the list of surfaces that enclose the current cell
                    mesh_cell.add_surfaces(surface_index)

                else:
                    # Surface is not present. Adding new surface object to the list of surfaces 
                    surface_verts = [ver1, ver2]
                    new_surface = surface(surface_count, surface_verts, self.vertices[ver1].loc,
                                                                        self.vertices[ver2].loc)
                    self.surfaces.append(new_surface)
                
                    # Adding the current cell as a cell that straddles the surface
                    self.surfaces[surface_count].add_straddling_cell(mesh_cell.index)
               
                    # Adding the surface to the list of surfaces that enclose the current cell
                    mesh_cell.add_surfaces(surface_count)
                
                    # Incrementing surface count
                    surface_count += 1
        print(green("Completed initialising surface topology data"))
        

        # INITIALISING CELL-SURFACE-SHARING NIEGHBOURS
        print(yellow("Initialising cell-surface-sharing topology data"))
        for mesh_cell in self.cells:
            for surf_index in mesh_cell.surface_indices:
                for possible_neighbour in self.surfaces[surf_index].straddling_cells:
                    if (possible_neighbour != mesh_cell.index):
                        mesh_cell.add_surface_neighbours(possible_neighbour)
        print(green("Completed initialising cell-surface-sharing topology data"))

        
        # ESTABLISHING THE BOUNDARIES
        print(yellow("Initialising boundary data"))
        boundarySurfaces = []
        # Looping list of surfaces and obtaining the surfaces that form the boundaries
        for i, val in enumerate(self.surfaces):
            if ( len(val.straddling_cells) == 1 ):
                # Setting the boundaryTag as true for the cell
                self.surfaces[i].boundaryTag = True
                boundarySurfaces.append(i)

                # Updating boundaryTag of cell in contact with boundary
                boundCell = self.surfaces[i].straddling_cells[0]
                self.cells[boundCell].boundaryTag = True
        
        # Tagging the inner boundary (our mesh is padded twice)
        for i, innBoundCell in enumerate(self.cells):
            if innBoundCell.boundaryTag:
                for j, ind in enumerate(innBoundCell.surface_neighbours):
                    if not self.cells[ind].boundaryTag:
                        self.cells[ind].innerBoundaryTag = True

        print(green("Completed initialising boundary data"))



    def surface_check(self, vert1, vert2):
        """
            Method that checks whether a surface, given 2 vertices, are preset in the list of
            added vertices

            Args:
                vert1: int; index of first vertex 
                vert2: int; index of second vertex

            Returns: 
                surface_index: int; index of the surface if present else returns -1
        """
    

        # Checking the list of surfaces to see if the surface is already present 
        for surf in self.surfaces:
            if ( ((surf.vertices[0] == vert1) & (surf.vertices[1] == vert2)) |\
                 ((surf.vertices[1] == vert1) & (surf.vertices[0] == vert2)) ):

                # Surface is present; returning the its index
                return surf.index

        # Surface is not present
        return -1


    def vanLeerLimx(self, cellIndex):
        """
            UNTESTED

            Function that implements van Leer's one-parameter family of the minmod limiters
            to limit (prevent) numerical oscillations for the x-direction. 

            Args:
                self (Mesh) : 
                cellIndex (int) : the index of the cell for which we want to calculate

            Returns:
                flowField (ndarray(4,)) : the "limited" partial derivative of flowField in the 
                                          x-direction
        """
        # Arbitrarily taking the value of theta to be 1.25
        theta = 1.25

        # field_ : flowField of the current cell
        # fieldN : flowField of the N cell
        # fieldS : flowField of the S cell 
        field_ = self.cells[cellIndex].flowField
        fieldN = self.cells[self.cells[cellIndex].surface_neighbours[2]].flowField
        fieldS = self.cells[self.cells[cellIndex].surface_neighbours[0]].flowField

        # Calculating input for minmod
        var1 = theta * (fieldN - field_)/self.dx
        var2 = (fieldN - fieldS)/(2*self.dx)
        var3 = theta * (field_ - fieldS)/self.dx

        return vll.minmod(var1, var2, var3)


    def vanLeerLimy(self, cellIndex):
        """
            UNTESTED

            Function that implements van Leer's one-parameter family of the minmod limiters
            to limit (prevent) numerical oscillations for the y-direction. 

            Args:
                self (Mesh) : 
                cellIndex (int) : the index of the cell for which we want to calculate
                dy (float) : mesh cell size

            Returns:
                flowField (ndarray(4,)) : the "limited" partial derivative of flowField in the 
                                          y-direction
        """
        # Arbitrarily taking the value of theta to be 1.25
        theta = 1.25

        # field_ : flowField of the current cell
        # fieldW : flowField of the W cell
        # fieldE : flowField of the E cell 
        field_ = self.cells[cellIndex].flowField
        fieldW = self.cells[self.cells[cellIndex].surface_neighbours[3]].flowField
        fieldE = self.cells[self.cells[cellIndex].surface_neighbours[1]].flowField

        # Calculating input for minmod
        var1 = theta * (fieldW - field_)/self.dy
        var2 = (fieldW - fieldE)/(2*self.dy)
        var3 = theta * (field_ - fieldE)/self.dy

        return vll.minmod(var1, var2, var3)

    def cellNorth(self, cellIndex):
        """
            Function returns the flowField of the cell above the cell of the cellIndex

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the present cell

            Returns:
                field (ndarray(4,)) : flowField of the north cell
        """
        # Derivative
        derivative = self.vanLeerLimy(cellIndex)
        # Calculating the final 
        field = self.cells[cellIndex].flowField + self.dy/2 * derivative

        return field 


    def cellSouth(self, cellIndex):
        """
            Function returns the flowField of the cell below the cell of the cellIndex

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the present cell

            Returns:
                field (ndarray(4,)) : flowField of the north cell
        """
        # Derivative
        derivative = self.vanLeerLimy(cellIndex)
        # Calculating the final 
        field = self.cells[cellIndex].flowField - self.dy/2 * derivative

        return field 


    def cellWest(self, cellIndex):
        """
            Function returns the flowField of the cell to the left of the cell of the cellIndex

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the present cell

            Returns:
                field (ndarray(4,)) : flowField of the north cell
        """
        # Derivative
        derivative = self.vanLeerLimx(cellIndex)
        # Calculating the final 
        field = self.cells[cellIndex].flowField - self.dx/2 * derivative

        return field 


    def cellEast(self, cellIndex):
        """
            Function returns the flowField of the cell to the right of the cell of the cellIndex

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the present cell

            Returns:
                field (ndarray(4,)) : flowField of the north cell
        """
        # Derivative
        derivative = self.vanLeerLimx(cellIndex)
        # Calculating the final 
        field = self.cells[cellIndex].flowField + self.dx/2 * derivative

        return field 


    def numFlux_x(self, cellIndex):
        """
            Function that calculated the intercell numerical flux along the right face 
            of a quadrilateral 2D mesh element for the 2D compressible Euler Equations. 

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the cell

            Returns:
                numFlux (ndarray(4,)) : numerical flux for the right interface
        """
        # Right neighbouring cell's cell-index
        nebCellIndex = self.cells[cellIndex].surface_neighbours[1]

        # Obtaining the flow fields of the required cells
        fieldEast = self.cellEast(cellIndex)
        fieldWest = self.cellWest(nebCellIndex)

        # Obtaining the maximum and the minimum interface speed of the interface
        maxIntf = speed.x(fieldEast, fieldWest)
        minIntf = speed.x(fieldEast, fieldWest, minim=True)

        # Calculating numerical flux
        if ((maxIntf == 0) & (minIntf == 0)):
            numFlux = 0
        else:
            numFlux = (maxIntf * frtd.fluxF(fieldEast) - minIntf * frtd.fluxF(fieldWest)) \
                        / (maxIntf - minIntf) \
                       + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldWest - fieldEast)

        return numFlux


    def numFlux_y(self, cellIndex):
        """
            Function that calculates the intercell numerical flux along the top face of a 
            quadrilateral 2D mesh element for the 2D compressible Euler Equations.

            Args:
                self (Mesh) : 
                cellIndex (int) : index of the cell

            Returns:
                numFlux (ndarray(4,)) : numerical flux for the top interface
        """
        
        # Top cell's cell-index
        nebCellIndex = self.cells[cellIndex].surface_neighbours[2] 

        # Obtaining the flow fields of the required cells
        fieldNorth = self.cellNorth(cellIndex)
        fieldSouth = self.cellSouth(nebCellIndex)

        # Obtaining the maximum and minimum interface speeds 
        maxIntf = speed.y(fieldNorth, fieldSouth)
        minIntf = speed.y(fieldNorth, fieldSouth, minim=True)


        # Calculating numerical flux
        if ((maxIntf == 0) & (minIntf == 0)):
            numFlux = 0
        else:
            numFlux = (maxIntf * frtd.fluxG(fieldNorth) - minIntf * frtd.fluxG(fieldSouth)) \
                        / (maxIntf - minIntf) \
                        + (maxIntf*minIntf)/(maxIntf - minIntf) * (fieldSouth - fieldNorth)

        return numFlux


    def updateCells(self, dt):
        """
            UNTESTED

            Method that updates the cells according to the semi-discrete scheme (slightly incorrect
            because I am unsure how to handle the derivative term on the left hand side (2.4, [1]))

            Args:
                self (Mesh) : 

            Returns:
                None
        """
        # Setting dt (WARNING::WARNING)
        #dt = 0.0005

        # Setting dx and dy
        dx = self.dx
        dy = self.dy

        # Updating the Mesh time 
        self.time += dt

        # Looping over all the cells and updating if not boundary tagged 
        for i, elem in enumerate(mesh.cells):
            if not elem.boundaryTag and not elem.innerBoundaryTag:
                # Obtaining the index of the left cell
                leftCellIndex = elem.surface_neighbours[3]

                # Obtaining the index of the bottom cell
                bottomCellIndex = elem.surface_neighbours[0]

                # Calculating fluxes
                fluxLeft = self.numFlux_x(leftCellIndex)
                fluxRight = self.numFlux_x(i)
                fluxUp = self.numFlux_y(i)
                fluxDown = self.numFlux_y(bottomCellIndex)
        
                # Saving temporary flow Field
                mesh.cells[i].flowFieldTemp = mesh.cells[i].flowField \
                                                       - dt*(fluxRight - fluxLeft)/dx \
                                                       - dt*(fluxUp - fluxDown)/dy
                
        for i , elem in enumerate(mesh.cells):
            if not elem.boundaryTag and not elem.innerBoundaryTag:
                mesh.cells[i].flowField = mesh.cells[i].flowFieldTemp
        


    def applyBoundaryConditions(self):
        """
            UNTESTED

            Method that applies the zero gradient boundary conditons on the cells
            CURRENT IMPLEMENTATION IS UNABLE TO HANDLE CORNER BOUNDARY TAGGED ELEMENTS
            AND CORNER INNER BOUNDARY TAGGED ELEMENTS
            Note that the current implementation of the mesh has 2 layer of boundary

            Args:
                self (Mesh) : 

            Returns:
                None
        """
        # Looping over all the cells and updating the boundary tagged cells
        for i, elem in enumerate(mesh.cells):
            if (elem.innerBoundaryTag):
                for j, ind in enumerate(elem.surface_neighbours):
                    if not self.cells[ind].boundaryTag and not self.cells[ind].innerBoundaryTag:
                        # Updating the inner boundary element 
                        self.cells[i].flowField = self.cells[ind].flowField
                        break
                for j, ind in enumerate(elem.surface_neighbours):
                    if self.cells[ind].boundaryTag:
                        # Updating the outer boundary element
                        self.cells[ind].flowField = self.cells[i].flowField
                        break


    def initialise(self, I, II, III, IV):
        """
            UNTESTED 

            Function that initialises the cells according to provided data

            Args:
                self (Mesh) : 
                I (ndarray(4,)) : 1st Quadrant of Mesh
                II (ndarray(4,)) : 2nd Quadrant of Mesh
                III (ndarray(4,)) : 3rd Quadrant of Mesh
                IV (ndarray(4,)) : 4th Quadrant of Mesh

            Returns:
                None
        """
        # Looping over all the cells and initialising them 
        for i, elem in enumerate(mesh.cells):
            if ( (elem.center[0] > (0.5+1e-6)) & (elem.center[1] > (0.5+1e-6)) ):
                self.cells[elem.index].flowField = I
            elif ( (elem.center[0] < (0.5-1e-6)) & (elem.center[1] > (0.5+1e-6)) ):
                self.cells[elem.index].flowField = II
            elif ( (elem.center[0] > (0.5+1e-6)) & (elem.center[1] < (0.5-1e-6)) ):
                self.cells[elem.index].flowField = IV
            elif( (elem.center[0] < (0.5-1e-6)) & (elem.center[1] < (0.5-1e-6)) ):
                self.cells[elem.index].flowField = III


    def save(self):
        """
            UNTESTED

            Method that saves the data corresponding to the current time in HDF5 format 
            that can be later read in my an appropriate post-processor. The data frame save the 
            cell-Index, the location of the cell center as the x-coord and the y-coord (separately), 
            density, x-comp of velocity, y-comp of velocity and the pressure of all the cells in the 
            mesh. 

            The way in which data is saved also allows it to be reloaded to continue from a certain
            time step. This is a critical feature. 
            
            HDF5 file is saved in a directory './results/' and writes over previous saves (no 
            appending). File is saved with the time for which the data is being saved as the 
            filename: $(self.time).h5
            
            Args: 
                self (Mesh) : 

            Returns:
                None
        """

        # Creating a new dataframe for the given time step
        df = pd.DataFrame(columns=['index', 'x', 'y', 'rho', 'u', 'v', 'p'])

        for i, cell in enumerate(mesh.cells):
            # Converting the flowField into its primitive variables 
            W = convert.conservedToPrimitive(cell.flowField)
                
            # Writing in the data to the data frame 
            df = df.append({'index': cell.index, 
                            'x': cell.center[0],
                            'y': cell.center[1],
                            'rho': W[0], 
                            'u': W[1], 
                            'v': W[2], 
                            'p': W[3]}, ignore_index=True)

        # Writing the file 
        fileName = str(round(self.time,5))
        df.to_hdf("./results/" + fileName + ".h5", "table", append=False)



if __name__ == "__main__":
    mesh = Mesh("./Meshes/2D101x101.msh", 0.1, 0.1)

    # Defining the different test configurations
    def config1():
        I = np.array([1.0, 0.0, 0.0, 1.0])
        II = np.array([0.5197, -0.7259, 0.0, 0.4])
        III = np.array([0.1072, -0.7259, -1.4045, 0.0439])
        IV = np.array([0.2579, 0.0, -1.4045, 0.15])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config2():
        I = np.array([1.0, 0.0, 0.0, 1.0])
        II = np.array([0.5197, -0.7259, 0.0, 0.4])
        III = np.array([1.0, -0.7259, -0.7259, 1.0])
        IV = np.array([0.5197, 0.0, -0.7259, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config3():
        I = np.array([1.5, 0.0, 0.0, 1.5])
        II = np.array([0.5323, 1.206, 0.0, 0.3])
        III = np.array([0.138, 1.206, 1.206, 0.029])
        IV = np.array([0.5323, 0.0, 1.206, 0.3])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config4():
        I = np.array([1.1, 0.0, 0.0, 1.1])
        II = np.array([0.5065, 0.8939, 0.0, 0.35])
        III = np.array([1.1, 0.8939, 0.8939, 1.1])
        IV = np.array([0.5065, 0.0, 0.8939, 0.35])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config5():
        I = np.array([1.0, -0.75, -0.75, 1.0])
        II = np.array([2.0, -0.75, 0.5, 1.0])
        III = np.array([1.0, 0.75, 0.5, 1.0])
        IV = np.array([3.0, 0.75, -0.5, 1.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config6():
        I = np.array([1.0, 0.75, -0.5, 1.0])
        II = np.array([2.0, 0.75, 0.5, 1.0])
        III = np.array([1.0, -0.75, 0.5, 1.0])
        IV = np.array([3.0, -0.75, -0.5, 1.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config7():
        I = np.array([1.0, 0.1, 0.1, 1.0])
        II = np.array([0.5197, -0.6259, 0.1, 0.4])
        III = np.array([0.8, 0.1, 0.1, 0.4])
        IV = np.array([0.5197, 0.1, -0.6259, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config8():
        I = np.array([0.5197, 0.1, 0.1, 0.4])
        II = np.array([1.0, -0.6259, 0.1, 1.0])
        III = np.array([0.8, 0.1, 0.1, 1.0])
        IV = np.array([1.0, 0.1, -0.6259, 1.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config9():
        I = np.array([1.0, 0.0, 0.3, 1.0])
        II = np.array([2.0, 0.0, -0.3, 1.0])
        III = np.array([1.039, 0.0, -0.8133, 0.4])
        IV = np.array([0.5197, 0.0, -0.4259, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config10():
        I = np.array([1.0, 0.0, 0.4297, 1.0])
        II = np.array([0.5, 0.0, 0.6076, 1.0])
        III = np.array([0.2281, 0.0, -0.6076, 0.3333])
        IV = np.array([0.4562, 0.0, -0.4297, 0.3333])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config11():
        I = np.array([1.0, 0.1, 0.0, 1.0])
        II = np.array([0.5313, 0.8276, 0.0, 0.4])
        III = np.array([0.8, 0.1, 0.0, 0.4])
        IV = np.array([0.5313, 0.1, 0.7276, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config12():
        I = np.array([0.5313, 0.0, 0.0, 0.4])
        II = np.array([1.0, 0.7276, 0.0, 1.0])
        III = np.array([0.8, 0.0, 0.0, 1.0])
        IV = np.array([1.0, 0.0, 0.7276, 1.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config13():
        I = np.array([1.0, 0.0, -0.3, 1.0])
        II = np.array([2.0, 0.0, 0.3, 1.0])
        III = np.array([1.0625, 0.0, 0.8145, 0.4])
        IV = np.array([0.5313, 0.0, 0.4276, 0.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config14():
        I = np.array([2.0, 0.0, -0.5606, 8.0])
        II = np.array([1.0, 0.0, -1.2172, 8.0])
        III = np.array([0.4736, 0.0, 1.2172, 2.6667])
        IV = np.array([0.9474, 0.0, 1.1606, 2.6667])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config15():
        I = np.array([1.0, 0.1, -0.3, 1.0])
        II = np.array([0.5197, -0.6259, -0.3, 0.4])
        III = np.array([0.8, 0.1, -0.3, 0.4])
        IV = np.array([0.5313, 0.1, 0.4276, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config16():
        I = np.array([0.5313, 0.1, 0.1, 0.4])
        II = np.array([1.0222, -0.6179, 0.1, 1.0])
        III = np.array([0.8, 0.1, 0.1, 1.0])
        IV = np.array([1.0, 0.1, 0.8276, 1.0])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config17():
        I = np.array([1.0, 0.0, -0.4, 1.0])
        II = np.array([2.0, 0.0, -0.3, 1.0])
        III = np.array([1.0625, 0.0, 0.2145, 0.4])
        IV = np.array([0.5197, 0.0, -1.1259, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config18():
        I = np.array([1.0, 0.0, 1.0, 1.0])
        II = np.array([2.0, 0.0, -0.3, 1.0])
        III = np.array([1.0625, 0.0, 0.2145, 0.4])
        IV = np.array([0.5197, 0.0, 0.2741, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV

    def config19():
        I = np.array([1.0, 0.0, 0.3, 1.0])
        II = np.array([2.0, 0.0, -0.3, 1.0])
        III = np.array([1.0625, 0.0, 0.2145, 0.4])
        IV = np.array([0.5197, 0.0, -0.4259, 0.4])

        I = convert.primitiveToConserved(I)
        II = convert.primitiveToConserved(II)
        III = convert.primitiveToConserved(III)
        IV = convert.primitiveToConserved(IV)

        return I, II, III, IV



    I, II, III, IV = config2()
    mesh.initialise(I, II, III, IV)
    
    mesh.save()

    for i in range(1000000):
        print(yellow("\nMarching to time: "), (i+1)*0.000001)
        mesh.updateCells(0.000001)
        if ( i % 10 == 0):
            mesh.save()
            print(green("saved"))
     

