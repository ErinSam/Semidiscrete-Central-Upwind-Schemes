"""
    Python script meant for prototyping the Mesh class for Semidiscrete Central Upwind
    Schemes
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor

from simple_colors import *

import sys
import numpy as np
import math as m
import gmshparser
import meshio

sys.path.append('./../../Mesh-Handlers-and-Pre-Processing/Topology')
import cartesian_geom as cgm



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
        
        Being used to store information for mesh consisting of quadrilateral elements
    """


    def __init__(self, mesh_file):
        # MAIN ATTRIBUTES
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


#        # INITIALISING CELL-VERTEX-SHARING NIEGHBOURS
#        print(yellow("Initialising cell-vertex-sharing topology data"))
#        for mesh_cell in self.cells:
#            for possible_neighbour in self.cells:
#                for vert_index in mesh_cell.vert_indices:
#                    if ( vert_index in possible_neighbour.vert_indices ):
#                        if ( (possible_neighbour.index != mesh_cell.index) & 
#                             (possible_neighbour.index not in mesh_cell.vertex_neighbours) ):
#                            mesh_cell.add_vertex_neighbours(possible_neighbour.index)
#        print(green("Completed initialising cell-vertex-sharing topology data"))
        
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

#        # Looping the list of boundarySurfaces to create all the boundaries
#        boundaryIndexCount = 0
#        while ( len(boundarySurfaces) > 0 ):
#            
#            newBoundary = boundary(boundaryIndexCount)
#            
#            # Looping until the newBoundary boundary object has been completely made
#            while( newBoundary.isOpen() ):
#                # Looping to find the next boundary element
#                for i, val in enumerate(boundarySurfaces):
#                    # Checking to see if the boundary element can be added
#                    if ( newBoundary.checkAppend(self, val) ):
#                        # the next boundary element for newBoundary has been found 
#                        # Removing this boundary element from boundarySurfaces
#                        del boundarySurfaces[i]
#                        break
#
#            # Adding the new boundary to the mesh's list of boundaries
#            self.boundaries.append(newBoundary)
#            boundaryIndexCount += 1
#        print(green("Completed initialising boundary data"))



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



if __name__ == "__main__":
    mesh = Mesh("./../Meshes/2D10x10.msh")

    print(yellow("\nFollowing are the boundary cells"))
    for i, mesh_cell in enumerate(mesh.cells):
        if ( mesh_cell.boundaryTag ):
            print(green("Cell index: "), mesh_cell.index)  
