
"""
    Mesh Generator for the IGA-DG code (2D)
    
Inputs:
    - Dimension of the grid (along x and y directions)
    - Order of the solution (P)
    - Number of processors the mesh will be run on
    - The domain of the mesh
    - Type of flow to be solved:
        - Periodic Vortex
        - Inviscid Channel

Outputs:
    - A mesh file and solution point file will be output that
    satisfy the given input parameters.
    - If configured (set CONST_Plot to True), a figure of the mesh will
    be shown using Matplotlib on a seperate window once the grid has been 
    generated.
    
Supported BCs:
    - Periodic
    - SlipWall
    - TotalTemperaturePressure
    - BackPressure

"""

import matplotlib.pyplot as plt
import math
    

#=================================================================

                        #INPUT PARAMETERS

CONST_DIMENSION_X = 6 #Number of elements along x coordinate direction
CONST_DIMENSION_Y = 1 #Number of elements along y coordinate direction

# This gives the order p of the solution approximation.
#   (ex: If p=2, then 3 solution points created along each
#   coordinate direction for the element)
CONST_P = 1

# Number of processors that will be used with the mesh.
# The number of processors must be a power of 2.
CONST_NumProcessors = 4

# The domain of the rectangular (undeformed) mesh
CONST_xMin = -4.0
CONST_xMax = 4.0
CONST_yMin = 0.0
CONST_yMax = 1.5

# Set to True if the user wants a visual representation of the mesh (needs matplotlib)
CONST_Plot = True
CONST_PlotXRange = [-4.2, 4.2]
CONST_PlotYRange = [0.0, 4]

# If the mesh will be deformed
CONST_Deform = True

# If a label needs to be added to the mesh's file name
CONST_Label = ""

CONST_CaseType = "InviscidChannel"  # PeriodicVortex, InviscidChannel
#=================================================================







#================================================================

                    #THE DEFORMATION FUNCTION

# The sinusoidal perturbation function that is used to perturb to
# deform the rectangular mesh. It has as parameters a (which is
# either x or y) . In addition, the parameters n, A and Lo define the
# number of periods and the curvature of the perturbation.

CONST_n = 4.
CONST_A = 0.4
CONST_Lo = 16.
def sinPertrubFunction(a):
    return CONST_A*math.sin((CONST_n/CONST_Lo)*math.pi*a)


CONST_a = 0.5
CONST_b = 0.0
CONST_c = 0.6
def gaussianBumpPerturbFunction(x, j, numJ):

    """
    j = j index (starting with 0 being bottom most node)
    """

    aVal = CONST_a*(numJ-j-1)/(numJ-1)

    return aVal*math.exp(-1.*(((x-CONST_b)**2)/(2.*CONST_c**2)))


# This method is called if the grid must be deformed. It will loop over
# all the solution points and node points of the grid.
def perturbGrid(Mesh, MeshNumPointsX, MeshNumPointsY):
    
    for i in range(MeshNumPointsX):
        for j in range(MeshNumPointsY):
            
            xOld = Mesh[i][j][0]
            yOld = Mesh[i][j][1]
            
            """
            dx = sinPertrubFunction(yOld)
            dy = sinPertrubFunction(xOld)
            """

            """
            dy = gaussianBumpPerturbFunction(xOld, j, MeshNumPointsY)
            dx = 0.0

            xNew = xOld + dx
            yNew = yOld + dy
            
            Mesh[i][j][0] = xNew
            Mesh[i][j][1] = yNew
            """

            exactPointTuple = gaussianBumpExact(i,j)
            Mesh[i][j][0] = exactPointTuple[0]
            Mesh[i][j][1] = exactPointTuple[1]

# The exact mesh case for the 1x6 mesh for the inviscid channel (P=1)

ExactPointsList = [
    [(-1.594143975246173, 0.187500000000000),(-1.062762650165489, 0.187500000000000),(-0.531381325083476, 0.187500000000000),(0.000000000000000, 0.187500000000000),(0.531381325080713, 0.187500000000000),(1.062762650162713, 0.187500000000000),(1.594143975246173, 0.187500000000000)],
    [(-1.594143975246173, 0.000000000000001),(-1.061329516525977, 0.000000047962662),(-0.528527827678252, 0.001902931798748),(0.000000000000000, 0.062500000000000),(0.528527827674931, 0.001902931798832),(1.061329516522372, 0.000000047962662),(1.594143975246173, 0.000000000000001)],
]

def gaussianBumpExact(i,j):

    """
    Set the exact value for the geometry node points for the P1 mesh
    j = 1 is at the highest point, and j = max is at the minimum y 
    point. However, j input is such that j minimum is at the minimum y point
    """

    jExact = 1-j

    return ExactPointsList[jExact][i]

#================================================================












# the thickness of each element along each coordinate direction
CONST_dx = (CONST_xMax - CONST_xMin)/(CONST_DIMENSION_X)
CONST_dy = (CONST_yMax - CONST_yMin)/(CONST_DIMENSION_Y)


# Note that the 1D GLL points must be placed in ascending sorted order
CONST_GaussLobattoRootsAndCoefficients = {
    2: [[-1,1],[0,0]],
    3: [[-1., 0, 1.], [0.3333333333333333333333, 1.333333333333333333333, 0.3333333333333333333333]],
    4: [[-1, -math.sqrt(5.)/5., math.sqrt(5.)/5., 1],
        [0.1666666666666666666667, 0.833333333333333333333, 0.833333333333333333333, 0.1666666666666666666667]],
    5: [[-1., -math.sqrt(21.)/7., 0.0, math.sqrt(21.)/7., 1.0],[0,0,0,0,0]]
}


# Creating the mesh's file name
CONST_Version = "V4.2"
meshfileName = str(CONST_DIMENSION_X) + "x" + \
    str(CONST_DIMENSION_Y) + "_" + "P" + str(CONST_P)+"_"

if(CONST_Deform == True):
    meshfileName = meshfileName + "Deform_" + str(CONST_CaseType) + "_" + CONST_Version \
        + ".msh"
else:
    meshfileName = meshfileName + "Rect_" + str(CONST_CaseType) + "_" + CONST_Version \
        + ".msh"

CONST_MeshFileName = meshfileName


# This is the data structure for holding all the information about an element.
# It will hold arrays for its solution points and vertices that will be "pointers"
# to a global 2D Mesh array
class Element(object):

    def __init__(self, Mesh, iMin, jMin):
        # iMin and jMin are the indeces for the point in
        # the mesh with the min (x,y) values of the element
        
        # A 2D matrix that holds the vertices of the element. As
        # usual, i=0,j=0 index of this matrix holds the xMin,yMin
        # vertex
        self.Vertices = []
        for i in range(2):
            rowArray = []
            for j in range(2):
                rowArray.append(None)
            self.Vertices.append(rowArray)
    
        for i in range(2):
            for j in range(2):
                # The location of the vertices for the element in
                # 2D mesh array
                iMesh = iMin + i*CONST_P
                jMesh = jMin + j*CONST_P
                
                self.Vertices[i][j] = Mesh[iMesh][jMesh]
        
        #A 2D matrix that holds the the solution points for the element.
        # Same as usual, the i=0, j=0 index holds the minX minY dof.
        self.GeometryNodePoints = []
        for i in range(CONST_P+1):
            rowArray = []
            for j in range(CONST_P+1):
                rowArray.append(None)
            self.GeometryNodePoints.append(rowArray)
                
                #fill the dof data
        for i in range(CONST_P+1):
            for j in range(CONST_P+1):
                iMesh = iMin + i
                jMesh = jMin + j
                
                self.GeometryNodePoints[i][j] = Mesh[iMesh][jMesh]

        self.createGeometryNodePoints()
        
        # Place the vertices/nodes into a 1D array. This is the order in which
        # the nodes will be placed into the triangle->node array in the CPR code.
        # The ordering of the nodes has been made to be consistent with what was
        # done with the old mesh generators.
        self.NodeArray = []

        self.NodeArray.append(self.Vertices[0][0]) # Bottom Left Node
        self.NodeArray.append(self.Vertices[1][0]) # Bottom Right Node
        self.NodeArray.append(self.Vertices[1][1]) # Top Right Node
        self.NodeArray.append(self.Vertices[0][1]) # Top Left Node

        # Store all the geometry node points into a 1D array for all j = 0, then j = 1, etc.
        self.GeometryNodeArray = []
        for j in range(CONST_P+1):
            for i in range(CONST_P+1):
                self.GeometryNodeArray.append(self.GeometryNodePoints[i][j])
        
        #For storing the indeces of the node (in the connectivity file)
        self.NodeArrayConnectivityFileIndeces = []
        

        self.IMEX = 0
        self.Partition = 0 #The partition number for the element
        

    # The method that is used for locating the positions of the solution points
    # (or DOFs) for the element based on the locations of the vertices.
    def createGeometryNodePoints(self):
        
        dxi = 2./CONST_P
        deta = 2./CONST_P

        #Create a 2D tensor product of the 1D GLL points.
        for i in range(CONST_P+1):
            for j in range(CONST_P+1):
                # The (xi,eta) value of the solution point in Computational Domain
                #xi = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][i]
                #eta = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][j]
                
                xi = i*dxi + -1.
                eta = j*deta + -1.

                # The mapped (x,y) value of the point in the Physical Domain
                (x,y) = self.mapParentToPhysicalRectangular(xi,eta)
    
                # Set the (x,y) value of the mapped solution point
                self.GeometryNodePoints[i][j][0] = x
                self.GeometryNodePoints[i][j][1] = y


    # Given a (xi,eta) value on the computational domain, this method
    # will give the (x,y) value of the point on the physical domain for
    # if the element was rectangular.
    def mapParentToPhysicalRectangular(self,xi,eta):
      
        # The width (dx) and length (dy) of the element (when it is in its
        # original rectangular state)
        deltaX = self.Vertices[1][1][0] - self.Vertices[0][0][0]  #xMax - xMin
        deltaY = self.Vertices[1][1][1] - self.Vertices[0][0][1]  #yMax - yMin
        
        xFactor = (xi-(-1.))/2.
        yFactor = (eta-(-1.))/2.
        
        #Mapped location of point on physical domain
        x = self.Vertices[0][0][0] + (deltaX)*(xFactor)
        y = self.Vertices[0][0][1] + (deltaY)*(yFactor)
        
        return (x,y)

    # Given two nodes, get the face of the element. As input, this
    # method takes the pointers to the Mesh's point/node.
    def getFaceIndex(self,node1,node2):
        """
        The convention used in the code for denoting what the index
        of a face should be based on the two nodes (n) it is in between
            
            f=1
        n2 ---- n1
      f=2 |    | f=0
        n3 ---- n0
            f=3
            
        """
        
        faceIndex = -1
        # search through the node array. Break when the first match is found
        # with one of the nodes
        index1 = self.NodeArray.index(node1)
        index2 = self.NodeArray.index(node2)
        
        if (index1 <= index2):
            faceIndex = index1
        else:
            faceIndex = index2
        
        # because the nodes loop, add a condition that if the two points
        # are the last and first of the array, then faceIndex shouldn't be 0,
        # it should be the index of the last point of the array
        
        if((index1 == 0 and index2 == (len(self.NodeArray)-1)) or \
           (index2 == 0 and index1 == (len(self.NodeArray)-1))):
            faceIndex = len(self.NodeArray)-1

        return faceIndex


#Takes as input the Element matrix and plots all the
# elements
def plotElements(ElementObjects, MeshVerticesArray):
    xVector = []
    yVector = []
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            elemObject = ElementObjects[i][j]
            
            for iDof in range(CONST_P+1):
                for jDof in range(CONST_P+1):
                    #loop through all the solution points for the element
                    xVector.append(elemObject.GeometryNodePoints[iDof][jDof][0])
                    yVector.append(elemObject.GeometryNodePoints[iDof][jDof][1])

    plt.scatter(xVector,yVector,s=10,c='b')
    plt.grid()

    # create the lines around all the elements
    for i1 in range(CONST_DIMENSION_X):
        for j1 in range(CONST_DIMENSION_Y):
            xVector = []
            yVector = []
            
            elemObject = ElementObjects[i1][j1]
            
            # bottom edge
            for i in range(CONST_P+1):
                xVector.append(elemObject.GeometryNodePoints[i][0][0])
                yVector.append(elemObject.GeometryNodePoints[i][0][1])
        
            # right edge
            for i in range(CONST_P+1):
                xVector.append(elemObject.GeometryNodePoints[CONST_P][i][0])
                yVector.append(elemObject.GeometryNodePoints[CONST_P][i][1])
            
            # top edge
            for i in range(CONST_P+1):
                xVector.append(elemObject.GeometryNodePoints[CONST_P-i][CONST_P][0])
                yVector.append(elemObject.GeometryNodePoints[CONST_P-i][CONST_P][1])
            
            # left edge
            for i in range(CONST_P+1):
                xVector.append(elemObject.GeometryNodePoints[0][CONST_P-i][0])
                yVector.append(elemObject.GeometryNodePoints[0][CONST_P-i][1])
                
            plt.plot(xVector,yVector)

    # Put annotation boxes next to all the vertices points
    deltaX = (CONST_xMax-CONST_xMin)*0.01
    deltaY = (CONST_yMax-CONST_yMin)*0.01

    placedLabels = []
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            elemObject = ElementObjects[i][j]

            for node in elemObject.NodeArray:
                if node not in placedLabels:
                    placedLabels.append(node)
                    indexVal = MeshVerticesArray.index(node)
                    textVal = str(indexVal)
                    plt.text(node[0]+deltaX, node[1]+deltaY, textVal)

    if CONST_PlotXRange is not None and CONST_PlotYRange is not None:
        plt.gca().set_xlim(CONST_PlotXRange)
        plt.gca().set_ylim(CONST_PlotYRange)

    plt.show(block=True)


# For printing the mesh file that will be read by the DG code. 
def PrintDGMeshFile(MeshNodesArray, MeshVerticesArray, ElementObjectsList, 
                    BoundaryConditionsDict):
    
    file = open(CONST_MeshFileName, "w")
    
    #Print some information in the header of the file
    file.write("Number of grid points: \n")
    file.write(str(len(MeshNodesArray)) + "\n")

    file.write("Number of vertices points: \n")
    file.write(str(len(MeshVerticesArray))+ "\n")
    
    file.write("Number of QUADS: \n")
    file.write(str(len(ElementObjectsList)) + "\n")

    file.write("P: \n")
    file.write(str(CONST_P) + "\n")

    file.write("Vertices Nodes Coordinates: \n")
    for point in MeshVerticesArray:
        file.write("%.15e %.15e \n" %(point[0], point[1]))

    file.write("Connectivity Vertices Nodes: \n")
    for elementObject in ElementObjectsList:
        connectivityString = ""
        for node in elementObject.NodeArray:
            indexValue = MeshVerticesArray.index(node)  # Index starts from 0
            connectivityString = connectivityString + \
                    str(indexValue) + " "

        connectivityString = connectivityString + "\n"
        file.write(connectivityString)

    file.write("Geometry Nodes Coordinates: \n")
    #print all the node coordinates with 15 digits of accuracy:
    for point in MeshNodesArray:
        file.write("%.15e %.15e \n" %(point[0], point[1]))

    file.write("Connectivity Geometry Nodes QUAD: \n")
    
    # Now, go through all the elements and print the information about the geometry node 
    # points. Points will be printed for all j = 0, then j = 1, ... Therefore, to check
    # vertex connectivity, simply access every P+1 point
    for elementObject in ElementObjectsList:
        connectivityString = ""
        for node in elementObject.GeometryNodeArray:
            indexValue = MeshNodesArray.index(node)  # Index starts from 0
            connectivityString = connectivityString + \
                    str(indexValue) + " "

        connectivityString = connectivityString + "\n"
        file.write(connectivityString)


    # Print all the boundary conditions now. The code will need to know
    # what the types of boundary conditions there are in the mesh
    file.write("Boundary Conditions : \n%d \n" % len(BoundaryConditionsDict))

    for k in BoundaryConditionsDict:
        file.write("%s \n%d \n" % (k, len(BoundaryConditionsDict[k])))

        for BCString in BoundaryConditionsDict[k]:
            file.write("%s \n" % BCString)

    file.close()


#For filling the geometry nodes of all the elements into a 1D array
def Fill1DGeometryNodeArray(Mesh, NodeArray1D, MeshNumPointsX, MeshNumPointsY):
    # loop through all the node points and place them into
    # a 1D array. Place the points starting from the bottom left, moving
    # up and then repeat by moving to the right.
    

    for j in range(MeshNumPointsY):
        for i in range(MeshNumPointsX):
            NodeArray1D.append(Mesh[i][j])


#Given a number of processors, this method is in charge of partitioning
#the grid
def PartitionGrid(ElementObjectsMatrix):
    
    """
        
    How the grid is partitioned in the case when the number of processors is 4
    
       2   3
       
       0   1
       
       
    """
    
    # get the number of divisions needed along each coordinate direction of the grid
    numDivisions = int(CONST_NumProcessors**(1./2.))
    
    # Note that this may be an underestimate if there are an odd number
    # of elements along a coordinate direction. To account for this, the
    # last block along the given direction will then have one more element.
    numElementsX = int(CONST_DIMENSION_X/numDivisions)
    numElementsY = int(CONST_DIMENSION_Y/numDivisions)
    
    
    # ni is the ith block in i direction and nj is the jth block in the j direction
    for ni in range(numDivisions):
        for nj in range(numDivisions):
            
            partitionNumber = ni + numDivisions*nj;
            
            # Calculate the number of elements along the i direction and j
            # direction to set their partition number
            
            iElemMin = ni*numElementsX
            jElemMin = nj*numElementsY
            
            iElemMax = iElemMin + numElementsX
            jElemMax = jElemMin + numElementsY
            
            # If we are on the last block along a coordinate direction, set the
            # partition number of all elements up to the edge of the mesh.
            if(ni == numDivisions-1):
                iElemMax = CONST_DIMENSION_X
            if(nj == numDivisions-1):
                jElemMax = CONST_DIMENSION_Y
            
            
            for i in range(iElemMin, iElemMax):
                for j in range(jElemMin, jElemMax):
                    ElementObjectsMatrix[i][j].Partition = partitionNumber



#This method is for computing the periodic Boundary conditions of the elements.
def SetupPeriodicVortexBoundaryConditions(BoundaryConditionsDict, \
                                    ElementObjectsMatrix, ElementObjectsList):
    
    """
    A sample grid
         1
        -----
    2 |      | 4
        -----
         3
    """
    
    PeriodicBoundaryConditionsList = []

    #compute the periodic connections between the bottom and top faces of mesh
    # (i.e. side 1 and 3).
    for i in range(CONST_DIMENSION_X):
        elementTopRow = ElementObjectsMatrix[i][CONST_DIMENSION_Y-1]
        elementTopRowIndex = ElementObjectsList.index(elementTopRow)
        elementBotRow = ElementObjectsMatrix[i][0]
        elementBotRowIndex = ElementObjectsList.index(elementBotRow)
        
        #Get the grid points that are on the "open" side of the element
        
        #get the top node points for elementTopRow
        elementTopRowGP1 =elementTopRow.Vertices[0][1]
        elementTopRowGP2 =elementTopRow.Vertices[1][1]
        
        #get the bottom node points for elementTopRow
        elementBotRowGP1 =elementBotRow.Vertices[0][0] #bottom left grid point
        elementBotRowGP2 =elementBotRow.Vertices[1][0]
        
        #get the face index values for the elements
        elementTopRowFaceIndex = elementTopRow.getFaceIndex(elementTopRowGP1,elementTopRowGP2)
        elementBotRowFaceIndex = elementBotRow.getFaceIndex(elementBotRowGP1, elementBotRowGP2)
        
        #store the data in a tuple
        dataTuple = (elementTopRowIndex, elementBotRowIndex, \
                     elementTopRowFaceIndex, elementBotRowFaceIndex)
                     
        PeriodicBoundaryConditionsList.append(dataTuple)
    
    
    #compute the periodic connections between the left and right faces of the mesh
    # (i.e. side 2 and 4).
    for j in range(CONST_DIMENSION_Y):
        elementLeftCol = ElementObjectsMatrix[0][j]
        elementLeftColIndex = ElementObjectsList.index(elementLeftCol)
        elementRightCol = ElementObjectsMatrix[CONST_DIMENSION_X-1][j]
        elementRightColIndex = ElementObjectsList.index(elementRightCol)
        
        
        #get the Left node points for elementLeftCol
        elementLeftColGP1 =elementLeftCol.Vertices[0][0]
        elementLeftColGP2 =elementLeftCol.Vertices[0][1]
        
        #get the right node points for elementRightCol
        elementRightColGP1 = elementRightCol.Vertices[1][0]
        elementRightColGP2 = elementRightCol.Vertices[1][1]
        
        #get the face indeces for the two elements
        elementLeftColFaceIndex = elementLeftCol.getFaceIndex(elementLeftColGP1, elementLeftColGP2)
        elementRightColFaceIndex = elementRightCol.getFaceIndex(elementRightColGP1, elementRightColGP2)
                                                              
        #store the data in a tuple
        dataTuple = (elementLeftColIndex, elementRightColIndex, elementLeftColFaceIndex, elementRightColFaceIndex)
                                                              
        #print dataTuple
        PeriodicBoundaryConditionsList.append(dataTuple)

    # Convert each tuple of BC information into the string and place all strings into a list that 
    # the dictionary keyword will point to
    PeriodicBoundaryConditionsStringList = []

    for tup in PeriodicBoundaryConditionsList:
        s = ""
        for t in tup:
            s = s + str(t) + " "
        PeriodicBoundaryConditionsStringList.append(s)

    BoundaryConditionsDict["Periodic"] = PeriodicBoundaryConditionsStringList



def Fill1DVerticesNodeArray(Mesh, MeshVerticesArray):

    for i in range(CONST_DIMENSION_X+1):
        for j in range(CONST_DIMENSION_Y+1):
            iNode = i*(CONST_P)
            jNode = j*(CONST_P)

            MeshVerticesArray.append(Mesh[iNode][jNode])


def SetupInviscidChannelBoundaryConditions(BoundaryConditionsDict, ElementObjectsMatrix, 
                                            ElementObjectsList):
    
    """
    Setup the Boundary Conditions for the Inviscid Channel case
    Will have the following boundary conditions:
        - SlipWall at the bottom and top of mesh
        - Total Temperature and Pressure at left of mesh
        - Back Pressure at right of mesh
    """

    SlipWallBC = []

    # Top Surface:
    for i in range(CONST_DIMENSION_X):
        elementTopRow = ElementObjectsMatrix[i][CONST_DIMENSION_Y-1]
        elementTopRowIndex = ElementObjectsList.index(elementTopRow)
        
        #Get the grid points that are on the "open" side of the element
        
        #get the top node points for elementTopRow
        elementTopRowGP1 =elementTopRow.Vertices[0][1]
        elementTopRowGP2 =elementTopRow.Vertices[1][1]
        
        #get the face index values for the elements
        elementTopRowFaceIndex = elementTopRow.getFaceIndex(elementTopRowGP1,elementTopRowGP2)
        
        # Get the string for the BC information:
        bcString = str(elementTopRowIndex) + " " + str(elementTopRowFaceIndex)
                     
        SlipWallBC.append(bcString)

    # Bottom Surface:
    for i in range(CONST_DIMENSION_X):
        elem = ElementObjectsMatrix[i][0]
        elemIndex = ElementObjectsList.index(elem)
        
        #Get the grid points that are on the "open" side of the element
        
        #get the bottom node points for elem
        elemGP1 =elem.Vertices[0][0]
        elemGP2 =elem.Vertices[1][0]
        
        #get the face index values for the elements
        elemFaceIndex = elem.getFaceIndex(elemGP1,elemGP2)
        
        # Get the string for the BC information:
        bcString = str(elemIndex) + " " + str(elemFaceIndex)
                     
        SlipWallBC.append(bcString)


    BoundaryConditionsDict["SlipWall"] = SlipWallBC

    TotalTemperaturePressureBC = []

    # Left Surface:
    for j in range(CONST_DIMENSION_Y):
        elem = ElementObjectsMatrix[0][j]
        elemIndex = ElementObjectsList.index(elem)
        
        #Get the grid points that are on the "open" side of the element
        
        #get the left node points for elem
        elemGP1 =elem.Vertices[0][0]
        elemGP2 =elem.Vertices[0][1]
        
        #get the face index values for the elements
        elemFaceIndex = elem.getFaceIndex(elemGP1,elemGP2)
        
        # Get the string for the BC information:
        bcString = str(elemIndex) + " " + str(elemFaceIndex)
                     
        TotalTemperaturePressureBC.append(bcString)

    BoundaryConditionsDict["TotalTemperaturePressure"] = TotalTemperaturePressureBC

    BackPressureBC = []

    # Right Surface:
    for j in range(CONST_DIMENSION_Y):
        elem = ElementObjectsMatrix[CONST_DIMENSION_X-1][j]
        elemIndex = ElementObjectsList.index(elem)
        
        #Get the grid points that are on the "open" side of the element
        
        #get the right node points for elem
        elemGP1 =elem.Vertices[1][0]
        elemGP2 =elem.Vertices[1][1]
        
        #get the face index values for the elements
        elemFaceIndex = elem.getFaceIndex(elemGP1,elemGP2)
        
        # Get the string for the BC information:
        bcString = str(elemIndex) + " " + str(elemFaceIndex)
                     
        BackPressureBC.append(bcString)

    BoundaryConditionsDict["BackPressure"] = BackPressureBC



# The Main method for the script.

# IMPORTANT: This code works by using implicit references to array elements.
# So, once the Mesh list is created and pointers to values in this list are
# created in the Element objects, they cannot be reassigned to another point
# but only changed. So the = operator cannot be used.
def main():
    
    # First, create the 2D matrix that will hold all the Mesh points.
    # This will include the vertices as well as the dofs
    
    MeshNumPointsX = CONST_DIMENSION_X*(CONST_P) + 1
    MeshNumPointsY = CONST_DIMENSION_Y*(CONST_P) + 1
    
    # Create the 2D Mesh object and initialize all the elements
    # in it to be [0,0]. The mesh object will have at i=0,j=0 the
    # xMin, YMin point and at i=Max, j=Max the xMax, YMax point.
    Mesh = []
    for i in range(MeshNumPointsX):
        rowArray = []
        for j in range(MeshNumPointsY):
            rowArray.append([0,0])
        Mesh.append(rowArray)
    
    for i in range(CONST_DIMENSION_X+1):
        for j in range(CONST_DIMENSION_Y+1):

            VertexXValue = i*CONST_dx + CONST_xMin
            VertexYValue = j*CONST_dy + CONST_yMin
            
            MeshIndexI = i*(CONST_P)
            MeshIndexJ = j*(CONST_P)
            
            Mesh[MeshIndexI][MeshIndexJ][0] = VertexXValue #x value of point
            Mesh[MeshIndexI][MeshIndexJ][1] = VertexYValue #y value of point


    # Create the element objects
    ElementObjectsMatrix = [] #The 2D array that will hold all the element objects
    
    # The i=0, j=0 index will refer to the element in the bottom
    # left of the grid.
    for i in range(CONST_DIMENSION_X):
        rowArray = []
        for j in range(CONST_DIMENSION_Y):
            rowArray.append(None)
        ElementObjectsMatrix.append(rowArray)

    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            # Starting from the bottom left of the grid, create all the
            # element objects
            
            # the index of the min (x,y) value of the vertex for the element in
            # question point on the Mesh element
            iMin = i*(CONST_P)
            jMin = j*(CONST_P)
            
            ElementObjectsMatrix[i][j] = Element(Mesh,iMin,jMin)

    #Fill a 1D array with all the geometry node points
    MeshNodesArray = []
    Fill1DGeometryNodeArray(Mesh, MeshNodesArray, MeshNumPointsX, MeshNumPointsY)

    # Fill a 1D array with the vertices node points
    MeshVerticesArray = []
    Fill1DVerticesNodeArray(Mesh, MeshVerticesArray)

    # Fill a 1D array with all the mesh elements:
    ElementObjectsList = []
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            # Creating pointers in this 1D list to the element objects. This is
            # the order in which the elements will be printed into the connectivity
            # file.
            ElementObjectsList.append(ElementObjectsMatrix[i][j])

    # Dictionary that will hold all the boundary conditions. It will have along
    # with each keyword the list of strings that will be printed for the boundary 
    # condition.
    BoundaryConditionsDict = {}

    # All that will change for the mesh will depend on what case is being run
    if CONST_CaseType == "PeriodicVortex":
        # Set up the periodic BC information for the periodicVortex case

        SetupPeriodicVortexBoundaryConditions(BoundaryConditionsDict, \
                                        ElementObjectsMatrix, ElementObjectsList)


    elif CONST_CaseType == "InviscidChannel":
        SetupInviscidChannelBoundaryConditions(BoundaryConditionsDict, \
                                        ElementObjectsMatrix, ElementObjectsList)
    
    else:
        print "Case Not Available"
        exit(0)
                        
    # Partition the grid properly:
    PartitionGrid(ElementObjectsMatrix)

    # If the mesh needs to be deformed
    if CONST_Deform:
        perturbGrid(Mesh, MeshNumPointsX, MeshNumPointsY)

    # Print all the data for the elements into the connectivity file
    PrintDGMeshFile(MeshNodesArray, 
                    MeshVerticesArray, 
                    ElementObjectsList, 
                    BoundaryConditionsDict)

    # Plot the Mesh
    if CONST_Plot:
        plotElements(ElementObjectsMatrix, MeshVerticesArray)

main()


"""
    
    The sign convention used in the arrays is that increasing i means
    moving to the right and increasing j means moving up in the element object
    matrix, solution points matrix, etc ...
    
    ^
    |
   j|
     ----> i
    
"""




