
import matplotlib.pyplot as plt

inputFile = "input.txt"



# The node points on the computational domain
"""
CONST_PointsMatrix = [[(-1.,1.),(1.,1.)],\
				[(-1.,-1.),(1.,-1.)]]
"""

CONST_PointsMatrix = [[(-1.,-1.),(-1.,0.),(-1.,1.)],\
					  [(0.,-1.),(0.,0.),(0.,1.)],\
					  [(1.,-1.),(1.,0.),(1.,1.)]]

CONST_Dim = 3

# this function computes the lagrange basis at a node
# i. It takes in a list of roots and the index value for
# the root at which the lagrange basis (l_i(x)) is being
# computed
def LagrangeBasis_i(listRoots,index_i,x):
	
	# Compute the numerator of the basis:

	productNumerator = 1.
	# go through each of the roots and take the x value and subtract
	# from it the value of the root. Dont subtract the ith root.
	for i in range(len(listRoots)):
		if (i != index_i):
			productNumerator *= (float(x) - listRoots[i])


	# Compute the denominator of the basis
	productDenominator = 1.
	root_i = listRoots[index_i]
	for i in range(len(listRoots)):
		if (i != index_i):
			productDenominator *= (root_i - listRoots[i])


	solution = float(productNumerator)/float(productDenominator)
	return solution



# The method that is used for computing the shape functions. It will
# then store these functions into a 2D list, where each spot will be the
# function for that corresponding node.
def computeShapeFunctions(ShapeFunctionsMatrix):
	for iNode in range(CONST_Dim):
		for jNode in range(CONST_Dim):
			
			# For this point in the nodes matrix, compute the
			# constant xi lagrange polynomial first
			
			xiNodes = []
			# change the i while keeping j fixed in computational
			# element
			for i in range(CONST_Dim):
				xiNodes.append(CONST_PointsMatrix[i][jNode][0])
			
			etaNodes = []
			# change the j while keeping i fixed in computational
			# element
			for j in range(CONST_Dim):
				etaNodes.append(CONST_PointsMatrix[iNode][j][1])
			
			# Compute the lagrange bases:
			LXi = lambda xi, iNode=iNode: (LagrangeBasis_i(xiNodes, iNode, xi))
			LEta = lambda eta, jNode=jNode: (LagrangeBasis_i(etaNodes, jNode, eta))
			
			# Compute the shape function at the node:
			M = lambda xi,eta, LXi=LXi,LEta = LEta: (LXi(xi)*LEta(eta))
	
			# store the shape function at the node location
			ShapeFunctionsMatrix[iNode][jNode] = M

"""

plt.scatter(elemGeoNodesX, elemGeoNodesY, c = 'b')
plt.scatter(elemSolNodesX, elemSolNodesY, c = 'r')
plt.grid()
plt.show(block=True)

"""

def main():


	with open(inputFile, 'r') as fp:
		l = fp.readline()
		numElems, numNodes = int(l.split()[0]), int(l.split()[1])

		elemGeoNodesX = []
		elemGeoNodesY = []

		elemSolNodesX = []
		elemSolNodesY = []

		individualElemNodes = []
		for elem in range(numElems):

			indElemNodes = []

			for node in range(numNodes):
				l = fp.readline()
				x,y = float(l.split()[-2]), float(l.split()[-1])

				indElemNodes.append((x,y))
				elemGeoNodesY.append(y)
				elemGeoNodesX.append(x)

			for node in range(numNodes):
				l = fp.readline()
				x,y = float(l.split()[-2]), float(l.split()[-1])

				elemSolNodesY.append(y)
				elemSolNodesX.append(x)

			individualElemNodes.append(indElemNodes)

	"""
	plt.scatter(elemGeoNodesX, elemGeoNodesY, c = 'b')
	plt.scatter(elemSolNodesX, elemSolNodesY, c = 'r')
	plt.grid()
	plt.show(block=True)
	"""

	"""
	GL Nodes:
	x1 = -0.774597
	x2 = 0
	x3 = 0.774597
	"""

	# Check the first element's geometry nodes and get the solution node positions
	# First, place the nodes into a 2D array in the correct order

	elem1GeoNodes = []
	for i in range(3):
		row = []
		for j in range(3):
			row.append(None)
		elem1GeoNodes.append(row)

	# Ordering of the nodes for some reason
	jOrder = [0,2,1]
	iOrder = [0,2,1]

	n = 0
	for j in range(3):
		for i in range(3):
			elem1GeoNodes[iOrder[i]][jOrder[j]] = individualElemNodes[0][n]
			n = n+1
		


	for j in range(3):
		for i in range(3):
			print "i,j = %d, %d -> %f, %f "%(i,j,elem1GeoNodes[i][j][0], elem1GeoNodes[i][j][1])


	# Create the matrix for the Shape functions
	ShapeFunctionsMatrix = []
	for i in range(CONST_Dim):
		rowArray = []
		for j in range(CONST_Dim):
			rowArray.append(None)
		ShapeFunctionsMatrix.append(rowArray)
			
	computeShapeFunctions(ShapeFunctionsMatrix)

	# Interpolate now the x,y position of the solution nodes
	xiGLL =  -0.774597
	etaGLL = -0.0

	xVal = 0
	yVal = 0

	for i in range(CONST_Dim):
		for j in range(CONST_Dim):
			xVal = xVal + ShapeFunctionsMatrix[i][j](xiGLL, etaGLL)*elem1GeoNodes[i][j][0]
			yVal = yVal + ShapeFunctionsMatrix[i][j](xiGLL, etaGLL)*elem1GeoNodes[i][j][1]
			
	print "xGLL, yGLL: %f %f "%(xVal, yVal)


main()








