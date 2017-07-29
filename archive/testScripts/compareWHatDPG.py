
import math

DPGFile = "DPGInput.txt"
MatlabFile = "MatlabInput.txt"

# For comparing the points
CONST_FloatComparison_Tolerance = 1E-3



# Method for reading in the files and loading the point data
def processFiles(fileName):

	# Fill the points into a list with the form:
	#	l = [ [x,y,[wVec]], ... ]

	retVec = []
	with open(fileName, "r") as fp:

		while(True):
			l1 = fp.readline()
			if l1 == "":
				break
			l2 = fp.readline()

			# Get the points:
			l1 = l1.rstrip("\n").split()
			xVal = float(l1[1].rstrip(","))
			yVal = float(l1[3])

			l2 = l2.rstrip("\n").rstrip(" ")
			l2 = l2.split(" = ")[1]
			l2 = l2.rstrip("]")
			l2 = l2.lstrip("[")
			l2 = l2.split(",")

			wVec = []
			for s in l2:
				wVec.append(float(s))
			
			retVec.append([xVal, yVal, wVec])


	return retVec


def compareOutputes(l1, l2):
	# return a list with the x,y value of the point and the L2 norm of the 
	# difference

	retList = []

	for a in l1:
		for b in l2:
			if (abs(b[0]-a[0])<CONST_FloatComparison_Tolerance and
				abs(b[1]-a[1])<CONST_FloatComparison_Tolerance):

				L2 = 0
				for i in range(4):
					L2 = L2 + (a[2][i]-b[2][i])**2

				L2 = math.sqrt(L2)

				retList.append([a[0], a[1], L2])

	return retList


def main():
	matlabPoints = processFiles(MatlabFile)
	DPGPoints = processFiles(DPGFile)

	diffVec = compareOutputes(matlabPoints, DPGPoints)

	diffVec.sort(key=lambda x: x[2])

	for l in diffVec:
		print l


main()






