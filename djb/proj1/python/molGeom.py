'''
    Molecular Geometry Analysis

    Author: Daniel Burrill
    Date: January 4, 2016

    Description: Perform geometry analysis as per Crawdad programming project 1
'''

# Imports
import numpy as np

# Classes
class atom():
    ## Constructor
    def __init__(self,atomicNum,xPos,yPos,zPos):
        self.atomicNumber = atomicNum
        self.xPos = xPos
        self.yPos = yPos
        self.zPos = zPos

class molecule():
    ## Constructor
    def __init__(self,inFileName):
        # Initialize atom list
        self.atomList = []

        # Load file
        self.atomicNum, self.xPos, self.yPos, self.zPos = np.loadtxt(inFileName,skiprows=1,unpack=True)

        # Assemble atoms
        for index,aN in enumerate(self.atomicNum):
            self.atomList.append(atom(aN,self.xPos[index],self.yPos[index],self.zPos[index]))

    def bondLength(self,atom1,atom2):
        '''
            Calculate bond length a1-a2

            INPUT
                atom1: (atom class)
                atom2: (atom class)

            OUTPUT
                Returns (float) in distance units
        '''

        # Define positions
        r1 = [atom1.xPos,atom1.yPos,atom1.zPos]
        r2 = [atom2.xPos,atom2.yPos,atom2.zPos]

        return np.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)

    def bondAngle(self,atom1,atom2,atom3):
        '''
            Calculate bond angle a1-a2-a3

            INPUT
                atom1: (atom class)
                atom2: (atom class) Central atom
                atom3: (atom class)

            OUTPUT
                Returns (float) in radians
        '''

        # Define positions
        r1 = [atom1.xPos,atom1.yPos,atom1.zPos]
        r2 = [atom2.xPos,atom2.yPos,atom2.zPos]
        r3 = [atom3.xPos,atom3.yPos,atom3.zPos]

        # Calculate bond lengths
        r12 = bondLength(atom1,atom2)
        r23 = bondLength(atom2,atom3)

        # Catch the case where two atoms are the same
        if ((r12 <= 0.00000001) or (r23 <= 0.00000001)):
            return 0.0

        # Calculate unit vectors
        e1 = [(r2[0]-r1[0])/r12,(r2[1]-r1[1])/r12,(r2[2]-r1[2])/r12]
        e2 = [(r3[0]-r2[0])/r23,(r3[1]-r2[1])/r23,(r3[2]-r2[2])/r23]

        return np.acos(np.dot(e1,e2))

    def OOPAngle(atom1,atom2,atom3,atom4):
        '''
            Calculate out-of-plane angle a1->a3-a2-a4

            INPUT
                atom1: (atom class) Out-of-plane atom
                atom2: (atom class)
                atom3: (atom class) Central atom
                atom4: (atom class)

            OUTPUT
                Returns (float) in radians
        '''

        # Define positions
        r1 = [atom1.xPos,atom1.yPos,atom1.zPos]
        r2 = [atom2.xPos,atom2.yPos,atom2.zPos]
        r3 = [atom3.xPos,atom3.yPos,atom3.zPos]
        r4 = [atom4.xPos,atom4.yPos,atom4.zPos]

        # Calculate bond lengths
        r13 = bondLength(atom1,atom3)
        r23 = bondLength(atom2,atom3)
        r34 = bondLength(atom4,atom3)

        # Catch the case where two atoms are the same
        if ((r23 <= 0.00000001) or (r34 <= 0.00000001)):
            return 0.0

        # Calculate unit vectors
        e23 = [(r3[0]-r2[0])/r23,(r3[1]-r2[1])/r23,(r3[2]-r2[2])/r23]
        e34 = [(r4[0]-r3[0])/r34,(r4[1]-r3[1])/r34,(r4[2]-r3[2])/r34]

        crossProd = np.cross(e23,e34)

# Main
if (__name__ == '__main__'):
    molClass = molecule("acetaldehyde.dat")

    # Print out bond lengths
    print "Index1\tIndex2\tBondLength"
    for index,atom in enumerate(molClass.atomList):
        for index2 in range(index,len(molClass.atomList)):
            # Filter out cases where the atom indices are the same
            if (index != index2):
                print str(index) + '\t' + str(index2) + '\t' + str(molClass.bondLength(atom,molClass.atomList[index2]))
