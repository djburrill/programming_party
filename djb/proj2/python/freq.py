'''
    Molecular Geometry Analysis

    Author: Daniel Burrill
    Date: January 7, 2016

    Description: Calculate vibrational frequency of h2o
'''

# Imports
import numpy as np
import periodic as per
import scipy.constants as sc

# Classes
class atom():
    ## Constructor
    def __init__(self,atomicNum,xPos,yPos,zPos):
        # Atomic number
        self.atomicNumber = atomicNum

        # Position
        self.xPos = xPos
        self.yPos = yPos
        self.zPos = zPos

        # Average mass
        self.mass = per.element(self.atomicNumber).mass

class molecule():
    ## Constructor
    def __init__(self,inFileName):
        # Set tolerance
        self.tol = 0.000001

        # Initialize atom list
        self.atomList = []

        # Load file
        self.atomicNum, self.xPos, self.yPos, self.zPos = np.loadtxt(inFileName,skiprows=1,unpack=True)

        # Assemble atoms
        for index,aN in enumerate(self.atomicNum):
            self.atomList.append(atom(aN,self.xPos[index],self.yPos[index],self.zPos[index]))

    def calcFreq(self,inFileName):
        '''
            Calculate vibrational frequency.

            INPUT
                inFileName: (string) Name of file containing Hessian data

            OUTPUT
                Returns (float) list [w1,w2,...] of frequencies.
        '''

        # Load Hessian
        hess = np.ones((9,9))
        inHess = np.matrix(np.loadtxt(inFileName,skiprows=1))

        # Reshape Hessian
        for index1 in range(inHess.shape[0]):
            for index2 in range(inHess.shape[1]):
                hess[(index2+index1)%9,]

        # Calculate mass weighted Hessian
        print hessian


# Main
if (__name__ == '__main__'):
    molClass = molecule("h2o.dat")

    # Print hessian
    molClass.calcFreq("h2o_hessian.dat")
