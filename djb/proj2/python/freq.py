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
                Returns [float] list [e1,e2,...] of eigenvalues (hartree/amu-bohr^2).
        '''

        # Load Hessian
        numAtoms = len(self.atomList)
        hessSize = 3*numAtoms
        hess = np.ones((hessSize,hessSize))
        inHess = np.matrix(np.loadtxt(inFileName,skiprows=1))

        # Reshape Hessian
        for index1 in range(inHess.shape[0]):
            for index2 in range(inHess.shape[1]):
                hess[index1/numAtoms,(index2+3*index1)%hessSize] = inHess[index1,index2]

        # Calculate mass weighted Hessian
        for row in range(hess.shape[0]):
            # Set mass for ith atom
            mi = self.atomList[row/3].mass

            for column in range(hess.shape[1]):
                # Set mass for jth atom
                mj = self.atomList[column/3].mass

                hess[row,column] /= np.sqrt(mi*mj)

        # Print mass weighted Hessian
        np.set_printoptions(linewidth=200,suppress=True)
        print 'Mass weighted Hessian'
        print hess
        print ''

        # Calculate eigenvalues
        eigs = np.linalg.eigvals(hess)
        eigs = np.sort(eigs)

        # Print eigenvalues
        print 'Eigenvalues'
        print eigs
        print ''

        return eigs

        '''
        Note that the actual vibrational frequencies are not calculated because (I think) all interesting programming has been finished.
        '''
# Main
if (__name__ == '__main__'):
    # H2O
    print 'H2O'
    molClass_h2o = molecule("h2o.dat")
    molClass_h2o.calcFreq("h2o_hessian.dat")

    # Benzene
    print 'Benzene'
    molClass_benzene = molecule("benzene.dat")
    molClass_benzene.calcFreq("benzene_hessian.dat")

    # 3-chloro-1-butene
    print '3-chloro-1-butene'
    molClass_3chloro1butene = molecule("3-chloro-1-butene.dat")
    molClass_3chloro1butene.calcFreq("3-chloro-1-butene_hessian.dat")
