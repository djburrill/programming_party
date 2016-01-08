'''
    Molecular Geometry Analysis

    Author: Daniel Burrill
    Date: January 4, 2016

    Description: Perform geometry analysis as per Crawdad programming project 1
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

    def translate(self,dx,dy,dz):
        '''
            Translate the position of the atom to [x+dx,y+dy,z+dz].

            INPUT
                dx: (float) Translation in x
                dy: (float) Translation in y
                dz: (float) Translation in z

            OUTPUT
                No output. Updates atom position.
        '''

        self.xPos += dx
        self.yPos += dy
        self.zPos += dz

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

    def unitVec(self,atom1,atom2):
        '''
            Calculate the unit vector from atom1 to atom2.

            INPUT
                atom1: (atom class)
                atom2: (atom class)

            OUTPUT
                Returns list (float).
        '''

        # Define atom positions
        r1 = [atom1.xPos,atom1.yPos,atom1.zPos]
        r2 = [atom2.xPos,atom2.yPos,atom2.zPos]

        # Calculate bond lengths
        r12 = self.bondLength(atom1,atom2)

        # Catch the case where two atoms are the same
        if (r12 <= self.tol):
            return [0.0,0.0,0.0]

        # Calculate unit vectors
        return [(r2[0]-r1[0])/r12,(r2[1]-r1[1])/r12,(r2[2]-r1[2])/r12]


    def bondLength(self,atom1,atom2):
        '''
            Calculate bond length a1-a2.

            INPUT
                atom1: (atom class)
                atom2: (atom class)

            OUTPUT
                Returns (float) in distance units.
        '''

        # Define positions
        r1 = [atom1.xPos,atom1.yPos,atom1.zPos]
        r2 = [atom2.xPos,atom2.yPos,atom2.zPos]

        return np.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)

    def bondAngle(self,atom1,atom2,atom3):
        '''
            Calculate bond angle a1-a2-a3.

            INPUT
                atom1: (atom class)
                atom2: (atom class) Central atom
                atom3: (atom class)

            OUTPUT
                Returns (float) in radians.
        '''

        # Calculate unit vectors
        e1 = self.unitVec(atom2,atom1)
        e2 = self.unitVec(atom2,atom3)

        # Calculate argument of arccos()
        sol = np.dot(e1,e2)

        # Check that returned value does not produce infinities
        if (sol < -1):
            return np.arccos(-1)
        elif (sol > 1):
            return np.arccos(1)
        else:
            return np.arccos(sol)

    def OOPAngle(self,atom1,atom2,atom3,atom4):
        '''
            Calculate out-of-plane angle a1->a3-a2-a4.

            INPUT
                atom1: (atom class) Out-of-plane atom
                atom2: (atom class)
                atom3: (atom class) Central atom
                atom4: (atom class)

            OUTPUT
                Returns (float) in radians.
        '''

        # Calculate unit vectors
        e32 = self.unitVec(atom3,atom2)
        e34 = self.unitVec(atom3,atom4)
        e31 = self.unitVec(atom3,atom1)

        crossProd = np.cross(e32,e34)
        leftVec = crossProd/np.sin(self.bondAngle(atom2,atom3,atom4))

        # Calculate argument of arcsin()
        sol = np.dot(leftVec,e31)

        # Check that arcsin does not produce infinities
        if (sol < -1):
            return np.arcsin(-1)
        elif (sol > 1):
            return np.arcsin(1)
        else:
            return np.arcsin(sol)

    def torsionAngle(self,atom1,atom2,atom3,atom4):
        '''
            Calculate torsion/dihedral angle a1-a2-a3-a4.

            INPUT
                atom1: (atom class)
                atom2: (atom class)
                atom3: (atom class)
                atom4: (atom class)

            OUTPUT
                Returns (float) in radians.
        '''

        # Calculate unit vectors
        e23 = self.unitVec(atom2,atom3)
        e34 = self.unitVec(atom3,atom4)
        e12 = self.unitVec(atom1,atom2)

        # Calculate cross products
        crossProd1 = np.cross(e12,e23)
        crossProd2 = np.cross(e23,e34)

        # Calculate bond angles
        bondAngle1 = self.bondAngle(atom1,atom2,atom3)
        bondAngle2 = self.bondAngle(atom2,atom3,atom4)

        # Compute argument of arccos (Check for infinity)
        if ((bondAngle1 != 0) and (bondAngle2 != 0.0)):
            sol = np.dot(crossProd1,crossProd2)/(np.sin(bondAngle1)*np.sin(bondAngle2))
        else:
            return 0.0

        # Check that returned value does not produce infinities
        if (sol < -1):
            return np.arccos(-1)
        elif (sol > 1):
            return np.arccos(1)
        else:
            return np.arccos(sol)

    def calcCOM(self):
        '''
            Calculate center of mass position.

            INPUT

            OUTPUT
                Returns (float) list [x,y,z] of position.
        '''

        COMpos = [0.0,0.0,0.0]
        totalMass = 0.0

        for atom in self.atomList:
            COMpos[0] += atom.mass*atom.xPos
            COMpos[1] += atom.mass*atom.yPos
            COMpos[2] += atom.mass*atom.zPos
            totalMass += atom.mass

        COMpos[0] /= totalMass
        COMpos[1] /= totalMass
        COMpos[2] /= totalMass

        return COMpos

    def translate(self,dx,dy,dz):
        '''
            Translate the position of the molecule to [x+dx,y+dy,z+dz]. Translates all constituent atoms.

            INPUT
                dx: (float) Translation in x
                dy: (float) Translation in y
                dz: (float) Translation in z

            OUTPUT
                No output. Updates atom positions.
        '''

        for atom in self.atomList:
            atom.translate(dx,dy,dz)

    def inertia(self):
        '''
            Calculate principal moments of interia

            INPUT

            OUTPUT
                Returns list (float) of the PMOI [Ia,Ib,Ic]
        '''

        # Initialize inertia matrix
        iMat = np.ones((3,3))

        # Calculate elements
        for atom in self.atomList:
            iMat[0,0] += atom.mass*(atom.yPos**2 + atom.zPos**2)
            iMat[1,1] += atom.mass*(atom.xPos**2 + atom.zPos**2)
            iMat[2,2] += atom.mass*(atom.xPos**2 + atom.yPos**2)

            iMat[0,1] += atom.mass*atom.xPos*atom.yPos
            iMat[0,2] += atom.mass*atom.xPos*atom.zPos
            iMat[1,2] += atom.mass*atom.yPos*atom.zPos

        # Apply symmetry
        iMat[1,0] += iMat[0,1]
        iMat[2,0] += iMat[0,2]
        iMat[2,1] += iMat[1,2]

        # Calculate eigenvalues of iMat & sort
        eigs = np.linalg.eigvals(iMat)
        eigs = np.sort(eigs)

        return [eigs[0],eigs[1],eigs[2]]

# Main
if (__name__ == '__main__'):
    molClass = molecule("acetaldehyde.dat")

    # Print out bond lengths
    print "Index1\tIndex2\tBondLength"
    for index,atom in enumerate(molClass.atomList):
        for index2 in range(index,len(molClass.atomList)):
            # Filter out cases where the atom indices are the same
            if (index != index2):
                ai = atom
                aj = molClass.atomList[index2]
                print str(index) + '\t' + str(index2) + '\t' + str(molClass.bondLength(ai,aj))

    # Print out bond angles
    print "Index1\tIndex2\tIndex3\tBondAngle"
    for index,atom in enumerate(molClass.atomList):
        for index2 in range(index):
            for index3 in range(index2):
                ai = atom
                aj = molClass.atomList[index2]
                ak = molClass.atomList[index3]

                # Filter out cases where the atom indices are the same
                if ((index != index2) and (index2 != index3) and (index != index3)):
                    # Only calculate bonded atoms
                    if ((molClass.bondLength(ai,aj) < 4.0) and (molClass.bondLength(aj,ak) < 4.0)):
                        print str(index) + '\t' + str(index2) + '\t' + str(index3) + '\t' + str(np.degrees(molClass.bondAngle(ai,aj,ak)))

    # Print OOP angles
    print "Index1\tIndex2\tIndex3\tIndex4\tOOPAngle"
    for index,atom in enumerate(molClass.atomList):
        for index2 in range(len(molClass.atomList)):
            for index3 in range(len(molClass.atomList)):
                for index4 in range(index3):
                    ai = atom
                    aj = molClass.atomList[index3]
                    ak = molClass.atomList[index2]
                    al = molClass.atomList[index4]

                    # Filter out cases where the atom indices are the same and use bonded set of atoms only
                    if ((index != index2) and (index != index3) and (index != index4) and (index2 != index3) and (index2 != index4) and (molClass.bondLength(ai,ak) < 4.0) and (molClass.bondLength(ak,aj) < 4.0) and (molClass.bondLength(ak,al) < 4.0)):
                        print str(index) + '\t' + str(index3) + '\t' + str(index2) + '\t' + str(index4) + '\t' + str(np.degrees(molClass.OOPAngle(ai,aj,ak,al)))

    # Print Torsion angles
    print "Index1\tIndex2\tIndex3\tIndex4\tTorsionAngle"
    for index,atom in enumerate(molClass.atomList):
        for index2 in range(index):
            for index3 in range(index2):
                for index4 in range(index3):
                    ai = atom
                    aj = molClass.atomList[index2]
                    ak = molClass.atomList[index3]
                    al = molClass.atomList[index4]

                    # Filter out cases where the atom indices are the same and use bonded set of atoms only
                    if ((molClass.bondLength(ai,aj) < 4.0) and (molClass.bondLength(aj,ak) < 4.0) and (molClass.bondLength(ak,al) < 4.0)):
                        print str(index) + '\t' + str(index2) + '\t' + str(index3) + '\t' + str(index4) + '\t' + str(np.degrees(molClass.torsionAngle(ai,aj,ak,al)))

    # Print center of mass and translate molecule
    COM = molClass.calcCOM()
    molClass.translate(-COM[0],-COM[1],-COM[2])
    print "Center of mass: " + str(COM[0]) + '\t' + str(COM[1]) + '\t' + str(COM[2])

    # Calculate principal moments of inertia
    PMOI = molClass.inertia()
    # Note that the values are slightly different! Which one is correct?!
    print PMOI

    # Calculate rotational constants (cm^-1)
    A = (sc.h/sc.c)*(1.0/(8.0*np.pi*np.pi)/(0.5291772e-10*0.5291772e-10*1.6605389e-27))
    B = (sc.h/sc.c)*(1.0/(8.0*np.pi*np.pi)/(0.5291772e-10*0.5291772e-10*1.6605389e-27))
    C = (sc.h/sc.c)*(1.0/(8.0*np.pi*np.pi)/(0.5291772e-10*0.5291772e-10*1.6605389e-27))
    print "Rotational Constants: A=" + str(A/PMOI[0]/100) + '\tB=' + str(B/PMOI[1]/100) + '\tC=' + str(C/PMOI[2]/100)
