#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
perturb_angle.py

Perturb geometry around an improper dihedral angle. The options are:
 1. Change valence angle without changing improper
 2. Change improper angle without changing valence angles

For use in generating geometries to parameterize valence and improper angles for an oemol

By: Victoria Lim and Jessica Maat

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import math
import sys
from openeye import oeomega
from openeye import oechem
from calc_improper import *
#from off_nitrogens.calc_improper import *


#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    Parameters
    ----------
    axis : tuple, list, or numpy array
        vector normal to plane in which rotation is performed
    theta : float
        angle of rotation amount in degrees

    Returns
    -------
    numpy array
        three-dimension rotation matrix

    Reference
    ---------
    https://tinyurl.com/yam5hu78
    """
    axis = np.asarray(axis)
    theta = np.radians(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])





def perturb_valence(atom0, atom1, atom2, atom3, theta, verbose=False):
    """
    Rotate atom3 in the plane of atom1, atom2, and atom3 by theta degrees.
    This involves a counterclockwise rotation relative to the normal vector.
    I.e., if normal vector is (0,0,1) the rotation will be counterclockwise
    looking from +z onto the xy plane. But if the normal vector is (0,0,-1)
    the rotation will be clockwise from the same viewpoint.

    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates TO BE MOVED
    theta : float
        how many degrees by which to rotate

    Returns
    -------
    atom* : numpy arrays
        the coordinates in the same order as given in parameters
    rot_mat : numpy array
        three-dimension rotation matrix, returned so that the same
        matrix can be applied to any atoms attached to atom3.
    """

    # outer atoms are atom1 atom2 atom3. get normal vector to that plane.
    v1 = np.asarray(atom2)-np.asarray(atom1)
    v2 = np.asarray(atom2) - np.asarray(atom3)
    w2 = np.cross(v1, v2)
    # calculate rotation matrix
    rot_mat = rotation_matrix(w2, theta)
    length = np.linalg.norm(atom0-atom3)
    atom3_rot = np.dot(rot_mat, atom3)
    new_length = np.linalg.norm(atom0-atom3_rot)

    return atom0, atom1, atom2, atom3_rot






def perturb_improper(atom0, atom1, atom2, atom3, theta, verbose=False):
    """
    Rotate atom3 out of the plane of atom1, atom2  by theta degrees while keeping valence angle the same.

    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates TO BE MOVED
    theta : float
        how many degrees by which to move atom upwards

    Returns
    -------
    atom* : numpy arrays
        the coordinates in the same order as given in parameters
    rot_mat : numpy array
        three-dimension rotation matrix, returned so that the same
        matrix can be applied to any atoms attached to atom3.

    """

    # get the vector between atom2 and atom1 and then rotate atom3 about that plane by angle theta.
    v1 = atom2-atom1
    # calculate rotation matrix
    rot_mat = rotation_matrix(v1, theta)
    atom3_rot = np.dot(rot_mat, atom3)

    return atom0, atom1, atom2, atom3_rot


def input2mol(xyz):
    """
    d e s c r i p t i o n :
    Takes an xyz file dictionary  and converts the coordinates to an eomol.
    Function returns a dictionary of oemols.

    p a r a m e t e r s :
    xyz: dictionary of xyz files


    """
    ifs = oechem.oemolistream()
    for key, value in xyz.items():
        ifs.open(value)
        for mol in ifs.GetOEGraphMols():
            xyz[key] = oechem.OEGraphMol(mol)
    return xyz

"""
ifs = oechem.oemolistream()
ofs = oechem.oemolostream()

if ifs.open("molClass_pyrnit_molecule_1.xyz"):
        if ofs.open("output.mol2"):
                    for mol in ifs.GetOEGraphMols():
                                    oechem.OEWriteMolecule(ofs, mol)
"""

def smileslist2mol(smilesList):
    """
    This function creates a dictionary of oemol with keys that are
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    oechem.OEAddExplicitHydrogens(mol)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega(mol)



def oemol_perturb(molList, angle_type, molClass, pertRange, pertIncr):
    """
    The function takes in a list of smiles strings, converts the smiles strings into OpenEye
    OEMols and then performs either a valence or improper perturbation to the molecules geometry
    by a specified angle "theta". The function creates an output .sdf file which contains all the
    molecules from the molList and a SD tag that contains the the indices of the atoms that
    are in the improper molecule center. The first index is the center of the improper and the last
    index indicates the atom that was perturbed in the geometry change. The function utilizes
    find_improper_angles from the calc_improper.py script to identify the improper locations in the molecule.

    The code will iterate through each trivalent center, and perturb the centers individually in the output
    .sdf file. The code will also generate 3 perturbed geometries for each nitrogen center where each of the
    individual constiutuents will be perturbed individually.


    Parameters
    ---------
    molList : Dictionary of oemol objects
    angle_tyle:
        True: Boolean
            True = Improper perturbation
            False = Valence perturbation
            Ex., True
    molClass : String
        Type of molecules
        Ex., pyrnitrogens
    PertRang : Int
        The range of degrees the attached atom will be perturbed by
    pertIncr: Int
        The increment of degrees that the molecule with be perturbed by

    Returns
    --------
    .sdf file for each molecule in smiles string that perturbs the improper or valence by increments specified in specified
    range with tags that include the frozen indices of the molecule

    """

    #Set perturbation type
    if angle_type == True:
        angle_type = perturb_improper
        perturbation = "improper"
    else:
        angle_type = perturb_valence
        perturbation = "valence"
    #convert molList into OEMols
    for key, mol in molList.items():
        impCent = find_improper_angles(mol)
        centcount = 0
        constcount = 0
        for center in impCent[1]:
            for constituent in center:
                constcount +=1
                #determine which improper angle on the atom is of interest and store the coordinates of the improper in  various variables
                #cmol is the parent mol which we will add conformers to
                cmol = oechem.OEMol(mol)
                mol_pert = str(constituent)
                center_atom = cmol.GetAtom(oechem.OEHasAtomIdx(center[0]))
                if constituent == center[0]:
                    continue
                move_atom = cmol.GetAtom(oechem.OEHasAtomIdx(constituent))
                center_coord = np.array(cmol.GetCoords(center_atom))
                move_coord = np.array(cmol.GetCoords(move_atom))


                #center_coord = center[0]
                #move_atom = constituent
                other_coords = list()

                #if constituent != center_coord:
                #    if constituent != move_atom:
                #        other_coords.append(constituent)


                for neighbor in center_atom.GetAtoms():
                    if neighbor.GetIdx() !=  constituent:
                        new_coord = np.array(cmol.GetCoords(neighbor))
                        other_coords.append(new_coord)

                # rotate atom by desired increment, and write out each perturbation to the .sdf file
                # DEBUGGING
                theta = 0
                print(pertRange)

                oemol_list = [cmol]

                #in the list adding tags for the indices around the nitrogen center, improver or valence move, and angle of perturbation
                ofile = oechem.oemolostream(str(molClass) + "_" + str(key) + "_constituent_" +  mol_pert +"_" + perturbation +'.sdf')
                count = 1
                theta = (360 - (pertRange/2))
                print("starting theta:" + str(theta))
                maxRange = ((pertRange/2)+360)
                print("Max range:" + str(maxRange))
                pertTheta = -(pertRange/2)
                while theta < maxRange:
                    #move in direction of theta
                    atom0, atom1, atom2, atom3_rot = angle_type(center_coord, other_coords[0], other_coords[1], move_coord, theta)
                    move_mol = oechem.OEMol(oemol_list[0])
                    move_mol.SetCoords(move_atom, oechem.OEFloatArray(atom3_rot))
                    move_mol.SetTitle(str(molClass) + "_" + str(key) + "_constituent_" +  mol_pert +"_" + perturbation)
                    oechem.OEAddSDData(move_mol, "Index list of atoms to freeze", str(center))
                    oechem.OEAddSDData(move_mol, "Angle OEMol is perturbed by", str(pertTheta))
                    oechem.OEAddSDData(move_mol, "Type of perturbation (improper or valence)", perturbation)
                    oemol_list.append(move_mol)
                    oechem.OEWriteConstMolecule(ofile, move_mol)

                    #seperate .mol2 for each perturbation
                    mfile = oechem.oemolostream(str(molClass) + "_" + str(key) + "_constituent_" +  mol_pert + "_" + perturbation + '_angle_'+ str(pertTheta) +'.mol2')
                    oechem.OEWriteConstMolecule(mfile, move_mol)
                    mfile.close()



                    #count for while loop
                    theta += pertIncr
                    count += 1
                    pertTheta += pertIncr


                ofile.close()
                print("end of loop")

    return




