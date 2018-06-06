#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
perturb_angle.py

Perturb geometry around an improper dihedral angle. The options are:
 1. Change valence angle without changing improper
 2. Change improper angle without changing valence angles

For use in generating geometries to parameterize valence and improper angles

By: Victoria Lim and Jessica Maat

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import math
import sys
from openeye import oeomega
from oeommtools.utils import openmmTop_to_oemol
from openeye import oechem

from calc_improper import *

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
    v1 = atom2-atom1
    v2 = atom2 - atom3
    w2 = np.cross(v1, v2)
    # calculate rotation matrix
    rot_mat = rotation_matrix(w2, theta)
    atom3_rot = np.dot(rot_mat, atom3)
    # print details of geometry
    if verbose:
        print("\n>>> Perturbing valence while maintaining improper angle...")
        # atom being moved is atom3.
        print("\natom0 original coords:\t",atom0)
        print("atom1 original coords:\t",atom1)
        print("atom2 original coords:\t",atom2)
        print("atom3 original coords:\t",atom3)
        print("atom3 {:.2f} deg rotated: {}".format(60.,atom3_rot))

        # check improper but make sure the central is first and moved is last
        print("\nImproper angle, before: ", calc_improper_angle(atom0, atom1, atom2, atom3))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom2, atom3, atom1))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom3, atom1, atom2))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom1, atom2, atom3_rot))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom2, atom3_rot, atom1))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom3_rot, atom1, atom2))

        ## check valences
        print("\nvalence angle 1, before: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, before: ", calc_valence_angle(atom0, atom1, atom3))
        print("valence angle 3, before: ", calc_valence_angle(atom0, atom2, atom3))
        print("valence angle 1, after: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, after: ", calc_valence_angle(atom0, atom1, atom3_rot))
        print("valence angle 3, after: ", calc_valence_angle(atom0, atom2, atom3_rot))
        print()

    return atom0, atom1, atom2, atom3_rot, rot_mat

def oemol_perturb(mol, central_atom, outer_atom, angle_type,  theta):
    """
    From an OpenEye OEMol, specify the improper angle and the specific atom of
    that improper that should be perturbed. The improper angles are obtained
    from the find_improper_angles function in the calc_improper script.

    Parameters
    ---------
    mol : OpenEye OEMol
        molecule from which to generate perturbed geometry
    central_atom : indices
        atom indice in the mol which is central to the improper of interest
        Ex., "3"
    outer_atom : indices
        atom indice in the mol which is to be rotated
        Ex., "7"
    True: Boolean
        True = Improper perturbation
        False = Valence perturbation
        Ex., True
    theta : float
        how many degrees by which to rotate

    [TODO]

    """
    #Set perturbation type
    if angle_type == True:
        angle_type = perturb_improper
    else:
        angle_type = perturb_valence


    #determine which improper angle on the atom is of interest and store the coordinates of the improper in  various variables
    print("161")
    cmol = oechem.OEMol(mol)
    print("before Conformer test")
    print("Conformer test")
    center_atom = cmol.GetAtom(oechem.OEHasAtomIdx(central_atom))
    move_atom = cmol.GetAtom(oechem.OEHasAtomIdx(outer_atom))
    print("165")
    print(center_atom)
    print(move_atom)
    center_coord = np.array(cmol.GetCoords(center_atom))
    move_coord = np.array(cmol.GetCoords(move_atom))
    print(center_coord)

    other_coords = list()
    for neighbor in center_atom.GetAtoms():
        if neighbor.GetName() != outer_atom:
            new_coord = np.array(cmol.GetCoords(neighbor))
            other_coords.append(new_coord)
    print(other_coords[0])

    atom0, atom1, atom2, atom3_rot, rot_mat = angle_type(center_coord, other_coords[0], other_coords[1], move_coord, theta)
    cmol.SetCoords(move_atom, oechem.OEFloatArray(atom3_rot))

    # Create an output file for visualization: .pdb file
    angle = str(theta)
    perturbation_type = str(angle_type)
    ofile = oechem.oemolostream(angle + '_' + perturbation_type + 'molecule.mol2')
    oechem.OEWriteMolecule(ofile, cmol)
    ofile.close()

    return cmol




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
    # print details of geometry
    if verbose:
        print("\n>>> Perturbing valence while maintaining improper angle...")
        # atom being moved is atom3.
        print("\natom0 original coords:\t",atom0)
        print("atom1 original coords:\t",atom1)
        print("atom2 original coords:\t",atom2)
        print("atom3 original coords:\t",atom3)
        print("atom3 {:.2f} deg rotated: {}".format(60.,atom3_rot))

        # check improper but make sure the central is first and moved is last
        print("\nImproper angle, before: ", calc_improper_angle(atom0, atom1, atom2, atom3))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom2, atom3, atom1))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom3, atom1, atom2))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom1, atom2, atom3_rot))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom2, atom3_rot, atom1))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom3_rot, atom1, atom2))

        ## check valences
        print("\nvalence angle 1, before: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, before: ", calc_valence_angle(atom0, atom1, atom3))
        print("valence angle 3, before: ", calc_valence_angle(atom0, atom2, atom3))
        print("valence angle 1, after: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, after: ", calc_valence_angle(atom0, atom1, atom3_rot))
        print("valence angle 3, after: ", calc_valence_angle(atom0, atom2, atom3_rot))
        print()

    return atom0, atom1, atom2, atom3_rot, rot_mat



# todo [8] write a test for your function
mol = oechem.OEMol()
oechem.OESmilesToMol(mol, 'FNCl')
oechem.OEAddExplicitHydrogens(mol)
omega = oeomega.OEOmega()
omega.SetMaxConfs(1)
omega(mol)
print('fining improper...')
a= find_improper_angles(mol)
print('a', a)


# test the perturbations at different angles
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  0)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  20)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  40)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  60)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  80)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  100)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  120)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  140)
oemol_perturb(mol, a[1][0][0], a[1][0][1], True,  160)


oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  0)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  20)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  40)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  60)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  80)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  100)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  120)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  140)
oemol_perturb(mol, a[1][0][0], a[1][0][1], False,  160)

