from off_nitrogens.calc_improper import *
#from calc_improper import *

import numpy as np

def test_two_vectors():
    """Test the angle between two simple vectors."""
    ang = angle_between((1, 0, 0), (0, 1, 0))
    if ang != 90.0:
        raise Exception("Wrong angle calculated between two vectors.")

def test_valence():
    """Test the valence angle calculation of NHFCl molecule."""
    atom0 = np.asarray([ 0.155, -0.088, -0.496])
    atom1 = np.asarray([-1.054, -0.776, -0.340])
    atom2 = np.asarray([-0.025,  0.906, -0.516])
    atom3 = np.asarray([ 1.689, -0.635, -1.263])
    ang1 = calc_valence_angle(atom0, atom1, atom2)
    ang2 = calc_valence_angle(atom0, atom2, atom3)
    ang3 = calc_valence_angle(atom0, atom3, atom1)
    if abs(ang1-109.38) > 0.01:
        raise Exception("Incorrect valence angle calculated for NHFCl: atoms 0 1 2")
    if abs(ang2-116.25) > 0.01:
        raise Exception("Incorrect valence angle calculated for NHFCl: atoms 0 2 3")
    if abs(ang3-129.36) > 0.01:
        raise Exception("Incorrect valence angle calculated for NHFCl: atoms 0 3 1")

def test_planar():
    """Planar test coordinates should return a zero improper measurement."""
    atom0 = np.asarray([0,0,0])
    atom1 = np.asarray([0,1,0])
    atom2 = np.asarray([-0.866,-0.5,0])
    atom3 = np.asarray([0.866,-0.5,0])
    ang = calc_improper_angle(atom0, atom1, atom2, atom3)
    if ang != 0.0:
        raise Exception("Planar set of coordinates is not returning a 0.0 improper dihedral angle.")

def test_above_plane():
    """Test coordinates in which central atom is above the plane with outer atoms."""
    atom0 = np.asarray([0,0,0])
    atom1 = np.asarray([0,1,-1])
    atom2 = np.asarray([-0.866,-0.5,-1])
    atom3 = np.asarray([0.866,-0.5,-1])
    ang = calc_improper_angle(atom0, atom1, atom2, atom3)
    if abs(ang-63.43) > 0.01:
        raise Exception("Improper dihedral with central atom above plane is incorrect.")

def test_below_plane():
    """Test coordinates in which central atom is below the plane with outer atoms."""
    atom0 = np.asarray([0,0,0])
    atom1 = np.asarray([0,1,1])
    atom2 = np.asarray([-0.866,-0.5,1])
    atom3 = np.asarray([0.866,-0.5,1])
    ang = calc_improper_angle(atom0, atom1, atom2, atom3)
    if abs(ang-116.56) > 0.01:
        raise Exception("Improper dihedral with central atom below plane is incorrect.")

def test_atom_order():
    """Test coordinates in which central atom is below the plane with outer atoms.
       Make sure that the results are the same with different order of outer atoms,
       as long as they maintain the same handedness, e.g., clockwise or counterclockwise."""
    atom0 = np.asarray([0,0,0])
    atom1 = np.asarray([0,1,1])
    atom2 = np.asarray([-0.866,-0.5,1])
    atom3 = np.asarray([0.866,-0.5,1])
    ang1 = calc_improper_angle(atom0, atom2, atom3, atom1)
    ang2 = calc_improper_angle(atom0, atom3, atom1, atom2)
    if abs(ang1-116.56) > 0.01 or abs(ang2-116.56) > 0.01:
        raise Exception("Changing order of outer atom specification (while maintaining handedness) gives wrong result.")

def test_oemol_nhfcl():
    """Test coordinates for NHFCl read in as OEMol."""

    import openeye.oechem as oechem
    import openeye.oeomega as oeomega
    # coordinates for N, F, H, Cl respectively
    # generated after minimization with improper phase of 150 degrees
    coordlist = [ 0.155,  -0.088,  -0.496,
                 -1.054,  -0.776,  -0.340,
                 -0.025,   0.906,  -0.516,
                  1.689,  -0.635,  -1.263]
    # create OEMol
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, 'FNCl')
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(False)
    status = omega(mol)
    oechem.OETriposAtomTypes(mol)
    oechem.OETriposAtomNames(mol)
    oechem.OEAddExplicitHydrogens(mol)
    # set provided coordinates
    mol.SetCoords(oechem.OEFloatArray(coordlist))
    # calculate and check improper angle
    crds, names = find_improper_angles(mol)
    ang = calc_improper_angle(crds[0][0], crds[0][1], crds[0][2], crds[0][3])
    if abs(ang-15.0) > 0.1 and abs(ang-165.0) > 0.1:
        raise Exception("Error calculating improper of test OEMol. Calculated {} degrees, but should be 15 or 165 degrees.".format(ang))

