from off_nitrogens.calc_improper import *
import numpy as np

def test_two_vectors():
    """Test the angle between two simple vectors."""
    ang = angle_between((1, 0, 0), (0, 1, 0))
    if ang != 90.0:
        raise Exception("Wrong angle calculated between two vectors.")

def test_planar():
    """Planar test coordinates should return a zero improper measurement."""
    atom1 = np.asarray([0,0,0])
    atom0 = np.asarray([0,1,0])
    atom2 = np.asarray([-0.866,-0.5,0])
    atom3 = np.asarray([0.866,-0.5,0])
    ang = calc_improper_angle(atom1, atom0, atom2, atom3)
    if ang != 0.0:
        raise Exception("Planar set of coordinates is not returning a 0.0 improper dihedral angle.")

def test_above_plane():
    """Test coordinates in which central atom is above the plane with outer atoms."""
    atom1 = np.asarray([0,0,0])
    atom0 = np.asarray([0,1,-1])
    atom2 = np.asarray([-0.866,-0.5,-1])
    atom3 = np.asarray([0.866,-0.5,-1])
    ang = calc_improper_angle(atom1, atom0, atom2, atom3)
    if abs(ang-63.43) > 0.01:
        raise Exception("Improper dihedral with central atom above plane is incorrect.")

def test_below_plane():
    """Test coordinates in which central atom is below the plane with outer atoms."""
    atom1 = np.asarray([0,0,0])
    atom0 = np.asarray([0,1,1])
    atom2 = np.asarray([-0.866,-0.5,1])
    atom3 = np.asarray([0.866,-0.5,1])
    ang = calc_improper_angle(atom1, atom0, atom2, atom3)
    if abs(ang-116.56) > 0.01:
        raise Exception("Improper dihedral with central atom below plane is incorrect.")

def test_atom_order():
    """Test coordinates in which central atom is below the plane with outer atoms.
       Make sure that the results are the same with different order of outer atoms,
       as long as they maintain the same handedness, e.g., clockwise or counterclockwise."""
    atom1 = np.asarray([0,0,0])
    atom0 = np.asarray([0,1,1])
    atom2 = np.asarray([-0.866,-0.5,1])
    atom3 = np.asarray([0.866,-0.5,1])
    ang1 = calc_improper_angle(atom1, atom2, atom3, atom0)
    ang2 = calc_improper_angle(atom1, atom3, atom0, atom2)
    if abs(ang1-116.56) > 0.01 or abs(ang2-116.56) > 0.01:
        raise Exception("Changing order of outer atom specification (while maintaining handedness) gives wrong result.")

def test_atom_order():
    """Test coordinates for test OEMol NHFCl."""

    import openeye.oechem as oechem
    import openeye.oeomega as oeomega
    # coordinates for N, F, H, Cl respectively
    coordlist = [ 0.155,  -0.088,  -0.496,
                 -1.054,  -0.776,  -0.340,
                 -0.025,   0.906,  -0.516,
                  1.689,  -0.635,  -1.263]
    # create OEMol
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, 'FNCl')
    omega = oeomega.OEOmega()
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(False)
    status = omega(mol)
    oechem.OETriposAtomTypes(mol)
    oechem.OETriposAtomNames(mol)
    oechem.OEAddExplicitHydrogens(mol)
    omega.SetMaxConfs(1)
    # set provided coordinates
    mol.SetCoords(oechem.OEFloatArray(coordlist))
    # calculate and check improper angle
    atom1, atom0, atom2, atom3 = find_improper_angles(mol)[0]
    ang = calc_improper_angle(atom1, atom0, atom2, atom3)
    if abs(ang-15.0) > 0.1:
        raise Exception("Error calculating improper angle in test OEMol system.")
