###WIP:
### This script will take an .xyz (scan-final.xyz) from the finished QM torsion scan and generate/process data to create a plot of QM and MM energies.


#Imports
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye.oechem import *
#import oenotebook as oenb
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from calc_improper import *
import numpy as np
import matplotlib.pyplot as plt


#Functions
def SDF2oemol(sdffile):
    """
    Description:
    Takes in the sdf file and converts it to an oemol

    Input:
    sdffile: .sdf file

    Returns:
    Oemol: Oemol object of the molecule stored in the input .sdf file
    """
    #open file
    ifs = oechem.oemolistream()
    oemol = oechem.OECreateOEGraphMol()
    ifs.open(sdffile)
    #empty oemol list
    molecules = list()
    oechem.OEReadMolecule(ifs, oemol)
    ifs.close()
    print(oemol.GetCoords())
    return oemol

def smiles2oemol(smiles):
    """
    Description:
    Takes in an SMILES string and returns an oemol.

    Input:
    smiles: SMILES string for molecule

    Returns:
    mol1: OEMol object
    """
    mol1 = OEMol()
    OESmilesToMol(mol1, smiles)

    #assign charges, necessary to keep?
    chargeEngine = OEAM1BCCCharges()
    OEAssignCharges(mol1, chargeEngine)

    return mol1


def QM2Oemol(xyzfile, oemol):
    """
    Description:
    Takes in an xyz file and creates oemol objects with the coordinates in the original .xyz file
    with QM energies. The .xyz files are formatted as an output geomeTRIC xyz file, scan-final.xyz.    The original energies are stored in Hartree but are converted to kcal/mol in the QM energy tag.
    Input:
    xyzfile:Take in scan-final.xyz file which contains the geometry and energy outputs from a
    geomeTRIC torsion scan
    oemol:The oemol of the molecule involved in the torsion scan of the .xyz file

    Returns:
    oemolList: A list of OEmols that contain the QM data tagged as "QMEng"
    """

    #open input xyz file
    ifs = oechem.oemolistream()
    ifs.open(xyzfile)

    #generate empty list of oemols
    oemolList=list()

    #iterate through oemols and set coordinates to original oemol
    for mol in ifs.GetOEGraphMols():
        #get coordinates
        coords = mol.GetCoords()
        #set coordinates of original oemol with correct bond connectivity
        newmol = copy.deepcopy(oemol)
        newmol.SetCoords(coords)
        #get the energies from the .xyz file title and save in oemol data
        title = mol.GetTitle()
        energy = float(str.split(title)[-1])
        energy_kcal= energy * 627.509
        #set the QM energy in a tag called qm
        newmol.SetData("QM", energy_kcal)

        #append the list with the new mol that has stored QM energies
        oemolList.append(newmol)

    return oemolList


def GetMM(oemol, FF, FFRmNit):
    """
    Description:
    Takes in an oemol and calculates the MM energies with two forcefields, one which has removed
    nitrogen improper parameters. The data from these calculations is stored in the oemol object
    as MM and MMRmNit. The units of energy are KCal/mol.

    Input:
    oemol: A single oemol object
    FF: .offxml file of smirnoff99Frosst.offxml
    FFRmNit: .offxml file of smirnoff99Frosst.offxml with the removed nitrogen improper parameter

    Return:
    oemol: An oemol with the MM energies stored in tags MM and MMRmNit in kcal/mol
    """

    #prep both force fields, create system and topology, get MM energies and store in data tag
    ff = ForceField(FF)
    topology = generateTopologyFromOEMol(oemol)
    system = ff.createSystem(topology, [oemol])
    positions = extractPositionsFromOEMol(oemol)
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy = True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("this is the original mm energy: " + str(energy))
    oemol.SetData("MM", energy)

    #repeat with removed nitrogen
    ffn = ForceField(FFRmNit)
    topology = generateTopologyFromOEMol(oemol)
    system_n = ffn.createSystem(topology, [oemol])
    positions = extractPositionsFromOEMol(oemol)
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    simulation = app.Simulation(topology, system_n, integrator)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy = True, getPositions=True)
    energy_rm = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("this is the original mm energy: " + str(energy))
    oemol.SetData("MMRmNit", energy_rm)

    return oemol

def findAngles(oemol, constraint):
    """
    Description:
    Calculates the improper (and potentially valence) angles of an oemol and stores the value in
    tag "improper" and "valence" for an oemol.

    Input:
    oemol: An oemol object
    constraint: An constraints.txt file used in an input for an geomeTRIC scan

    Return:
    oemol: An oemol with the calculated improper or valence angles with corresponding tags
    "improper" or "valence"
    """

    #determine the constraints
    #open the constraint file
    constraintfile = open(constraint, "r")
    f = open(constraint)
    lines = f.readlines()
    f.close()
    coords = oemol.GetCoords()

    for l in lines:
        if "dihedral" in l:
            split = l.split()
            atoms_imp = (int(split[1])-1, int(split[2])-1, int(split[3])-1, int(split[4])-1)
            print(atoms_imp)
            crd1 = np.asarray(coords[atoms_imp[0]])
            crd2 = np.asarray(coords[atoms_imp[1]])
            crd3 = np.asarray(coords[atoms_imp[2]])
            crd4 = np.asarray(coords[atoms_imp[3]])
            angle_imp = calc_improper_angle(crd1, crd2, crd3, crd4, True)
            print(angle_imp)
            if angle_imp < 90 and angle_imp > 0:
                oemol.SetData("improper", angle_imp)
            if angle_imp < 0:
                if angle_imp < -90:
                    #angle_imp = 180 + angle_imp
                    oemol.SetData("improper", angle_imp)
            if angle_imp > 90:
                angle_imp = angle_imp - 180
                oemol.SetData("improper", angle_imp)
            print(angle_imp)

        #calculate the valence angle and store in tag (2-d scan)
        if "angle" in l:
            split = l.split()
            atoms_val = (int(split[1])-1, int(split[2])-1, int(split[3])-1)
            crd1 = np.asarray(coords[atoms_val[0]])
            crd2 = np.asarray(coords[atoms_val[1]])
            crd3 = np.asarray(coords[atoms_val[2]])
            angle_val= calc_valence_angle(crd1, crd2, crd3)
            oemol.SetData("valence", angle_val)

    return oemol



#TODO
def makeOEB(oemolList, tag):
    """
    Description:
    Takes in an oemol list and creates an output OEB file.

    Input:
    oemolList: A list of oemols
    tag: The title of the OEB file.

    Return:

    """
    ofile = oemolostream(tag+'.oeb')
    for mol in oemolList:
        OEWriteConstMolecule(ofile, mol)
    ofile.close()
    return





def adjust_energy(oemolList, qm_tag, mm_tag):
    """
    Description:
    Iterates through oemols in a list and normalizes the energies to the lowest
    QM energy in the list. The corresponding MM energy geometry is subtracted from the
    MM energies tagged in the oemol.

    Input:
    oemolList: A list of oemols
    qm_tag: The tag for the qm energy
    mm_tag: The tag for the mm energy being normalized

    Return:
    oemolList with the updated energies.
    """

    #sort the oemols based on the lowest QM energy, this oemol is now "low_mol"
    low_mol = sorted(oemolList, key=lambda x: x.GetData(qm_tag))[0]
    #get corresponding lowest qm energy
    low_qm = low_mol.GetData(qm_tag)

    #get correspoding lowest mm energy
    #low_mol_mm = sorted(oemolList, key=lambda x: x.GetData(mm_tag))[0]
    low_mm = low_mol.GetData(mm_tag)
    print("this is the low mm value:" + str(low_mm))

    #iterate through the list and subtract lowest mm and qm energies
    for m in oemolList:
        qm = m.GetData(qm_tag) - low_qm
        m.SetData(qm_tag, qm)
        mm = m.GetData(mm_tag) - low_mm
        m.SetData(mm_tag, mm)

    return oemolList


def plotResults(oemolList):
    """
    Description: This function plots the QM and MM data for 1d torsion scans

    Input:
    oemolList: List of eomols with stored angles, MM, and QM energies
    """
    #sort oemolList based on the improper angles
    sort_oemol = sorted(oemolList, key=lambda x: x.GetData('improper'))

    xs = [m.GetData('improper') for m in sort_oemol]
    qm = [m.GetData('QM') for m in sort_oemol]
    mm = [m.GetData('MM') for m in sort_oemol]
    mm_nit = [m.GetData('MMRmNit') for m in sort_oemol]

    plt.plot(xs, qm, color="royalblue", label='QM', marker='^',markersize="8", linewidth="3")
    plt.plot(xs, mm, color="red", label='MM', marker='^',markersize="8", linewidth="3")
    plt.plot(xs, mm_nit, color="orange", label='MM removed Imp.', marker='^',markersize="8", linewidth="3")
    #plt.plot(xs, qm, color="royalblue", label='Carbon 1', marker='^',markersize="12", linewidth="5")
    plt.legend()
    plt.show()


def oeb2mollist(oeb):
    """
    Description:
    Takes in oeb file and creates oemolList

    Input:
    oeb: oeb file

    Return:
    oemolList: a list of eomols contained in the .oeb file
    """

    #open input xyz file
    ifs = oechem.oemolistream()
    ifs.open(oeb)

    #generate empty list of oemols
    oemolList=list()

    #iterate through oemols and set coordinates to original oemol
    for mol in ifs.GetOEGraphMols():
        oemolList.append(oechem.OEGraphMol(mol))

    return oemolList



def compareGromacs(oemolList, FF, top, gro):
    """
    Description:
    Compares energies of openMM to gromacs

    input:
    oemolList: List of oemol objects to compare
    FF: openmm force field, .offxml file
    top: GROMACS topology file name to write
    gro: GROMACS coordinate file name (.gro format) to write


    """
    for mol in oemolList:
        ff = ForceField(FF)
        topology = generateTopologyFromOEMol(mol)
        system = ff.createSystem(topology, [mol])
        positions = extractPositionsFromOEMol(mol)
        save_system_to_gromacs( topology, system, positions, top, gro )




#Main
def processData(smiles, xyzfile, constraintFile, moltitle, FF, FFRmNit):
    """
    Description:
    Takes in innitial .sdf file, xyzfile from geomeTRIC output, constraint file and a title to
    create an .oeb file with oemols with QM, and MM data.
    Also generates plot comparing QM and MM energies.


    input:
    #UPDATE to .sdf or .mol2 file
    smiles: Smiles string for molecule
    xyzfile: output .xyz file from geomeTRIC scan, scan-final.xyz that contains QM energies in title
    constraintFile: constraint.txt for geomeTRIC input that contains the indicies constrained in scan
    molTitle: title for the output .oeb file
    FF: Forcefield
    FFRmNit: Force field with removed nitrogens

    Return:
    none, generates .oeb file and plot of data
    """

    #generate lsit of oemols with stored QM energies
    oemolList = QM2Oemol(xyzfile, SDF2oemol(smiles))

    #get MM energies and improper angles from geometries (and potentially valence depending no the scan)
    for mol in oemolList:
        GetMM(mol, FF, FFRmNit)
        findAngles(mol, constraintFile)

    #normalize eneriges:
    adjust_energy(oemolList, 'QM', 'MM')
    adjust_energy(oemolList, 'QM', 'MMRmNit')

    #write out oeb file
    makeOEB(oemolList, moltitle)
    #plot results
    plotResults(oemolList)

    return







#plotResults(oeb2mollist('results.oeb'))

compareGromacs(QM2Oemol('scan-final.xyz', SDF2oemol('mol6.sdf')), 'smirnoff99Frosst.offxml', 'topology', 'coordfile')

#test
#processData('mol6.sdf', 'scan-final.xyz', 'constraints.txt', 'fixresults', 'smirnoff99Frosst.offxml', 'nitRem.offxml')



