### This script will take an .xyz (scan-final.xyz) from the finished QM torsion scan and generate/process data to create a plot of QM and MM energies.
#Imports
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye.oechem import *
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from calc_improper import find_improper_angles, calc_improper_angles, angle_betwen
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


def GetMM(oemol, FF, tag):
    """
    Description:
    Takes in an oemol and calculates the MM energies with two forcefields, one which has removed
    nitrogen improper parameters. The data from these calculations is stored in the oemol object
    as MM and MMRmNit. The units of energy are KCal/mol.

    Input:
    oemol: A single oemol object
    FF: .offxml file of smirnoff99Frosst.offxml
    tag: The name to tag the energy as, ex: "MM"

    Return:
    oemol: An oemol with the MM energies stored in tag in units kcal/mol
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
    oemol.SetData(tag, energy)

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
            crd1 = np.asarray(coords[atoms_imp[0]])
            crd2 = np.asarray(coords[atoms_imp[1]])
            crd3 = np.asarray(coords[atoms_imp[2]])
            crd4 = np.asarray(coords[atoms_imp[3]])
            angle_imp = calc_improper_angle(crd1, crd2, crd3, crd4, True)
            if angle_imp < 90 and angle_imp > 0:
                oemol.SetData("improper", angle_imp)
            if angle_imp < 0:
                if angle_imp < -90:
                    #angle_imp = 180 + angle_imp
                    oemol.SetData("improper", angle_imp)
            if angle_imp > 90:
                angle_imp = angle_imp - 180
                oemol.SetData("improper", angle_imp)

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

    #iterate through the list and subtract lowest mm and qm energies
    for m in oemolList:
        qm = m.GetData(qm_tag) - low_qm
        m.SetData(qm_tag, qm)
        mm = m.GetData(mm_tag) - low_mm
        m.SetData(mm_tag, mm)

    return oemolList




def plotResultsMulti(oemolList, names, tagList, colorList, markers):
    """
    Description: This function plots data for multiple molecules in a tag list.

    Input:
    oemolList: List of oemols with tags of energy data
    names: List of molecule names to title plot
    tagList: List of tagged data of interest
    colorList: List of colors for plot corresponding with the tags
    markers: Specified marker types for the plot
    """
    for mols, title in zip(oemolList, names):
        sort_oemol = sorted(mols, key=lambda x: x.GetData('improper'))
        xs = [m.GetData('improper') for m in sort_oemol]
        qm = [m.GetData('QM') for m in sort_oemol]
        plt.plot(xs, qm, color="royalblue", label='QM', marker= "^", markersize="8", linewidth="3")
        for tag, col, mar in zip(tagList, colorList, markers):
            trend = str(tag)
            name = [m.GetData(tag) for m in sort_oemol]
            plt.plot(xs, name, color=col, label=trend, marker=mar,markersize="8", linewidth="3")
        plt.title(title)
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



def compareGromacs(oemolList, FF, topfile, grofile):
    """
    Description:
    Compares energies of openMM to gromacs

    input:
    oemolList: List of oemol objects to compare
    FF: openmm force field, .offxml file
    topfile: GROMACS topology file name to write
    grofile: GROMACS coordinate file name (.gro format) to write


    """

    for idx, mol in enumerate(oemolList):
        name = str(idx) + "_" + topfile
        name_gro = str(idx) + "_" + grofile
        ff = ForceField(FF)
        topology = generateTopologyFromOEMol(mol)
        system = ff.createSystem(topology, [mol])
        positions = extractPositionsFromOEMol(mol)
        save_system_to_gromacs(topology, system, positions, name, name_gro)


        top = parmed.load_file(name)
        gromacssys = top.createSystem(nonbondedMethod= app.NoCutoff, constraints = None, implicitSolvent = None)
        gro = parmed.load_file(name_gro)

        #smirnoff
        smirfftop, smirffsys, smirffpos = create_system_from_molecule(ff, mol, verbose = False)

        print(oechem.OEMolToSmiles(mol))
        #compare
        try:
            groups0, groups1, energy0, energy1 = compare_system_energies( smirfftop, top.topology, smirffsys, gromacssys, positions, verbose = False)
            print(groups0, groups1, energy0, energy1)
        except:
            print(oechem.OEMolToSmiles(mol))
            print("Failed")





def sdf2mol2(sdf, output):
    """
    Takes in .sdf file and writes to a .mol2 file
    input:
    sdf: .sdf file
    output: name of output .mol2 file
    """

    ifs = oechem.oemolistream(sdf)
    ofs = oechem.oemolostream(output)

    ifs.SetFormat(oechem.OEFormat_SDF)
    ofs.SetFormat(oechem.OEFormat_MOL2)

    for mol in ifs.GetOEGraphMols():
            oechem.OEWriteMolecule(ofs, mol)



def oeb2oemol(oebfile):
    """
    Takes in oebfile and generates oemolList
    input:
    oebfile: Oebfile (.oeb)

    return:
    mollist: List of oemols

    """

    ifs = oechem.oemolistream(oebfile)
    mollist = []

    for mol in ifs.GetOEGraphMols():
            mollist.append(oechem.OEGraphMol(mol))

    return mollist


def processData(smiles, xyzfile, constraintFile, moltitle, FF, tags, color):
    """
    Description:
    Takes in innitial .sdf file, xyzfile from geomeTRIC output, constraint file and a title to
    create an .oeb file with oemols with QM, and MM data.


    input:
    #UPDATE to .sdf or .mol2 file
    smiles: Smiles string for molecule
    xyzfile: output .xyz file from geomeTRIC scan, scan-final.xyz that contains QM energies in title
    constraintFile: constraint.txt for geomeTRIC input that contains the indicies constrained in scan
    molTitle: title for the output .oeb file
    FF: Forcefield list
    tags: List of tags corresponding to force field
    color: List of colors used for plot

    Return:
    none, generates .oeb file and plot of data
    """

    #generate lsit of oemols with stored QM energies
    oemolList = QM2Oemol(xyzfile, SDF2oemol(smiles))

    #get MM energies from force fields
    for f, tag in zip(FF,tags):
        for mol in oemolList:
            print("THIS IS THE FORCE FIELD" + str(f))
            GetMM(mol, f, tag)
            findAngles(mol, constraintFile)

    #normalize eneriges:
    for tag in tags:
        adjust_energy(oemolList, 'QM', tag)

    #write out oeb file
    makeOEB(oemolList, moltitle)

    return oemolList


def makeInputs(molName):
    '''
    This function makes the input files in a format to process MM and QM energies of groups of molecules

    Input:
    molName: List of molecules name in a set of molecules

    Return
    molName: Original molName list
    SDFList: List of sdf files for molecule in list
    contList: List of constraint files for molecule
    XYZList: List of .xyz scan files
    '''
    SDFList = []
    contList = []
    XYZList = []

    for name in molName:
        sdfName = str(name) + '.sdf'
        SDFList.append(sdfName)

        contName = 'c-' + str(name) + '.txt'
        contList.append(contName)

        XYZName = 'scan-final-' + str(name) + '.xyz'
        XYZList.append(XYZName)

    return molName, SDFList, contList, XYZList


def processSet(filesList, molNames, FFList, tagNames, colorList, markerList):
    '''
    This function takes in information for a set of molecules and generates
    subplots for all of the molecules from MM energies from the specified FFs.

    input:
    filesList: A list of molNames, SDF file of starting molecule, constraints file list, xyz files with QM scans
    molNames: List of oemol names
    FFList : List of FFs of interest
    tagNames: name of tags associated with FFs
    colorList: List of colors for the plot
    markerList: Type of markers for pots
    '''

    names = filesList[0]
    sdf = filesList[1]
    cont = filesList[2]
    xyz = filesList[3]

    mols = []

    i = 0
    while i < len(filesList[0]):
        doneList = processData(sdf[i], xyz[i], cont[i], names[i], FFList, tagNames, colorList)
        mols.append(doneList)
        i += 1

    #plot the list of list of oemols
    plotResultsMulti(mols, molNames, tagNames, colorList, markerList)




