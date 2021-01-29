#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
QM_min.py


[Script 1 in Quipp Pipeline]
This script takes in a SMILES string and submits a QM minimization calculation for the initial.
"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import math
import sys
import shutil
from openeye import oeomega
from openeye import oechem

#Improper Quanformer, needs update for generalizabaility
import sys
sys.path.insert(0, '/export/home/jmaat/psi4_clone/01_scripts/indexedit/off_psi4/quanformer')
import confs2psi



#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================


#main
def QM_min(smileslist, molName, method, basisset):
    """
    Description: Takes in a dictionary of molecules in a smiles format and returns the minimized geometry coordinates in an .sdf file.

    Parameters:
    smileslist (dictionary): list of smiles strings
    molName (string): some classifer for the group of molecules
    method: The method for the QM calculation
    basissset: The basis set for the desired QM calculation
    """
    #create molName.sdf file
    smiles2sdf(smileslist,molName)

    #generate the
    confs2psi.confs2psi(molName + ".sdf", method, basisset, False, "1.5 Gb")

    #submit the jobs within the directory with the slurm script
    #copy one_psi.slurm into directory
    name =
    os.popen('cp one_psi.slurm /')


    return




def smiles2sdf(smiles, molClass):
    """
    Description:
    Takes in a list of smiles strings and creates an .sdf file of the omega generated conformers.
    This functions purpose is for input qm calculations to minimize.

    Parameters:
    smiles (dictionary): list of smiles strings
    molClass (string): some classifer for the group of molecules

    Returns:
    .sdf file of all the molecules from the smiles dictionary.
    The title of the molecules in the .sdf file is their key from the smiles dictionary.
    """

    #create .sdf output file
    ofile = oechem.oemolostream(molClass + '.sdf')

    #write smiles strings to molecules in the .sdf file
    for key, molecule in smiles.items():
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, molecule)
        oechem.OEAddExplicitHydrogens(mol)
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(1)
        omega(mol)
        mol.SetTitle("molClass_" + molClass + "_molecule_" + str(key))
        oechem.OEWriteConstMolecule(ofile, mol)
    ofile.close()
    return


#testing
test = {1: "CNC"}
QM_min(test, 'pyr', 'mp2', 'def2-SV(P)')






