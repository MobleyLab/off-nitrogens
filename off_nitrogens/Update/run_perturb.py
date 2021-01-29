#!/usr/bin/env python
# I m p o r t s
from perturb_angle import *



# run the ~script!!~
mol2_files = {3:"molClass_pyrnit_molecule_3.xyz", 4:"molClass_pyrnit_molecule_4.xyz", 5:"molClass_pyrnit_molecule_5.xyz", 6:"molClass_pyrnit_molecule_6.xyz",  7:"molClass_pyrnit_molecule_7.xyz",  8:"molClass_pyrnit_molecule_8.xyz",  9:"molClass_pyrnit_molecule_9.xyz"}

input = input2mol(mol2_files)
print(input)
oemol_perturb(input, True, "pyrnit", 20, 1)


