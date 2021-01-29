import sys
sys.path.insert(0, '/export/home/jmaat/psi4_clone/01_scripts/')

import confs2psi
confs2psi.confs2psi('/export/home/jmaat/off_nitrogens/QM_inputs/pyrnit.sdf','mp2','def2-SV(P)',False,"1.5 Gb")
