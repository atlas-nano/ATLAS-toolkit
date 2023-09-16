import sys
import os
import numpy as np
#from os.path import stripext
from ase.io import read, write

def stripext(p):
    return os.path.splitext(p)[0]

if len(sys.argv) < 5:
    print "Usage: xyz2vasp.py Lx Ly Lz file1.xyz [file2.xyz [...]]"

try:
    periodicLengths = [float(x) for x in sys.argv[1:4]]
    pv = np.diag(periodicLengths)
except ValueError:
    print "Usage: xyz2vasp.py Lx Ly Lz file1.xyz [file2.xyz [...]]"

for input in sys.argv[4:]:
    output = stripext(input)[0] + '.vasp'
    mol = read(input, format='xyz')
    mol.set_pbc(True)
    mol.set_cell(pv)
    write(output, mol, format='vasp')
