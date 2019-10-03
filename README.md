#Second-Order Taylor Expansion About Mean-field Solution

<img src = "https://travis-ci.com/hrgrimsl/STEAMS.svg?token=y5H9g77PxszWJHZmEWzC&branch=master">

Module to use as a python API- e.g.

import STEAMS
geometry = """
    symmetry c1
    0 1
    H 0 0 0
    Cl 0 0 0
"""
basis = 'aug-cc-pvdz'
mol = STEAMS.Molecule(geometry, basis, RHF = True, UNS = False
print(mol.conj_grad())

