Coupled Cluster Newton Step

<img src = "https://travis-ci.com/hrgrimsl/STEAMS.svg?token=y5H9g77PxszWJHZmEWzC&branch=master">

Module to use as a python API- e.g.
```python
import STEAMS
geometry = """
    symmetry c1
    0 1
    H 0 0 0
    Cl 0 0 0
"""
basis = 'aug-cc-pvdz'
mol = STEAMS.molecule(geometry, basis, rhf = True, uns = False)
print(mol.conj_grad())
```
Program uses a Psi4 backend to compute electron integrals, etc.  Note that if the UNS flag is set to False, one obtains the CEPA(0) energy with single and double excitations included.  See 

https://doi.org/10.1063/1.454125

and 

https://doi.org/10.1063/1.455824

for discussion on this topic.

See also the documentation at https://htmlpreview.github.io/?github.com/hrgrimsl/CCNS/blob/master/docs/build/index.html
