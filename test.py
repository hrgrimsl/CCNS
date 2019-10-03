from STEAMS import *

def test_1():
    geometry = """
        0 1
        H 0 0 0
        Cl 0 0 1
        symmetry c1
    """
    basis = 'cc-pvdz'
    mol = Molecule(geometry, basis, RHF = True, UNS = False)
    assert abs(mol.conj_grad()+460.17119717)<1e-7

def test_2():
    geometry = """
        0 2
        Cl 0 0 1
        symmetry c1
    """
    basis = 'sto-3g'
    mol = Molecule(geometry, basis, RHF = False, UNS = True)
    assert abs(mol.conj_grad()+2.16972360)<1e-7
