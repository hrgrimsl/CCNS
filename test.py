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
    assert abs(mol.conj_grad()+460.1711971769602201)<1e-8
