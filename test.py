from STEAMS import *
def test_1():
    geometry = """
        0 1
        H 0 0 0
        Cl 0 0 1
        symmetry c1
    """
    basis = 'cc-pvdz'
    mol = molecule(geometry, basis, rhf = True, uns = False)
    assert abs(mol.conj_grad()+460.17119717)<1e-7


