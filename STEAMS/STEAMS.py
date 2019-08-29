#Input: Molecular data, basis, etc.
#Output: Second-order Taylor Expansion About Mean-field Solution

import numpy as np
from opt_einsum import contract
import copy
import psi4

class Molecule:
    def __init__(self, geometry, basis):
        #Get Fock matrices fa/fb, as well as exchange matrices, e.g. Labab
        molecule = psi4.geometry(geometry)
        psi4.core.set_output_file('output.dat', False)
        psi4.set_options({'basis': str(basis), 'scf_type': 'pk', 'reference': 'rhf'})
        self.hf_energy, wfn = psi4.energy('scf', return_wfn=True)
        mints = psi4.core.MintsHelper(wfn.basisset())
        Ca = wfn.Ca()
        Cb = wfn.Cb()
        self.multiplicity = molecule.multiplicity
        if molecule.multiplicity == 1:
            Cb = Ca
        self.fa = wfn.Fa()
        self.fb = wfn.Fb()
        self.Jaaaa = np.array(mints.mo_eri(Ca, Ca, Ca, Ca))
        self.Jabab = np.array(mints.mo_eri(Ca, Ca, Cb, Cb))
        self.Jbaba = np.array(mints.mo_eri(Cb, Cb, Ca, Ca))
        self.Jbbbb = np.array(mints.mo_eri(Cb, Cb, Cb, Cb))
        self.fa.transform(Ca)
        self.fa = np.array(self.fa)
        self.fb = np.array(self.fb)
        self.Jaaaa = self.Jaaaa.swapaxes(1, 2)
        self.Jabab = self.Jabab.swapaxes(1, 2)
        self.Jbaba = self.Jbaba.swapaxes(1, 2)
        self.Jbbbb = self.Jbbbb.swapaxes(1, 2)
        self.Kaaaa = copy.copy(self.Jaaaa)
        self.Kabab = copy.copy(self.Jabab)
        self.Kbaba = copy.copy(self.Jbaba)
        self.Kbbbb = copy.copy(self.Jbbbb)
        self.Kaaaa = self.Kaaaa.swapaxes(2, 3)
        self.Kabab = self.Kabab.swapaxes(2, 3)
        self.Kbaba = self.Kbaba.swapaxes(2, 3)
        self.Kbbbb = self.Kbbbb.swapaxes(2, 3)
        self.Laaaa = self.Jaaaa - self.Kaaaa
        self.Labab = self.Jabab - self.Kabab
        self.Lbaba = self.Jbaba - self.Kbaba
        self.Lbbbb = self.Jbbbb - self.Kbbbb
        self.NOa = wfn.Ca_subset("AO", "ACTIVE_OCC").shape[1]
        self.NOb = wfn.Cb_subset("AO", "ACTIVE_OCC").shape[1]
        self.NVa = wfn.Ca_subset("AO", "ACTIVE_VIR").shape[1]
        self.NVb = wfn.Cb_subset("AO", "ACTIVE_VIR").shape[1]
        #Initialize trial vector t and gradient vector g
        self.Build_Trial()
        self.Build_Gradient()

    def Build_Trial(self):
        self.taa = np.zeros((self.NOa, self.NVa))
        self.tbb = np.zeros((self.NOb, self.NVb))
        self.taaaa = np.zeros((self.NOa, self.NOa, self.NVa, self.NVa))
        self.tabab = np.zeros((self.NOa, self.NOb, self.NVa, self.NVb))
        self.tbaba = np.zeros((self.NOb, self.NOa, self.NVb, self.NVa))
        self.tbbbb = np.zeros((self.NOb, self.NOb, self.NVb, self.NVb))

    def Build_Gradient(self):
        self.gaa = np.zeros((self.NOa, self.NVa))
        self.gbb = np.zeros((self.NOb, self.NVb))
        AAOO = np.triu_indices(self.NOa, k = 1)
        AAVV = np.triu_indices(self.NVa, k = 1)
        self.gaaaa = 2**.5*self.Laaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:][AAOO][(slice(None),) + AAVV].flatten()
        BBOO = np.triu_indices(self.NOb, k = 1)
        BBVV = np.triu_indices(self.NVb, k = 1)
        self.gbbbb = 2**.5*self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,self.NOb:][BBOO][(slice(None),) + BBVV].flatten()
        AO = np.arange(self.NOa)
        BO = np.arange(self.NOb)
        AV = np.arange(self.NOa,self.NVa+self.NOa)
        BV = np.arange(self.NOb,self.NVb+self.NOb)
        self.gabab = 2**.5*self.Jabab[AO,:,:,:][:,BO,:,:][:,:,AV,:][:,:,:,BV].flatten()
        self.gbaba = 2**.5*self.Jbaba[BO,:,:,:][:,AO,:,:][:,:,BV,:][:,:,:,AV].flatten()

if __name__ == '__main__':
    geometry = """
        H 0 0 0
        H 0 0 1
        H 0 0 2
        H 0 0 3
    symmetry c1
    """
    basis = 'STO-3G'
    mol = Molecule(geometry, basis)