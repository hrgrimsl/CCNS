#Input: Molecular data, basis, etc.
#Output: Second-order Taylor Expansion About Mean-field Solution

import numpy as np
#opt_einsum.contract() is just np.linalg.einsum() but more efficient
from opt_einsum import contract
import copy
import psi4

class Molecule:
    def __init__(self, geometry, basis):
        #compute and store all 1- and 2- electron integrals; currently memory-limiting
        #needs to be implemented with on-the-fly calculations or with density fitting
        molecule = psi4.geometry(geometry)
        psi4.core.set_output_file('output.dat', False)
        psi4.set_options({'basis': str(basis), 'scf_type': 'pk', 'reference': 'rhf'})
        if molecule.multiplicity != 1:
            psi4.set_options({'reference': 'rohf'})
        self.hf_energy, wfn = psi4.energy('scf', return_wfn=True)
        mints = psi4.core.MintsHelper(wfn.basisset())
        Ca = wfn.Ca()
        Cb = wfn.Cb()
        self.multiplicity = molecule.multiplicity
        if molecule.multiplicity == 1:
            Cb = Ca
        #f_ = fock matrix
        #J_ _ _ _ = 2-electron matrix if electrons of opposite spin <ij|ab>
        #L_ _ _ _ = 2-electron matrix if electrons of same spin <ij||ab>
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
        Kaaaa = copy.copy(self.Jaaaa)
        Kabab = copy.copy(self.Jabab)
        Kbaba = copy.copy(self.Jbaba)
        Kbbbb = copy.copy(self.Jbbbb)
        Kaaaa = Kabab.swapaxes(2, 3)
        Kabab = Kabab.swapaxes(2, 3)
        Kbaba = Kbaba.swapaxes(2, 3)
        Kbbbb = Kbbbb.swapaxes(2, 3)
        self.Laaaa = self.Jaaaa - Kaaaa
        self.Labab = self.Jabab - Kabab
        self.Lbaba = self.Jbaba - Kbaba
        self.Lbbbb = self.Jbbbb - Kbbbb
        #Number_Occupied/Virtual_alpha/beta:
        self.NOa = wfn.Ca_subset("AO", "ACTIVE_OCC").shape[1]
        self.NOb = wfn.Cb_subset("AO", "ACTIVE_OCC").shape[1]
        self.NVa = wfn.Ca_subset("AO", "ACTIVE_VIR").shape[1]
        self.NVb = wfn.Cb_subset("AO", "ACTIVE_VIR").shape[1]
        #Initialize trial vector t and gradient vector g
        self.Build_Trial()
        self.Build_Gradient()
        assert(self.taa.shape == self.gaa.shape)
        assert(self.tbb.shape == self.gbb.shape)
        assert(self.taaaa.shape == self.gaaaa.shape)
        assert(self.tabab.shape == self.gabab.shape)
        assert(self.tbaba.shape == self.gbaba.shape)
        assert(self.tbbbb.shape == self.gbbbb.shape)

    def Build_Trial(self):
        #
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
        self.gaaaa = 2**.5*self.Laaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:]
        BBOO = np.triu_indices(self.NOb, k = 1)
        BBVV = np.triu_indices(self.NVb, k = 1)
        self.gbbbb = 2**.5*self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,self.NOb:]
        AO = np.arange(self.NOa)
        BO = np.arange(self.NOb)
        AV = np.arange(self.NOa,self.NVa+self.NOa)
        BV = np.arange(self.NOb,self.NVb+self.NOb)
        self.gabab = 2**.5*self.Jabab[AO,:,:,:][:,BO,:,:][:,:,AV,:][:,:,:,BV]
        self.gbaba = 2**.5*self.Jbaba[BO,:,:,:][:,AO,:,:][:,:,BV,:][:,:,:,AV]

    def Hessian_Action(self, vector):
        #This will consist of 6 sets of 6 contractions, and should be of identical shape to the gradient and trial vectors.
        haa = self.HAA(vector)
        hbb = self.HBB(vector)
        haaaa = self.HAAAA(vector)
        #habab = self.HABAB(vector)
        #hbaba = self.HBABA(vector)
        hbbbb = self.HBBBB(vector)
    def HAA(self, v):
        #Sums across every interaction with the taa terms
        h = self.AA_AA(v['aa'])
        h += self.AA_BB(v['bb'])
        h += self.AA_AAAA(v['aaaa'])
        h += self.AA_ABAB(v['abab'])
        h += self.AA_BABA(v['baba'])
        h += self.AA_BBBB(v['bbbb'])
        return h
        #h += self.BBBB_AA(v)
    def HBB(self, v):
        h = self.BB_AA(v['aa'])
        h += self.BB_BB(v['bb'])
        h += self.BB_AAAA(v['aaaa'])
        h += self.BB_ABAB(v['abab'])
        h += self.BB_BABA(v['baba'])
        h += self.BB_BBBB(v['bbbb'])
        return h
    def HAAAA(self, v):
        h = self.AAAA_AA(v['aa'])
        h = self.AAAA_BB(v['bb'])
        h = self.AAAA_AAAA(v['aaaa'])
        h = self.AAAA_ABAB(v['abab'])
        h = self.AAAA_BABA(v['baba'])
        h = self.AAAA_BBBB(v['bbbb'])
    def HBBBB(self, v):
        h = self.BBBB_AA(v['aa'])
        h += self.BBBB_BB(v['bb'])
        h += self.BBBB_AAAA(v['aaaa'])
        h += self.BBBB_ABAB(v['abab'])
        h += self.BBBB_BABA(v['baba'])
        h += self.BBBB_BBBB(v['bbbb'])
        return h
        #h += self.BBBB_AA(v)

    def AA_AA(self,v):
        h = np.zeros(self.taa.shape)

        #HAB term
        h += contract('ijab,ia->jb', self.Laaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:], v)

        h += contract('ba,ia->ib', self.fa[self.NOa:, self.NOa:], v)
        h -= contract('ij,ia->ja', self.fa[:self.NOa, :self.NOa], v)
        h += contract('ibaj,ia->jb', self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], v)
        return h

    def AA_BB(self,v):
        h = np.zeros(self.taa.shape)

        #HAB term
        h += contract('ijab,ia->jb', self.Jbaba[:self.NOb,:self.NOa,self.NOb:,self.NOa:], v)

        h += contract('ibaj,ia->jb', self.Jbaba[:self.NOb,self.NOa:,self.NOb:,:self.NOa], v)
        return h
    def AA_AAAA(self,v):
        h = np.zeros(self.taa.shape)
        h += contract('cjab,ijab->ic', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,self.NOa:], v)
        h -= contract('ijkb,ijab->ka', self.Laaaa[:self.NOa,:self.NOa,:self.NOa,self.NOa:], v)
        return h
    def AA_ABAB(self,v):
        h = np.zeros(self.taa.shape)
        h += contract('cjab,ijab->ic', self.Jabab[self.NOa:,:self.NOb,self.NOa:,self.NOb:], v)
        h -= contract('ijkb,ijab->ka', self.Jabab[:self.NOa,:self.NOb,:self.NOa,self.NOb:], v)
        return h
    def AA_BABA(self,v):
        h = np.zeros(self.taa.shape)

        return h
    def AA_BBBB(self,v):
        h = np.zeros(self.taa.shape)
        return h
    def BB_AA(self, v):
        h = np.zeros(self.taa.shape)

        # HAB term
        h += contract('ijab,ia->jb', self.Jabab[:self.NOa, :self.NOb, self.NOa:, self.NOb:], v)

        h += contract('ibaj,ia->jb', self.Jabab[:self.NOa, self.NOb:, self.NOa:, :self.NOb], v)
        return h
    def BB_AA(self,v):
        h = np.zeros(self.tbb.shape)

        #HAB term
        h += contract('ijab,ia->jb', self.Jabab[:self.NOa,:self.NOb,self.NOa:,self.NOb:], v)

        h += contract('ibaj,ia->jb', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], v)
        return h
    def BB_BB(self,v):
        h = np.zeros(self.tbb.shape)

        #HAB term
        h += contract('ijab,ia->jb', self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,self.NOb:], v)

        h += contract('ba,ia->ib', self.fb[self.NOb:, self.NOb:], v)
        h -= contract('ij,ia->ja', self.fb[:self.NOb, :self.NOb], v)
        h += contract('ibaj,ia->jb', self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], v)
        return h
    def BB_ABAB(self,v):
        h = np.zeros(self.tbb.shape)
        return h
    def BB_BABA(self,v):
        h = np.zeros(self.tbb.shape)
        h += contract('cjab,ijab->ic', self.Jbaba[self.NOb:,:self.NOa,self.NOb:,self.NOa:], v)
        h -= contract('ijkb,ijab->ka', self.Jbaba[:self.NOb,:self.NOa,:self.NOb,self.NOa:], v)
        return h
    def BB_AAAA(self,v):
        h = np.zeros(self.tbb.shape)
        return h
    def BB_BBBB(self,v):
        h = np.zeros(self.tbb.shape)
        h += contract('cjab,ijab->ic', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,self.NOb:], v)
        h -= contract('ijkb,ijab->ka', self.Lbbbb[:self.NOb,:self.NOb,:self.NOb,self.NOb:], v)
        return h

    def AAAA_AA(self,v):
        h = np.zeros(self.taaaa.shape)
        h -= contract('ibkj,ia->kjab', self.Laaaa[:self.NOa,self.NOa:,:self.NOa,:self.NOa],v)
        h += contract('iakj,ib->kjab', self.Laaaa[:self.NOa,self.NOa:,:self.NOa,:self.NOa],v)
        h += contract('cbaj,ia->ijcb', self.Laaaa[self.NOa:,self.NOa:,self.NOa:,:self.NOa],v)
        h -= contract('cbai,ja->ijcb', self.Laaaa[self.NOa:,self.NOa:,self.NOa:,:self.NOa],v)
        return h

    def AAAA_BB(self, v):
        h = np.zeros(self.taaaa.shape)
        return h

    def AAAA_AAAA(self,v):
        h = np.zeros(self.taaaa.shape)
        h += .5*contract('klab,ijab->klab', self.Laaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:],v)
        h += .5*contract('cdab,ijab->ijcd', self.Laaaa[self.NOa:,self.NOa:,self.NOa:,self.NOa:],v)
        h -= contract('cjal,ijab->ilcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa],v)
        h += contract('cjai,ljab->ilcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa],v)
        h += contract('bjal,ijac->ilcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa],v)
        h -= contract('bjai,ljac->ilcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa],v)
        h += contract('jcbk,ijab->ikac', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, :self.NOa], v)
        h -= contract('jabk,ijcb->ikac', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, :self.NOa], v)
        h -= contract('jcbi,kjab->ikac', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, :self.NOa], v)
        h += contract('jabi,kjcb->ikac', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, :self.NOa], v)
        h -= contract('jl,ijab->ilab', self.fa[:self.NOa,:self.NOa],v)
        h += contract('ji,ljab->ilab', self.fa[:self.NOa,:self.NOa],v)
        h += contract('cb,ijab->ijac', self.fa[self.NOa:,self.NOa:],v)
        h -= contract('ab,ijcb->ijac', self.fa[self.NOa:,self.NOa:],v)
        return h

    def AAAA_ABAB(self, v):
        h = np.zeros(self.taaaa.shape)
        return h

    def AAAA_BABA(self,v):
        h = np.zeros(self.taaaa.shape)
        return h

    def AAAA_BBBB(self, v):
        h = np.zeros(self.taaaa.shape)
        return h

    def BBBB_AA(self,v):
        pass
    def BBBB_BB(self,v):
        h = np.zeros(self.tbbbb.shape)
        h -= contract('ibkj,ia->kjab', self.Lbbbb[:self.NOb,self.NOb:,:self.NOb,:self.NOb],v)
        h += contract('iakj,ib->kjab', self.Lbbbb[:self.NOb,self.NOb:,:self.NOb,:self.NOb],v)
        h += contract('cbaj,ia->ijcb', self.Lbbbb[self.NOb:,self.NOb:,self.NOb:,:self.NOb],v)
        h -= contract('cbai,ja->ijcb', self.Lbbbb[self.NOb:,self.NOb:,self.NOb:,:self.NOb],v)
        return h
    def BBBB_AAAA(self,v):
        pass
    def BBBB_ABAB(self,v):
        pass
    def BBBB_BABA(self,v):
        pass
    def BBBB_BBBB(self,v):
        h = np.zeros(self.bbbb.shape)
        h += .5*contract('klab,ijab->klab', self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,self.NOb:],v)
        h += .5*contract('cdab,ijab->ijcd', self.Lbbbb[self.NOb:,self.NOb:,self.NOb:,self.NOb:],v)
        h -= contract('cjal,ijab->ilcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb],v)
        h += contract('cjai,ljab->ilcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb],v)
        h += contract('bjal,ijac->ilcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb],v)
        h -= contract('bjai,ljac->ilcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb],v)
        h += contract('jcbk,ijab->ikac', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, :self.NOb], v)
        h -= contract('jabk,ijcb->ikac', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, :self.NOb], v)
        h -= contract('jcbi,kjab->ikac', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, :self.NOb], v)
        h += contract('jabi,kjcb->ikac', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, :self.NOb], v)
        h -= contract('jl,ijab->ilab', self.fa[:self.NOb,:self.NOb],v)
        h += contract('ji,ljab->ilab', self.fa[:self.NOb,:self.NOb],v)
        h += contract('cb,ijab->ijac', self.fa[self.NOb:,self.NOb:],v)
        h -= contract('ab,ijcb->ijac', self.fa[self.NOb:,self.NOb:],v)
        return h

if __name__ == '__main__':
    geometry = """
        0 2
        H 0 0 0
        H 0 0 1
        H 0 0 2
    symmetry c1
    """
    basis = 'STO-3G'
    mol = Molecule(geometry, basis)
    trial = {'aa': mol.taa, 'bb': mol.tbb, 'aaaa': mol.taaaa, 'abab': mol.tabab, 'baba': mol.tbaba, 'bbbb': mol.tbbbb}
    mol.Hessian_Action(trial)
