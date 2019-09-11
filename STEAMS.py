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
        self.AA = self.taa
        self.BB = self.tbb
        self.AAAA = self.taaaa
        self.ABAB = self.tabab
        self.BABA = self.tbaba
        self.BBBB = self.tbbbb

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
        self.AA = np.zeros(self.taa.shape)
        self.BB = np.zeros(self.tbb.shape)
        self.AAAA = np.zeros(self.taaaa.shape)
        self.ABAB = np.zeros(self.tabab.shape)
        self.BABA = np.zeros(self.tbaba.shape)
        self.BBBB = np.zeros(self.tbbbb.shape)
        self.SS(vector)
        self.SD(vector)
        self.DS(vector)
        self.DD(vector)
        print(self.AAAA)
        print(self.ABAB)
        print(self.BABA)
        print(self.BBBB)
        return {'aa': self.AA, 'bb': self.BB, 'aaaa': self.AAAA, 'abab': self.ABAB, 'baba': self.BABA, 'bbbb': self.BBBB}

    def SS(self, vector):
        #HAB
        self.AA += .5*contract('ijab,ia->jb',self.Jaaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:], vector['aa'])
        self.BB += .5*contract('ijab,ia->jb',self.Labab[:self.NOa,:self.NOb,self.NOa:,self.NOb:], vector['aa'])
        self.BB += .5*contract('ijab,ia->jb',self.Jaaaa[:self.NOb,:self.NOb,self.NOb:,self.NOb:], vector['bb'])
        self.AA += .5*contract('ijab,ia->jb',self.Labab[:self.NOb,:self.NOa,self.NOb:,self.NOa:], vector['bb'])

        #BHA
        self.AA += contract('ba,ia->ib',self.fa[self.NOa:,self.NOa:], vector['aa'])
        self.BB += contract('ba,ia->ib',self.fa[self.NOb:,self.NOb:], vector['bb'])
        self.AA -= contract('ij,ia->ja',self.fa[:self.NOa,:self.NOa], vector['aa'])
        self.BB -= contract('ij,ia->ja',self.fa[:self.NOb,:self.NOb], vector['bb'])

        self.AA += contract('ibaj,ia->jb',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], vector['aa'])
        self.BB += contract('ibaj,ia->jb',self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['aa'])
        self.BB += contract('ibaj,ia->jb',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], vector['bb'])
        self.AA += contract('ibaj,ia->jb',self.Jbaba[:self.NOb,self.NOa:,self.NOb:,:self.NOa], vector['bb'])

    def SD(self, vector):
        #(Doubles Amplitude)
        self.AA -= .5*contract('ijak,ijab->kb',self.Laaaa[:self.NOa,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AA += .5*contract('ijka,ijba->kb',self.Jabab[:self.NOa,:self.NOb,:self.NOa,self.NOb:], vector['abab'])
        self.AA -= .5*contract('ijak,ijab->kb',self.Jbaba[:self.NOb,:self.NOa,self.NOb:,:self.NOa], vector['aaaa'])
        self.BB -= .5*contract('ijak,ijab->kb',self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
        self.BB += .5*contract('ijka,ijba->kb',self.Jbaba[:self.NOb,:self.NOa,:self.NOb,self.NOa:], vector['baba'])
        self.BB -= .5*contract('ijak,ijab->kb',self.Jabab[:self.NOa,:self.NOb,self.NOa:,:self.NOb], vector['bbbb'])

        self.AA += .5*contract('icab,ijab->jc',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,self.NOb:], vector['aaaa'])
        self.AA += .5*contract('icab,ijab->jc',self.Jbaba[:self.NOb,self.NOa:,self.NOb:,self.NOb:], vector['baba'])
        self.AA -= .5*contract('icba,ijba->jc',self.Jbaba[:self.NOb,self.NOa:,self.NOb:,self.NOb:], vector['baba'])
        self.BB += .5*contract('icab,ijab->jc',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,self.NOa:], vector['bbbb'])
        self.BB += .5*contract('icab,ijab->jc',self.Jabab[:self.NOa,self.NOb:,self.NOa:,self.NOa:], vector['abab'])
        self.BB -= .5*contract('icba,ijba->jc',self.Jabab[:self.NOa,self.NOb:,self.NOa:,self.NOa:], vector['abab'])

    def DS(self, vector):
        self.AAAA -= contract('ibkj,ia->kjab',self.Laaaa[:self.NOa,self.NOa:,:self.NOa,:self.NOa], vector['aa'])
        self.BABA -= contract('ibjk,ia->kjba',self.Jabab[:self.NOa,self.NOb:,:self.NOa,:self.NOb], vector['aa'])
        self.ABAB -= contract('ibkj,ia->kjab',self.Jabab[:self.NOa,self.NOb:,:self.NOa,:self.NOb], vector['aa'])
        self.AAAA += contract('iakj,ib->kjab',self.Laaaa[:self.NOa,self.NOa:,:self.NOa,:self.NOa], vector['aa'])
        self.BBBB -= contract('ibkj,ia->kjab',self.Lbbbb[:self.NOb,self.NOb:,:self.NOb,:self.NOb], vector['bb'])
        self.ABAB -= contract('ibjk,ia->kjba',self.Jbaba[:self.NOb,self.NOa:,:self.NOb,:self.NOa], vector['bb'])
        self.BABA -= contract('ibkj,ia->kjab',self.Jbaba[:self.NOb,self.NOa:,:self.NOb,:self.NOa], vector['bb'])
        self.BBBB += contract('iakj,ib->kjab',self.Lbbbb[:self.NOb,self.NOb:,:self.NOb,:self.NOb], vector['bb'])

        self.AAAA += contract('cbaj,ia->ijcb',self.Laaaa[self.NOa:,self.NOa:,self.NOa:,:self.NOa], vector['aa'])
        self.AAAA -= contract('cbai,ja->ijcb', self.Laaaa[self.NOa:, self.NOa:, self.NOa:, :self.NOa], vector['aa'])
        self.ABAB += contract('cbaj,ia->ijcb', self.Jabab[self.NOa:, self.NOb:, self.NOa:, :self.NOb], vector['aa'])
        self.ABAB += contract('bcaj,ia->ijbc', self.Jabab[self.NOa:, self.NOb:, self.NOa:, :self.NOb], vector['aa'])
        self.BBBB += contract('cbaj,ia->ijcb',self.Lbbbb[self.NOb:,self.NOb:,self.NOb:,:self.NOb], vector['bb'])
        self.BBBB -= contract('cbai,ja->ijcb', self.Lbbbb[self.NOb:, self.NOb:, self.NOb:, :self.NOb], vector['bb'])
        self.BABA += contract('cbaj,ia->ijcb', self.Jbaba[self.NOb:, self.NOa:, self.NOb:, :self.NOa], vector['bb'])
        self.BABA += contract('bcaj,ia->ijbc', self.Jbaba[self.NOb:, self.NOa:, self.NOb:, :self.NOa], vector['bb'])

    def DD(self,vector):
        self.AAAA += contract('ca,ijab->ijcb',self.fa[self.NOa:,self.NOa:],vector['aaaa'])
        self.AAAA -= contract('ba,ijba->ijcb',self.fa[self.NOa:,self.NOa:],vector['aaaa'])
        self.ABAB += contract('ca,jiab->jicb',self.fa[self.NOa:,self.NOa:],vector['abab'])
        self.ABAB += contract('ca,ijab->ijcb',self.fa[self.NOa:,self.NOa:],vector['abab'])

        self.AAAA += .5*contract('cdab,ijab->ijcd',self.Laaaa[self.NOa:,self.NOa:,self.NOa:,self.NOa:], vector['aaaa'])
        self.ABAB += .5*contract('cdab,ijab->ijcd',self.Jabab[self.NOa:,self.NOb:,self.NOa:,self.NOb:], vector['abab'])
        self.ABAB += .5*contract('cdba,ijba->ijcd',self.Jabab[self.NOa:,self.NOb:,self.NOa:,self.NOb:], vector['abab'])
        self.ABAB += .5*contract('dcab,ijab->ijdc',self.Jabab[self.NOa:,self.NOb:,self.NOa:,self.NOb:], vector['abab'])
        self.ABAB += .5*contract('dcba,ijba->ijdc',self.Jabab[self.NOa:,self.NOb:,self.NOa:,self.NOb:], vector['abab'])
        self.BBBB += .5*contract('cdab,ijab->ijcd',self.Lbbbb[self.NOb:,self.NOb:,self.NOb:,self.NOb:], vector['bbbb'])
        self.BABA += .5*contract('cdab,ijab->ijcd',self.Jbaba[self.NOb:,self.NOa:,self.NOb:,self.NOa:], vector['baba'])
        self.BABA += .5*contract('cdba,ijba->ijcd',self.Jbaba[self.NOb:,self.NOa:,self.NOb:,self.NOa:], vector['baba'])
        self.BABA += .5*contract('dcab,ijab->ijdc',self.Jbaba[self.NOb:,self.NOa:,self.NOb:,self.NOa:], vector['baba'])
        self.BABA += .5*contract('dcba,ijba->ijdc',self.Jbaba[self.NOb:,self.NOa:,self.NOb:,self.NOa:], vector['baba'])

        self.AAAA += .5*contract('ijkl,ijab->klab',self.Laaaa[:self.NOa,:self.NOa,:self.NOa,:self.NOa], vector['aaaa'])
        self.ABAB += .5*contract('ijkl,ijab->klab',self.Jabab[:self.NOa,:self.NOb,:self.NOa,:self.NOb], vector['abab'])
        self.ABAB += .5*contract('jikl,jiab->klab',self.Jabab[:self.NOa,:self.NOb,:self.NOa,:self.NOb], vector['abab'])
        self.ABAB += .5*contract('jilk,jiab->lkab',self.Jabab[:self.NOa,:self.NOb,:self.NOa,:self.NOb], vector['abab'])
        self.ABAB += .5*contract('ijlk,ijab->lkab',self.Jabab[:self.NOa,:self.NOb,:self.NOa,:self.NOb], vector['abab'])
        self.BBBB += .5*contract('ijkl,ijab->klab',self.Lbbbb[:self.NOb,:self.NOb,:self.NOb,:self.NOb], vector['bbbb'])
        self.BABA += .5*contract('ijkl,ijab->klab',self.Jbaba[:self.NOb,:self.NOa,:self.NOb,:self.NOa], vector['baba'])
        self.BABA += .5*contract('jikl,jiab->klab',self.Jbaba[:self.NOb,:self.NOa,:self.NOb,:self.NOa], vector['baba'])
        self.BABA += .5*contract('jilk,jiab->lkab',self.Jbaba[:self.NOb,:self.NOa,:self.NOb,:self.NOa], vector['baba'])
        self.BABA += .5*contract('ijlk,ijab->lkab',self.Jbaba[:self.NOb,:self.NOa,:self.NOb,:self.NOa], vector['baba'])

        self.AAAA -= contract('cjak,ijab->ikab',self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA += contract('ciak,kjab->ikab',self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA += contract('bjak,ijac->ikab',self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA -= contract('biak,kjac->ikab',self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.ABAB -= contract('cjak,ijab->ikab',self.Jabab[self.NOa:,:self.NOb,self.NOa:,:self.NOb], vector['abab'])
        self.BBBB -= contract('cjak,ijab->ikab',self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB += contract('ciak,kjab->ikab',self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB += contract('bjak,ijac->ikab',self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB -= contract('biak,kjac->ikab',self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
        self.BABA -= contract('cjak,ijab->ikab',self.Jbaba[self.NOb:,:self.NOa,self.NOb:,:self.NOa], vector['baba'])

        self.AAAA += contract('icak,ijab->kjcb',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA -= contract('ibak,ijac->kjcb',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA -= contract('icaj,ikab->kjcb',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA += contract('ibaj,ikac->kjcb',self.Laaaa[:self.NOa,self.NOa:,self.NOa:,:self.NOa], vector['aaaa'])
        self.ABAB -= contract('ciak,jiab->jkcb',self.Jabab[self.NOa:,:self.NOb,self.NOa:,:self.NOb], vector['abab'])
        self.BBBB += contract('icak,ijab->kjcb',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB -= contract('ibak,ijac->kjcb',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB -= contract('icaj,ikab->kjcb',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], vector['bbbb'])
        self.BBBB += contract('ibaj,ikac->kjcb',self.Lbbbb[:self.NOb,self.NOb:,self.NOb:,:self.NOb], vector['bbbb'])
        self.BABA -= contract('ciak,jiab->jkcb',self.Jbaba[self.NOb:,:self.NOa,self.NOb:,:self.NOa], vector['baba'])

    def CG(self):
        b = {'aa': self.gaa, 'bb': self.gbb, 'aaaa': self.gaaaa, 'abab': self.gabab, 'baba': self.gbaba, 'bbbb': self.gbbbb}
        trial = {'aa': self.taa, 'bb': self.tbb, 'aaaa': self.taaaa, 'abab':self.tabab, 'baba': self.tbaba, 'bbbb': self.gbbbb}
        x = trial
        Ax0 = self.Hessian_Action(x)
        r = TLAdd(Ax0,TLMult(b,-1))
        p = TLMult(r,-1)
        r_k_norm = TLDot(r,r)
        k = 0
        print('%5s|%20.16s' % (('Iter.', 'Residual Norm')))
        while r_k_norm > 1e-16:
            Ap = self.Hessian_Action(p)
            alpha = r_k_norm/TLDot(p,Ap)
            x = TLAdd(x, TLMult(p, alpha))
            r = TLAdd(r, TLMult(Ap, alpha))
            r_kplus1_norm = TLDot(r,r)
            beta = r_kplus1_norm/r_k_norm
            r_k_norm = r_kplus1_norm
            p = TLAdd(TLMult(p,beta),TLMult(r,-1))
            print('-'*30)
            k += 1
            print('{}'.format(k).ljust(5) + '|' + '%20.16f' % (r_k_norm))
        E = TLDot(b, x)+.5*TLDot(x, self.Hessian_Action(x))
        print('Converged energy:'.ljust(20) + '%20.16f Eh\n' % (E + self.hf_energy))


def TLDot(A,B):
    v = 0
    v += contract('ia,ia', A['aa'], B['aa'])
    v += contract('ia,ia', A['bb'], B['bb'])
    v += contract('ijab,ijab', A['aaaa'], B['aaaa'])
    v += contract('ijab,ijab', A['abab'], B['abab'])
    v += contract('ijab,ijab', A['baba'], B['baba'])
    v += contract('ijab,ijab', A['bbbb'], B['bbbb'])
    return v

def TLAdd(A,B):
    C = copy.copy(A)
    for c in C.keys():
        C[c]+=B[c]
    return C

def TLMult(A,a):
    C = copy.copy(A)
    for c in C.keys():
        C[c]*=a
    return C

if __name__ == '__main__':
    geometry = """
        0 1
        H 0 0 0
        H 0 0 1
        
    symmetry c1
    """
    basis = 'STO-3G'
    mol = Molecule(geometry, basis)
    trial = {'aa': mol.taa, 'bb': mol.tbb, 'aaaa': mol.taaaa, 'abab': mol.tabab, 'baba': mol.tbaba, 'bbbb': mol.tbbbb}
    mol.CG()
