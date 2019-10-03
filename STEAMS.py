#Input: Molecular data, basis, etc.
#Output: Second-order Taylor Expansion About Mean-field Solution

import numpy as np
#opt_einsum.contract() is just np.linalg.einsum() but more efficient
from opt_einsum import contract
import copy
import psi4

class Molecule:
    def __init__(self, geometry, basis, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        #compute and store all 1- and 2- electron integrals; currently memory-limiting
        molecule = psi4.geometry(geometry)
        psi4.core.set_output_file('output.dat', False)
        psi4.set_options({'basis': str(basis), 'scf_type': 'pk', 'reference': 'rhf', 'd_convergence': 1e-12})
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
        self.fb.transform(Cb)
        self.fa = np.array(self.fa)
        self.fb = np.array(self.fb)
        self.Jaaaa = self.Jaaaa.swapaxes(1, 2)
        self.Jabab = self.Jabab.swapaxes(1, 2)
        self.Jbaba = self.Jbaba.swapaxes(1, 2)
        self.Jbbbb = self.Jbbbb.swapaxes(1, 2)
        Kaaaa = copy.copy(self.Jaaaa)
        Kbbbb = copy.copy(self.Jbbbb)
        Kaaaa = Kaaaa.swapaxes(2, 3)
        Kbbbb = Kbbbb.swapaxes(2, 3)
        self.Laaaa = self.Jaaaa - Kaaaa
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
        self.AAAA = self.taaaa
        self.ABAB = self.tabab
        self.BB = self.tbb
        self.BBBB = self.tbbbb
        print(self.gaaaa)
        print(self.gabab)
        print(self.gbbbb)

    def Build_Trial(self):
        self.taa = np.zeros((self.NOa, self.NVa))
        self.taaaa = np.zeros((self.NOa, self.NOa, self.NVa, self.NVa))
        self.tabab = np.zeros((self.NOa, self.NOb, self.NVa, self.NVb))
        self.tbb = np.zeros((self.NOb, self.NVb))
        self.tbbbb = np.zeros((self.NOb, self.NOb, self.NVb, self.NVb))

    def Build_Gradient(self):
        self.gaa = np.zeros((self.NOa, self.NVa))
        self.gaaaa = 2**.5*self.Laaaa[:self.NOa,:self.NOa,self.NOa:,self.NOa:]
        self.gabab = 2**.5*self.Jabab[:self.NOa,:self.NOb,self.NOa:,self.NOb:]
        self.gbb = np.zeros((self.NOb, self.NVb))
        self.gbbbb = 2**.5*self.Lbbbb[:self.NOb,:self.NOb,self.NOb:,self.NOb:]

    def Hessian_Action(self, vector):
        self.AA = np.zeros(self.taa.shape)
        self.AAAA = np.zeros(self.taaaa.shape)
        self.ABAB = np.zeros(self.tabab.shape)
        self.BB = np.zeros(self.tbb.shape)
        self.BBBB = np.zeros(self.tbbbb.shape)
        self.SS(vector)
        self.SD(vector)
        self.DS(vector)
        self.DD(vector)
        return {'aa': self.AA, 'bb': self.BB, 'aaaa': self.AAAA, 'abab': self.ABAB, 'bbbb': self.BBBB}


    def SS(self, vector):
        #HBA
        #Particle Fock Hamiltonian
        self.AA += contract('ba,ia->ib', self.fa[self.NOa:, self.NOa:], vector['aa'])
        if self.RHF == True:
            self.BB = self.AA
        else:
            self.BB += contract('ba,ia->ib', self.fb[self.NOb:, self.NOb:], vector['bb'])

        #Hole Fock Hamiltonian
        self.AA -= contract('ij,ia->ja', self.fa[:self.NOa, :self.NOa], vector['aa'])
        if self.RHF == True:
            self.BB = self.AA
        else:
            self.BB -= contract('ij,ia->ja', self.fb[:self.NOb, :self.NOb], vector['bb'])

        #Ring Potential Hamiltonian
        self.AA += contract('ibaj,ia->jb', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, :self.NOa], vector['aa'])
        self.AA += contract('bija,ia->jb', self.Jabab[self.NOa:, :self.NOb, :self.NOa, self.NOb:], vector['bb'])
        if self.RHF == True:
            self.BB = self.AA
        else:
            self.BB += contract('ibaj,ia->jb', self.Jabab[:self.NOa, self.NOb:, self.NOa:, :self.NOb], vector['aa'])
            self.BB += contract('ibaj,ia->jb', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, :self.NOb], vector['bb'])

        if self.UNS == True:
            #Double Ring Hamiltonian
            self.AA += contract('ijab,ia->jb', self.Laaaa[:self.NOa, :self.NOa, self.NOa:, self.NOa:], vector['aa'])
            self.AA += contract('jiba,ia->jb', self.Jabab[:self.NOa, :self.NOb, self.NOa:, self.NOb:], vector['bb'])
            if self.RHF != True:
                self.BB += contract('ijab,ia->jb', self.Jabab[:self.NOa, :self.NOb, self.NOa:, self.NOb:], vector['aa'])
                self.BB += contract('ijab,ia->jb', self.Lbbbb[:self.NOb, :self.NOb, self.NOb:, self.NOb:], vector['bb'])
            else:
                self.BB = self.AA

    def SD(self, vector):
        #Particle Potential Hamiltonian
        self.AA += .5 * contract('icab,ijab->jc', self.Laaaa[:self.NOa, self.NOa:, self.NOa:, self.NOa:], vector['aaaa'])
        self.AA += .5 * contract('ciab,jiab->jc', self.Jabab[self.NOa:, :self.NOb, self.NOa:, self.NOb:], vector['abab'])
        self.AA += .5 * contract('ciba,jiba->jc', self.Jabab[self.NOa:, :self.NOb, self.NOa:, self.NOb:], vector['abab'])
        if self.RHF != True:
            self.BB += .5 * contract('icba,ijba->jc', self.Jabab[:self.NOa, self.NOb:, self.NOa:, self.NOb:], vector['abab'])
            self.BB += .5 * contract('icab,ijab->jc', self.Jabab[:self.NOa, self.NOb:, self.NOa:, self.NOb:], vector['abab'])
            self.BB += .5 * contract('icab,ijab->jc', self.Lbbbb[:self.NOb, self.NOb:, self.NOb:, self.NOb:], vector['bbbb'])
        else:
            self.BB = self.AA
        #Hole Potential Hamiltonian
        self.AA -= .5 * contract('ijak,ijab->kb', self.Laaaa[:self.NOa, :self.NOa, self.NOa:, :self.NOa], vector['aaaa'])
        self.AA -= .5 * contract('jika,jiba->kb', self.Jabab[:self.NOa, :self.NOb, :self.NOa, self.NOb:], vector['abab'])
        self.AA -= .5 * contract('ijka,ijba->kb', self.Jabab[:self.NOa, :self.NOb, :self.NOa, self.NOb:], vector['abab'])
        if self.RHF != True:
            self.BB -= .5 * contract('jiak,jiab->kb', self.Jabab[:self.NOa, :self.NOb, self.NOa:, :self.NOb], vector['abab'])
            self.BB -= .5 * contract('ijak,ijab->kb', self.Jabab[:self.NOa, :self.NOb, self.NOa:, :self.NOb], vector['abab'])
            self.BB -= .5 * contract('ijak,ijab->kb', self.Lbbbb[:self.NOb, :self.NOb, self.NOb:, :self.NOb], vector['bbbb'])
        else:
            self.BB = self.AA
    def DS(self, vector):
        #One of the AAAA terms here is the problem
        #Particle Potential Hamiltonian
        self.AAAA += contract('cbaj,ia->ijcb', self.Laaaa[self.NOa:, self.NOa:, self.NOa:, :self.NOa], vector['aa'])
        self.ABAB += contract('cbaj,ia->ijcb', self.Jabab[self.NOa:, self.NOb:, self.NOa:, :self.NOb], vector['aa'])
        self.AAAA -= contract('cbai,ja->ijcb', self.Laaaa[self.NOa:, self.NOa:, self.NOa:, :self.NOa], vector['aa'])
        self.ABAB += contract('bcja,ia->jibc', self.Jabab[self.NOa:, self.NOb:, :self.NOa, self.NOb:], vector['bb'])
        if self.RHF != True:
            self.BBBB += contract('cbaj,ia->ijcb', self.Lbbbb[self.NOb:, self.NOb:, self.NOb:, :self.NOb], vector['bb'])
            self.BBBB -= contract('cbai,ja->ijcb', self.Lbbbb[self.NOb:, self.NOb:, self.NOb:, :self.NOb], vector['bb'])
        else:
            self.BBBB = self.AAAA
        #Hole Potential Hamiltonian
        self.AAAA -= contract('ibkj,ia->kjab', self.Laaaa[:self.NOa, self.NOa:, :self.NOa, :self.NOa], vector['aa'])
        self.ABAB -= contract('ibkj,ia->kjab', self.Jabab[:self.NOa, self.NOb:, :self.NOa, :self.NOb], vector['aa'])
        self.AAAA += contract('iakj,ib->kjab', self.Laaaa[:self.NOa, self.NOa:, :self.NOa, :self.NOa], vector['aa'])
        self.ABAB -= contract('bijk,ia->jkba', self.Jabab[self.NOa:, :self.NOb, :self.NOa, :self.NOb], vector['bb'])
        if self.RHF != True:
            self.BBBB += contract('iakj,ib->kjab', self.Lbbbb[:self.NOb, self.NOb:, :self.NOb, :self.NOb], vector['bb'])
            self.BBBB -= contract('ibkj,ia->kjab', self.Lbbbb[:self.NOb, self.NOb:, :self.NOb, :self.NOb], vector['bb'])
        else:
            self.BBBB = self.AAAA


    def DD(self,vector):
        #Particle/Particle Potential Hamiltonian
        self.ABAB += .5*contract('dcba,jiba->jidc', self.Jabab[self.NOa:,self.NOb:,self.NOa:,self.NOb:], vector['abab'])
        self.ABAB += .5 * contract('dcab,ijab->ijdc', self.Jabab[self.NOa:, self.NOb:, self.NOa:, self.NOb:],
                                   vector['abab'])

        Oa = np.triu_indices(self.NOa, k = 1)
        Ob = np.triu_indices(self.NOb, k = 1)
        Va = np.triu_indices(self.NVa, k = 1)

        taaaa = vector['aaaa'][Oa][(slice(None),) + Va]

        Va = tuple([np.array(Va[0])+self.NOa, np.array(Va[1])+self.NOa])

        laaaa = self.Laaaa[Va][(slice(None),) + Va]

        raaaa = contract('ab, ib -> ia', laaaa, taaaa)

        Va = np.triu_indices(self.NVa, k = 1)

        self.AAAA[Oa[0][:, None], Oa[1][:, None], Va[0], Va[1]] += raaaa
        self.AAAA[Oa[0][:, None], Oa[1][:, None], Va[1], Va[0]] += -raaaa
        self.AAAA[Oa[1][:, None], Oa[0][:, None], Va[0], Va[1]] += -raaaa
        self.AAAA[Oa[1][:, None], Oa[0][:, None], Va[1], Va[0]] += raaaa
        if self.RHF != True:
            Vb = np.triu_indices(self.NVb, k=1)
            tbbbb = vector['bbbb'][Ob][(slice(None),) + Vb]
            Vb = tuple([np.array(Vb[0]) + self.NOb, np.array(Vb[1]) + self.NOb])
            lbbbb = self.Lbbbb[Vb][(slice(None),) + Vb]
            rbbbb = contract('ab, ib -> ia', lbbbb, tbbbb)
            Vb = np.triu_indices(self.NVb, k = 1)
            self.BBBB[Ob[0][:, None], Ob[1][:, None], Vb[0], Vb[1]] += rbbbb
            self.BBBB[Ob[0][:, None], Ob[1][:, None], Vb[1], Vb[0]] += -rbbbb
            self.BBBB[Ob[1][:, None], Ob[0][:, None], Vb[0], Vb[1]] += -rbbbb
            self.BBBB[Ob[1][:, None], Ob[0][:, None], Vb[1], Vb[0]] += rbbbb
        else:
            self.BBBB = self.AAAA

        #Hole/Hole Potential Hamiltonian
        self.AAAA += .5*contract('ijkl,ijab->klab', self.Laaaa[:self.NOa,:self.NOa,:self.NOa,:self.NOa], vector['aaaa'])
        self.ABAB += .5*contract('ijlk,ijab->lkab', self.Jabab[:self.NOa,:self.NOb,:self.NOa,:self.NOb], vector['abab'])
        self.ABAB += .5 * contract('jilk,jiab->lkab', self.Jabab[:self.NOa, :self.NOb, :self.NOa, :self.NOb],
                                   vector['abab'])
        if self.RHF != True:
            self.BBBB += .5*contract('ijkl,ijab->klab', self.Lbbbb[:self.NOb,:self.NOb,:self.NOb,:self.NOb], vector['bbbb'])
        else:
            self.BBBB = self.AAAA
        #Particle/Hole Potential Hamiltonian
        self.AAAA -= contract('cjak,ijab->ikcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA += contract('bjak,ijac->ikcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA += contract('cjai,kjab->ikcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA -= contract('bjai,kjac->ikcb', self.Laaaa[self.NOa:,:self.NOa,self.NOa:,:self.NOa], vector['aaaa'])
        self.AAAA -= contract('cjka,ijba->ikcb', self.Jabab[self.NOa:,:self.NOb,:self.NOa,self.NOb:], vector['abab'])
        self.AAAA += contract('cjia,kjba->ikcb', self.Jabab[self.NOa:,:self.NOb,:self.NOa,self.NOb:], vector['abab'])
        self.AAAA += contract('bjka,ijca->ikcb', self.Jabab[self.NOa:,:self.NOb,:self.NOa,self.NOb:], vector['abab'])
        self.AAAA -= contract('bjia,kjca->ikcb', self.Jabab[self.NOa:,:self.NOb,:self.NOa,self.NOb:], vector['abab'])
        self.ABAB -= contract('jcak,ijab->ikbc', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['aaaa'])
        self.ABAB -= contract('cjak,jiab->kicb', self.Laaaa[self.NOa:, :self.NOa, self.NOa:, :self.NOa], vector['abab'])
        self.ABAB -= contract('cjak,ijab->ikcb', self.Jabab[self.NOa:, :self.NOb, self.NOa:, :self.NOb], vector['abab'])
        self.ABAB -= contract('jcka,jiba->kibc', self.Jabab[:self.NOa, self.NOb:, :self.NOa, self.NOb:], vector['abab'])
        self.ABAB -= contract('cjak,ijba->ikbc', self.Lbbbb[self.NOb:, :self.NOb, self.NOb:, :self.NOb], vector['abab'])
        self.ABAB -= contract('cjka,ijab->kicb', self.Jabab[self.NOa:, :self.NOb, :self.NOa, self.NOb:], vector['bbbb'])
        if self.RHF != True:
            self.BBBB -= contract('cjak,ijab->ikcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
            self.BBBB += contract('bjak,ijac->ikcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
            self.BBBB += contract('cjai,kjab->ikcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
            self.BBBB -= contract('bjai,kjac->ikcb', self.Lbbbb[self.NOb:,:self.NOb,self.NOb:,:self.NOb], vector['bbbb'])
            self.BBBB -= contract('jcak,jiab->ikcb', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['abab'])
            self.BBBB += contract('jcai,jkab->ikcb', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['abab'])
            self.BBBB += contract('jbak,jiac->ikcb', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['abab'])
            self.BBBB -= contract('jbai,jkac->ikcb', self.Jabab[:self.NOa,self.NOb:,self.NOa:,:self.NOb], vector['abab'])
        else:
            self.AAAA = self.BBBB
        #Particle Fock Hamiltonian
        self.AAAA += contract('ca,ijab->ijcb', self.fa[self.NOa:,self.NOa:], vector['aaaa'])
        self.AAAA -= contract('ba,ijac->ijcb', self.fa[self.NOa:,self.NOa:], vector['aaaa'])
        self.ABAB += contract('ca,ijab->ijcb', self.fa[self.NOa:,self.NOa:], vector['abab'])
        self.ABAB += contract('ca,jiba->jibc', self.fb[self.NOb:, self.NOb:], vector['abab'])
        if self.RHF != True:
            self.BBBB += contract('ca,ijab->ijcb', self.fb[self.NOb:,self.NOb:], vector['bbbb'])
            self.BBBB -= contract('ba,ijac->ijcb', self.fb[self.NOb:,self.NOb:], vector['bbbb'])
        else:
            self.BBBB = self.AAAA

        #Hole Fock Hamiltonian
        self.AAAA -= contract('ik,ijab->kjab', self.fa[:self.NOa,:self.NOa], vector['aaaa'])
        self.AAAA += contract('ij,ikab->kjab', self.fa[:self.NOa,:self.NOa], vector['aaaa'])
        self.ABAB -= contract('ik,ijab->kjab', self.fa[:self.NOa,:self.NOa], vector['abab'])
        self.ABAB -= contract('ik,jiba->jkba', self.fb[:self.NOb, :self.NOb], vector['abab'])
        if self.RHF != True:
            self.BBBB -= contract('ik,ijab->kjab', self.fb[:self.NOb,:self.NOb], vector['bbbb'])
            self.BBBB += contract('ij,ikab->kjab', self.fb[:self.NOb,:self.NOb], vector['bbbb'])
        else:
            self.BBBB = self.AAAA

    def conj_grad(self):


        b = {'aa': -self.gaa, 'bb': -self.gbb, 'aaaa': -self.gaaaa, 'abab': -self.gabab, 'bbbb': -self.gbbbb}
        trial = {'aa': self.taa, 'bb': self.tbb, 'aaaa': self.taaaa, 'abab':self.tabab, 'bbbb': self.tbbbb}
        gradient = vec_lc(-1,b,0,b)

        x = trial
        ax0 = self.Hessian_Action(x)
        r = vec_lc(1,ax0,1,gradient)
        p = vec_lc(-1,r,0,r)
        r_k_norm = vec_dot(r,r)
        k = 0
        print('%5s|%20.16s' % (('Iter.', 'Residual Norm')))
        while r_k_norm > 1e-16:
            print('-'*30)
            print('{}'.format(k).ljust(5) + '|' + '%20.16f' % (r_k_norm))
            ap = self.Hessian_Action(p)
            alpha = r_k_norm/vec_dot(p,ap)
            palpha = vec_lc(alpha, p, 0, p)
            x = vec_lc(1, x, 1, palpha)
            apalpha = vec_lc(alpha, ap, 0, ap)
            r = vec_lc(1, r, 1, apalpha)
            r_kplus1_norm = vec_dot(r,r)
            beta = r_kplus1_norm/r_k_norm
            r_k_norm = r_kplus1_norm
            p = vec_lc(beta,p,-1,r)
            k += 1
        print('-' * 30)
        print('{}'.format(k).ljust(5) + '|' + '%20.16f' % (r_k_norm))
        E = vec_dot(gradient, x)+.5*vec_dot(x,self.Hessian_Action(x))
        print('Converged energy:'.ljust(20) + '%20.16f Eh\n' % (E + self.hf_energy))
        energy = E + self.hf_energy
        return energy

def vec_lc(scalar_1,tensor_1,scalar_2,tensor_2):
     vec = {}
     for key in tensor_1:
         vec[key] = tensor_1[key]*scalar_1+tensor_2[key]*scalar_2
     return vec

def vec_dot(A, B):
    dot = 0
    dot += contract('ij,ij', A['aa'],(B['aa']))
    dot += contract('ij,ij', A['bb'],(B['bb']))
    dot += contract('ijab,ijab', A['abab'],(B['abab']))
    dot += .25*contract('ijab,ijab', A['aaaa'],(B['aaaa']))
    dot += .25*contract('ijab,ijab', A['bbbb'],(B['bbbb']))
    return dot


if __name__ == '__main__':
    geometry = """
        0 2
        H 0 0 0 
        H 0 0 1 
        H 0 0 2 
        symmetry c1  
    """
    basis = 'sto-3g'
    mol = Molecule(geometry, basis, RHF = False, UNS = True)
    mol.conj_grad()
