import numpy as np
from opt_einsum import contract
import copy
import psi4

class molecule:
    """
        Returns molecule object.

        Contains tools for computing Newton step approximations to the energy of a chemical system.

        :geometry: Matrix for which sqare root is computed
        :basis: Basis for Psi4 to use
        :uns: Optional, defaults to False.  Whether to include ring diagram from T_1^2
        :reference: Optional, defaults to rhf.  Reference for SCF calculation
    """

    def __init__(self, geometry, basis, **kwargs):
        import psi4
        psi4.core.clean()
        self.uns = False
        self.reference = 'rhf'
        for key, value in kwargs.items():
            setattr(self, key, value)
        if self.reference == 'rhf':
            self.rhf = True
        else:
            self.rhf = False
        molecule = psi4.geometry(geometry)
        psi4.core.be_quiet()
        psi4.set_options({'basis': basis, 'd_convergence': 1e-12, 'scf_type': 'pk'})
        self.hf_energy, wfn = psi4.energy('scf', return_wfn=True)
        print("HF energy:".ljust(30)+("{0:20.16f}".format(self.hf_energy)))
        mints = psi4.core.MintsHelper(wfn.basisset())
        ca = wfn.Ca()
        cb = wfn.Cb()
        self.fa = wfn.Fa()
        self.fb = wfn.Fb()
        self.j_aaaa = np.array(mints.mo_eri(ca, ca, ca, ca))
        self.j_abab = np.array(mints.mo_eri(ca, ca, cb, cb))
        self.j_baba = np.array(mints.mo_eri(cb, cb, ca, ca))
        self.j_bbbb = np.array(mints.mo_eri(cb, cb, cb, cb))
        
        self.fa.transform(ca)
        if self.rhf == False:
            self.fb.transform(cb)
        self.fa = np.array(self.fa)
        self.fb = np.array(self.fb)
        self.j_aaaa = self.j_aaaa.swapaxes(1, 2)
        self.j_abab = self.j_abab.swapaxes(1, 2)
        self.j_baba = self.j_baba.swapaxes(1, 2)
        self.j_bbbb = self.j_bbbb.swapaxes(1, 2)
        k_aaaa = copy.copy(self.j_aaaa)
        k_bbbb = copy.copy(self.j_bbbb)
        k_aaaa = k_aaaa.swapaxes(2, 3)
        k_bbbb = k_bbbb.swapaxes(2, 3)
        self.l_aaaa = self.j_aaaa - k_aaaa
        self.l_bbbb = self.j_bbbb - k_bbbb
        self.noa = wfn.Ca_subset("AO", "ACTIVE_OCC").shape[1]
        self.nob = wfn.Cb_subset("AO", "ACTIVE_OCC").shape[1]
        self.nva = wfn.Ca_subset("AO", "ACTIVE_VIR").shape[1]
        self.nvb = wfn.Cb_subset("AO", "ACTIVE_VIR").shape[1]
        self.build_trial()
        self.build_gradient()
        self.r_aa = self.taa
        self.r_aaaa = self.taaaa
        self.r_abab = self.tabab
        self.r_bb = self.tbb
        self.r_bbbb = self.tbbbb

    def build_trial(self):
        """
            No return value.

            Initializes trial vector for CG solver
        """
        self.taa = np.zeros((self.noa, self.nva))
        self.taaaa = np.zeros((self.noa, self.noa, self.nva, self.nva))
        self.tabab = np.zeros((self.noa, self.nob, self.nva, self.nvb))
        self.tbb = np.zeros((self.nob, self.nvb))
        self.tbbbb = np.zeros((self.nob, self.nob, self.nvb, self.nvb))

    def build_gradient(self):
        """
            No return value.

            Computes gradient vector for CG solver
        """
        self.gaa = np.zeros((self.noa, self.nva))
        self.gaaaa = 2**.5*self.l_aaaa[:self.noa,:self.noa,self.noa:,self.noa:]
        self.gabab = 2**.5*self.j_abab[:self.noa,:self.nob,self.noa:,self.nob:]
        self.gbb = np.zeros((self.nob, self.nvb))
        self.gbbbb = 2**.5*self.l_bbbb[:self.nob,:self.nob,self.nob:,self.nob:]

    def hessian_action(self, vector):
        """
            Returns Hessian action on an arbitrary vector.

            Computes Hessian action r required for CG solver

            :vector:  Dictionary of tensors.
        """
        self.r_aa = np.zeros(self.taa.shape)
        self.r_aaaa = np.zeros(self.taaaa.shape)
        self.r_abab = np.zeros(self.tabab.shape)
        self.r_bb = np.zeros(self.tbb.shape)
        self.r_bbbb = np.zeros(self.tbbbb.shape)
        self.ss(vector)
        self.sd(vector)
        self.ds(vector)
        self.dd(vector)
        return {'aa': self.r_aa, 'bb': self.r_bb, 'aaaa': self.r_aaaa, 'abab': self.r_abab, 'bbbb': self.r_bbbb}

    def ss(self, vector):
        """
            Returns nothing.

            Computes effect of singles on T_1.

            :vector:  Dictionary of tensors.
        """
        #HBA
        #Particle Fock Hamiltonian
        self.r_aa += contract('ba,ia->ib', self.fa[self.noa:, self.noa:], vector['aa'])
        if self.rhf == True:
            self.r_bb = self.r_aa
        else:
            self.r_bb += contract('ba,ia->ib', self.fb[self.nob:, self.nob:], vector['bb'])

        #Hole Fock Hamiltonian
        self.r_aa -= contract('ij,ia->ja', self.fa[:self.noa, :self.noa], vector['aa'])
        if self.rhf == True:
            self.r_bb = self.r_aa
        else:
            self.r_bb -= contract('ij,ia->ja', self.fb[:self.nob, :self.nob], vector['bb'])

        #Ring Potential Hamiltonian
        self.r_aa += contract('ibaj,ia->jb', self.l_aaaa[:self.noa, self.noa:, self.noa:, :self.noa], vector['aa'])
        self.r_aa += contract('bija,ia->jb', self.j_abab[self.noa:, :self.nob, :self.noa, self.nob:], vector['bb'])
        if self.rhf == True:
            self.r_bb = self.r_aa
        else:
            self.r_bb += contract('ibaj,ia->jb', self.j_abab[:self.noa, self.nob:, self.noa:, :self.nob], vector['aa'])
            self.r_bb += contract('ibaj,ia->jb', self.l_bbbb[:self.nob, self.nob:, self.nob:, :self.nob], vector['bb'])

        if self.uns == True:
            #Double Ring Hamiltonian
            self.r_aa += contract('ijab,ia->jb', self.l_aaaa[:self.noa, :self.noa, self.noa:, self.noa:], vector['aa'])
            self.r_aa += contract('jiba,ia->jb', self.j_abab[:self.noa, :self.nob, self.noa:, self.nob:], vector['bb'])
            if self.rhf != True:
                self.r_bb += contract('ijab,ia->jb', self.j_abab[:self.noa, :self.nob, self.noa:, self.nob:], vector['aa'])
                self.r_bb += contract('ijab,ia->jb', self.l_bbbb[:self.nob, :self.nob, self.nob:, self.nob:], vector['bb'])
            else:
                self.r_bb = self.r_aa

    def sd(self, vector):
        """
            Returns nothing.

            Computes effect of singles on T_2.

            :vector:  Dictionary of tensors.
        """
        #Particle Potential Hamiltonian
        self.r_aa += .5 * contract('icab,ijab->jc', self.l_aaaa[:self.noa, self.noa:, self.noa:, self.noa:], vector['aaaa'])
        self.r_aa += .5 * contract('ciab,jiab->jc', self.j_abab[self.noa:, :self.nob, self.noa:, self.nob:], vector['abab'])
        self.r_aa += .5 * contract('ciba,jiba->jc', self.j_abab[self.noa:, :self.nob, self.noa:, self.nob:], vector['abab'])
        if self.rhf != True:
            self.r_bb += .5 * contract('icba,ijba->jc', self.j_abab[:self.noa, self.nob:, self.noa:, self.nob:], vector['abab'])
            self.r_bb += .5 * contract('icab,ijab->jc', self.j_abab[:self.noa, self.nob:, self.noa:, self.nob:], vector['abab'])
            self.r_bb += .5 * contract('icab,ijab->jc', self.l_bbbb[:self.nob, self.nob:, self.nob:, self.nob:], vector['bbbb'])
        else:
            self.r_bb = self.r_aa
        #Hole Potential Hamiltonian
        self.r_aa -= .5 * contract('ijak,ijab->kb', self.l_aaaa[:self.noa, :self.noa, self.noa:, :self.noa], vector['aaaa'])
        self.r_aa -= .5 * contract('jika,jiba->kb', self.j_abab[:self.noa, :self.nob, :self.noa, self.nob:], vector['abab'])
        self.r_aa -= .5 * contract('ijka,ijba->kb', self.j_abab[:self.noa, :self.nob, :self.noa, self.nob:], vector['abab'])
        if self.rhf != True:
            self.r_bb -= .5 * contract('jiak,jiab->kb', self.j_abab[:self.noa, :self.nob, self.noa:, :self.nob], vector['abab'])
            self.r_bb -= .5 * contract('ijak,ijab->kb', self.j_abab[:self.noa, :self.nob, self.noa:, :self.nob], vector['abab'])
            self.r_bb -= .5 * contract('ijak,ijab->kb', self.l_bbbb[:self.nob, :self.nob, self.nob:, :self.nob], vector['bbbb'])
        else:
            self.r_bb = self.r_aa

    def ds(self, vector):
        """
            Returns nothing.

            Computes effect of doubles on T_1.

            :vector:  Dictionary of tensors.
        """
        #Particle Potential Hamiltonian
        self.r_aaaa += contract('cbaj,ia->ijcb', self.l_aaaa[self.noa:, self.noa:, self.noa:, :self.noa], vector['aa'])
        self.r_abab += contract('cbaj,ia->ijcb', self.j_abab[self.noa:, self.nob:, self.noa:, :self.nob], vector['aa'])
        self.r_aaaa -= contract('cbai,ja->ijcb', self.l_aaaa[self.noa:, self.noa:, self.noa:, :self.noa], vector['aa'])
        self.r_abab += contract('bcja,ia->jibc', self.j_abab[self.noa:, self.nob:, :self.noa, self.nob:], vector['bb'])
        if self.rhf != True:
            self.r_bbbb += contract('cbaj,ia->ijcb', self.l_bbbb[self.nob:, self.nob:, self.nob:, :self.nob], vector['bb'])
            self.r_bbbb -= contract('cbai,ja->ijcb', self.l_bbbb[self.nob:, self.nob:, self.nob:, :self.nob], vector['bb'])
        else:
            self.r_bbbb = self.r_aaaa
        #Hole Potential Hamiltonian
        self.r_aaaa -= contract('ibkj,ia->kjab', self.l_aaaa[:self.noa, self.noa:, :self.noa, :self.noa], vector['aa'])
        self.r_abab -= contract('ibkj,ia->kjab', self.j_abab[:self.noa, self.nob:, :self.noa, :self.nob], vector['aa'])
        self.r_aaaa += contract('iakj,ib->kjab', self.l_aaaa[:self.noa, self.noa:, :self.noa, :self.noa], vector['aa'])
        self.r_abab -= contract('bijk,ia->jkba', self.j_abab[self.noa:, :self.nob, :self.noa, :self.nob], vector['bb'])
        if self.rhf != True:
            self.r_bbbb += contract('iakj,ib->kjab', self.l_bbbb[:self.nob, self.nob:, :self.nob, :self.nob], vector['bb'])
            self.r_bbbb -= contract('ibkj,ia->kjab', self.l_bbbb[:self.nob, self.nob:, :self.nob, :self.nob], vector['bb'])
        else:
            self.r_bbbb = self.r_aaaa


    def dd(self,vector):
        """
            Returns nothing.

            Computes effect of doubles on T_2.

            :vector:  Dictionary of tensors.
        """
        #Particle/Particle Potential Hamiltonian
        self.r_abab += .5 * contract('dcba,jiba->jidc', self.j_abab[self.noa:, self.nob:, self.noa:, self.nob:], vector['abab'])
        self.r_abab += .5 * contract('dcab,ijab->ijdc', self.j_abab[self.noa:, self.nob:, self.noa:, self.nob:], vector['abab'])
        oa = np.triu_indices(self.noa, k = 1)
        ob = np.triu_indices(self.nob, k = 1)
        va = np.triu_indices(self.nva, k = 1)
        taaaa = vector['aaaa'][oa][(slice(None),) + va]
        va = tuple([np.array(va[0])+self.noa, np.array(va[1])+self.noa])
        laaaa = self.l_aaaa[va][(slice(None),) + va]
        raaaa = contract('ab, ib -> ia', laaaa, taaaa)
        va = np.triu_indices(self.nva, k = 1)
        self.r_aaaa[oa[0][:, None], oa[1][:, None], va[0], va[1]] += raaaa
        self.r_aaaa[oa[0][:, None], oa[1][:, None], va[1], va[0]] += -raaaa
        self.r_aaaa[oa[1][:, None], oa[0][:, None], va[0], va[1]] += -raaaa
        self.r_aaaa[oa[1][:, None], oa[0][:, None], va[1], va[0]] += raaaa
        if self.rhf != True:
            vb = np.triu_indices(self.nvb, k=1)
            tbbbb = vector['bbbb'][ob][(slice(None),) + vb]
            vb = tuple([np.array(vb[0]) + self.nob, np.array(vb[1]) + self.nob])
            lbbbb = self.l_bbbb[vb][(slice(None),) + vb]
            rbbbb = contract('ab, ib -> ia', lbbbb, tbbbb)
            vb = np.triu_indices(self.nvb, k = 1)
            self.r_bbbb[ob[0][:, None], ob[1][:, None], vb[0], vb[1]] += rbbbb
            self.r_bbbb[ob[0][:, None], ob[1][:, None], vb[1], vb[0]] += -rbbbb
            self.r_bbbb[ob[1][:, None], ob[0][:, None], vb[0], vb[1]] += -rbbbb
            self.r_bbbb[ob[1][:, None], ob[0][:, None], vb[1], vb[0]] += rbbbb
        else:
            self.r_bbbb = self.r_aaaa
        #Hole/Hole Potential Hamiltonian
        self.r_aaaa += .5 * contract('ijkl,ijab->klab', self.l_aaaa[:self.noa, :self.noa, :self.noa, :self.noa], vector['aaaa'])
        self.r_abab += .5 * contract('ijlk,ijab->lkab', self.j_abab[:self.noa, :self.nob, :self.noa, :self.nob], vector['abab'])
        self.r_abab += .5 * contract('jilk,jiab->lkab', self.j_abab[:self.noa, :self.nob, :self.noa, :self.nob], vector['abab'])
        if self.rhf != True:
            self.r_bbbb += .5*contract('ijkl,ijab->klab', self.l_bbbb[:self.nob,:self.nob,:self.nob,:self.nob], vector['bbbb'])
        else:
            self.r_bbbb = self.r_aaaa
        #Particle/Hole Potential Hamiltonian
        self.r_aaaa -= contract('cjak,ijab->ikcb', self.l_aaaa[self.noa:,:self.noa,self.noa:,:self.noa], vector['aaaa'])
        self.r_aaaa += contract('bjak,ijac->ikcb', self.l_aaaa[self.noa:,:self.noa,self.noa:,:self.noa], vector['aaaa'])
        self.r_aaaa += contract('cjai,kjab->ikcb', self.l_aaaa[self.noa:,:self.noa,self.noa:,:self.noa], vector['aaaa'])
        self.r_aaaa -= contract('bjai,kjac->ikcb', self.l_aaaa[self.noa:,:self.noa,self.noa:,:self.noa], vector['aaaa'])
        self.r_aaaa -= contract('cjka,ijba->ikcb', self.j_abab[self.noa:,:self.nob,:self.noa,self.nob:], vector['abab'])
        self.r_aaaa += contract('cjia,kjba->ikcb', self.j_abab[self.noa:,:self.nob,:self.noa,self.nob:], vector['abab'])
        self.r_aaaa += contract('bjka,ijca->ikcb', self.j_abab[self.noa:,:self.nob,:self.noa,self.nob:], vector['abab'])
        self.r_aaaa -= contract('bjia,kjca->ikcb', self.j_abab[self.noa:,:self.nob,:self.noa,self.nob:], vector['abab'])
        self.r_abab -= contract('jcak,ijab->ikbc', self.j_abab[:self.noa,self.nob:,self.noa:,:self.nob], vector['aaaa'])
        self.r_abab -= contract('cjak,jiab->kicb', self.l_aaaa[self.noa:, :self.noa, self.noa:, :self.noa], vector['abab'])
        self.r_abab -= contract('cjak,ijab->ikcb', self.j_abab[self.noa:, :self.nob, self.noa:, :self.nob], vector['abab'])
        self.r_abab -= contract('jcka,jiba->kibc', self.j_abab[:self.noa, self.nob:, :self.noa, self.nob:], vector['abab'])
        self.r_abab -= contract('cjak,ijba->ikbc', self.l_bbbb[self.nob:, :self.nob, self.nob:, :self.nob], vector['abab'])
        self.r_abab -= contract('cjka,ijab->kicb', self.j_abab[self.noa:, :self.nob, :self.noa, self.nob:], vector['bbbb'])
        if self.rhf != True:
            self.r_bbbb -= contract('cjak,ijab->ikcb', self.l_bbbb[self.nob:,:self.nob,self.nob:,:self.nob], vector['bbbb'])
            self.r_bbbb += contract('bjak,ijac->ikcb', self.l_bbbb[self.nob:,:self.nob,self.nob:,:self.nob], vector['bbbb'])
            self.r_bbbb += contract('cjai,kjab->ikcb', self.l_bbbb[self.nob:,:self.nob,self.nob:,:self.nob], vector['bbbb'])
            self.r_bbbb -= contract('bjai,kjac->ikcb', self.l_bbbb[self.nob:,:self.nob,self.nob:,:self.nob], vector['bbbb'])
            self.r_bbbb -= contract('jcak,jiab->ikcb', self.j_abab[:self.noa,self.nob:,self.noa:,:self.nob], vector['abab'])
            self.r_bbbb += contract('jcai,jkab->ikcb', self.j_abab[:self.noa,self.nob:,self.noa:,:self.nob], vector['abab'])
            self.r_bbbb += contract('jbak,jiac->ikcb', self.j_abab[:self.noa,self.nob:,self.noa:,:self.nob], vector['abab'])
            self.r_bbbb -= contract('jbai,jkac->ikcb', self.j_abab[:self.noa,self.nob:,self.noa:,:self.nob], vector['abab'])
        else:
            self.r_aaaa = self.r_bbbb
        #Particle Fock Hamiltonian
        self.r_aaaa += contract('ca,ijab->ijcb', self.fa[self.noa:,self.noa:], vector['aaaa'])
        self.r_aaaa -= contract('ba,ijac->ijcb', self.fa[self.noa:,self.noa:], vector['aaaa'])
        self.r_abab += contract('ca,ijab->ijcb', self.fa[self.noa:,self.noa:], vector['abab'])
        self.r_abab += contract('ca,jiba->jibc', self.fb[self.nob:, self.nob:], vector['abab'])
        if self.rhf != True:
            self.r_bbbb += contract('ca,ijab->ijcb', self.fb[self.nob:,self.nob:], vector['bbbb'])
            self.r_bbbb -= contract('ba,ijac->ijcb', self.fb[self.nob:,self.nob:], vector['bbbb'])
        else:
            self.r_bbbb = self.r_aaaa

        #Hole Fock Hamiltonian
        self.r_aaaa -= contract('ik,ijab->kjab', self.fa[:self.noa,:self.noa], vector['aaaa'])
        self.r_aaaa += contract('ij,ikab->kjab', self.fa[:self.noa,:self.noa], vector['aaaa'])
        self.r_abab -= contract('ik,ijab->kjab', self.fa[:self.noa,:self.noa], vector['abab'])
        self.r_abab -= contract('ik,jiba->jkba', self.fb[:self.nob, :self.nob], vector['abab'])
        if self.rhf != True:
            self.r_bbbb -= contract('ik,ijab->kjab', self.fb[:self.nob,:self.nob], vector['bbbb'])
            self.r_bbbb += contract('ij,ikab->kjab', self.fb[:self.nob,:self.nob], vector['bbbb'])
        else:
            self.r_bbbb = self.r_aaaa

    def conj_grad(self):
        """
            Returns energy of the system.

            Uses a conjugate gradient approach to iteratively solve the second-order Taylor energy expansion.

            :vector:  Dictionary of tensors.
        """
        b = {'aa': -self.gaa, 'bb': -self.gbb, 'aaaa': -self.gaaaa, 'abab': -self.gabab, 'bbbb': -self.gbbbb}
        trial = {'aa': self.taa, 'bb': self.tbb, 'aaaa': self.taaaa, 'abab':self.tabab, 'bbbb': self.tbbbb}
        gradient = vec_lc(-1,b,0,b)
        x = trial
        N = self.noa+self.nob
        Ec = 0
        shift = 0
        ax0 = self.hessian_action(x)
        ax0 = vec_lc(1, ax0, shift, x)
        r = vec_lc(1,ax0,1,gradient)
        p = vec_lc(-1,r,0,r)
        r_k_norm = vec_dot(r,r)
        k = 0
        print('Conjugate Gradient Tracking:')
        while r_k_norm > 1e-16:
            ap = self.hessian_action(p)
            ap = vec_lc(1, ap, shift, p)
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
            energy = self.hf_energy + vec_dot(gradient, x)+.5*vec_dot(x,vec_lc(1, self.hessian_action(x), shift, x))
            Ec = energy-self.hf_energy
            shift = 0
            print('Iter. '+str(k)+': '+str(r_k_norm)+'|E: '+str(energy))

        print("UNS energy:".ljust(30) + ("{0:20.16f}".format(energy)))
        return energy

def vec_lc(scalar_1,tensor_1,scalar_2,tensor_2):
    """
        Returns scalar_1*tensor_1+scalar_2*tensor_2.

        QOL function for manipulating dictionaries of tensors

        :scalar_1:
        :scalar_2:
        :tensor_1:
        :tensor_2:
    """
    vec = {}
    for key in tensor_1:
         vec[key] = tensor_1[key]*scalar_1+tensor_2[key]*scalar_2
    return vec

def vec_dot(v1, v2):
    """
        Returns dot product of v1, v2..

        QOL function for manipulating dictionaries of tensors

        :v1:
        :v2:
    """
    dot = 0
    dot += contract('ij,ij', v1['aa'],(v2['aa']))
    dot += contract('ij,ij', v1['bb'],(v2['bb']))
    dot += contract('ijab,ijab', v1['abab'],(v2['abab']))
    dot += .25*contract('ijab,ijab', v1['aaaa'],(v2['aaaa']))
    dot += .25*contract('ijab,ijab', v1['bbbb'],(v2['bbbb']))
    return dot


if __name__ == '__main__':
    geometry = """
        0 1
        H 0 0 0
        Cl 0 0 1
        symmetry c1
    """
    basis = 'cc-pvdz'
    mol = molecule(geometry, basis, reference = 'rhf', uns = False)
    mol.conj_grad()
