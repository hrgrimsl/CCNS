import numpy as np
from opt_einsum import contract
import psi4
import copy

import random
from sys import argv
from timeit import default_timer as timer

def get_gradient( n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    gradient = []
    #doubles
    for i in range(0, n_occs):
        for j in range(i+1, n_occs):
            for a in range(n_occs, n_orbitals):
                for b in range(a+1, n_orbitals):
                    gradient.append(2**.5*(g[i,j,b,a]-g[i,j,a,b]))
    #singles
    for i in range(0, n_occs):
        for a in range(n_occs, n_orbitals):
            gradient.append(0)

    gradient = np.array(gradient)
    return gradient

def Hessian_Action(dx, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    dx = Populate_Tensor(dx, n_occs, n_noccs, n_orbitals, f, g, hf_energy)    
    gvvvv = g[n_occs:,n_occs:,n_occs:,n_occs:]
    goooo = g[:n_occs,:n_occs,:n_occs,:n_occs]
    goovv = g[:n_occs,:n_occs,n_occs:,n_occs:]
    gvvoo = g[n_occs:,n_occs:,:n_occs,:n_occs]
    gvovo = g[n_occs:,:n_occs,n_occs:,:n_occs]
    gvoov = g[n_occs:,:n_occs,:n_occs,n_occs:]
    govov = g[:n_occs,n_occs:,:n_occs,n_occs:]
    govvo = g[:n_occs,n_occs:,n_occs:,:n_occs]
    gvvov = g[n_occs:,n_occs:,:n_occs,n_occs:]
    gvvvo = g[n_occs:,n_occs:,n_occs:,:n_occs]
    govoo = g[:n_occs,n_occs:,:n_occs,:n_occs]
    govvv = g[:n_occs,n_occs:,n_occs:,n_occs:]
    gvooo = g[n_occs:,:n_occs,:n_occs,:n_occs]
    goovo = g[:n_occs,:n_occs,n_occs:,:n_occs]
    gooov = g[:n_occs,:n_occs,:n_occs,n_occs:]
    foo = f[:n_occs,:n_occs]
    fvv = f[n_occs:,n_occs:]

    Hdx = 0
    Hdx2 = 0
    
    #Singles/Singles
    
    Hdx -= contract('ji,bi->bj',foo,dx[0])
    Hdx += contract('ba,aj->bj',fvv,dx[0]) 

    Hdx += contract('ibaj,ai->bj',govvo,dx[0])
    Hdx -= contract('ibja,ai->bj',govov,dx[0])
    
    Hdx += contract('ijab,ai->bj',goovv,dx[0])
    Hdx -= contract('ijba,ai->bj',goovv,dx[0])
    

    #Singles/Doubles - Out of Order

    Hdx2 += contract('cbaj,ak->bckj',gvvvo,dx[0])
    Hdx2 -= contract('cbja,ak->bckj',gvvov,dx[0])
    Hdx2 -= contract('cbak,aj->bckj',gvvvo,dx[0])
    Hdx2 += contract('cbka,aj->bckj',gvvov,dx[0])

    Hdx2 -= contract('ibkj,ci->bckj',govoo,dx[0])
    Hdx2 += contract('ibjk,ci->bckj',govoo,dx[0])
    Hdx2 += contract('ickj,bi->bckj',govoo,dx[0])
    Hdx2 -= contract('icjk,bi->bckj',govoo,dx[0]) 

    
    #Doubles/Singles
    Hdx -= .5*contract('ibac,acij->bj',govvv,dx[1])
    Hdx += .5*contract('ibca,acij->bj',govvv,dx[1])

    Hdx += .5*contract('ikaj,abik->bj',goovo,dx[1])
    Hdx -= .5*contract('ikja,abik->bj',gooov,dx[1])
    

    #Doubles/Doubles
    Hdx2 += contract('ab,acjk->bcjk',fvv,dx[1]) 
    Hdx2 -= contract('ac,abjk->bcjk',fvv,dx[1]) 

    Hdx2 -= contract('ij,bcik->bcjk',foo,dx[1]) 
    Hdx2 += contract('ik,bcij->bcjk',foo,dx[1])
 
    Hdx2 += .5*contract('bcad,adjk->bcjk',gvvvv,dx[1])
    Hdx2 -= .5*contract('bcda,adjk->bcjk',gvvvv,dx[1])
 
    Hdx2 += .5*contract('iljk,adil->adjk',goooo,dx[1])
    Hdx2 -= .5*contract('ilkj,adil->adjk',goooo,dx[1])
    
    Hdx2 -= contract('icjd,bdik->cbkj',govov,dx[1])
    Hdx2 += contract('icdj,bdik->cbkj',govvo,dx[1])
    Hdx2 += contract('ickd,bdij->cbkj',govov,dx[1])
    Hdx2 -= contract('icdk,bdij->cbkj',govvo,dx[1])
    Hdx2 += contract('ibjd,cdik->cbkj',govov,dx[1])
    Hdx2 -= contract('ibdj,cdik->cbkj',govvo,dx[1])
    Hdx2 -= contract('ibkd,cdij->cbkj',govov,dx[1])
    Hdx2 += contract('ibdk,cdij->cbkj',govvo,dx[1])

    return(Populate_Vector([Hdx2, Hdx], n_occs, n_noccs, n_orbitals, f, g, hf_energy))


def Populate_Tensor(t_vec, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    t_tensor = [np.zeros((n_noccs, n_occs)), np.zeros((n_noccs, n_noccs, n_occs, n_occs))]
    ind = 0
    for i in range(0, n_occs):
        for j in range(i+1, n_occs):
            for a in range(0, n_noccs):
                for b in range(a+1, n_noccs):
                    t_tensor[1][a,b,i,j] = t_vec[ind]
                    t_tensor[1][a,b,j,i] = -t_vec[ind]
                    t_tensor[1][b,a,i,j] = -t_vec[ind]
                    t_tensor[1][b,a,j,i] = t_vec[ind]
                    ind += 1

    for a in range(0, n_noccs):
        for i in range(0, n_occs):
            t_tensor[0][a,i] = t_vec[ind]
            ind += 1             

    return t_tensor

def Populate_Vector(t_tensor, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
   t_vector = []
   for i in range(0, n_occs):
         for j in range(i+1, n_occs):
             for a in range(0, n_noccs):
                 for b in range(a+1, n_noccs):
                     t_vector.append(t_tensor[0][a,b,i,j])
 
   for a in range(0, n_noccs):
       for i in range(0, n_occs):
             t_vector.append(t_tensor[1][a,i])

   return t_vector




def CG_Solver(t_vec, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
     k = 0
     b = -get_gradient(n_occs, n_noccs, n_orbitals, f, g, hf_energy)
     x = [t_vec]
     Ax0 = np.array(Hessian_Action(x[0], n_occs, n_noccs, n_orbitals, f, g, hf_energy))
     r = Ax0-b
     p = -r
     r_k_norm = np.dot(r,r)
     print('%5s|%20.16s|%20.16s'%(('Iter.', 'Residual Norm', 'Energy (a.u.)')))
     while r_k_norm > 1e-10:

         Ap = np.array(Hessian_Action(p, n_occs, n_noccs, n_orbitals, f, g, hf_energy))

         alpha = r_k_norm/np.dot(p,Ap)
         x += alpha * p
         r += alpha * Ap
         r_kplus1_norm = np.dot(r,r)
         beta = r_kplus1_norm/r_k_norm
         r_k_norm = r_kplus1_norm
         p = beta * p - r
         print('-'*48)
         k += 1
         print('{}'.format(k).ljust(5)+'|'+'%20.16f|%20.16f'%(r_k_norm,hf_energy+Energy(x, b, n_occs, n_noccs, n_orbitals, f, g, hf_energy)))

     end = timer()
     print('\nConverged in {} iterations  ({} seconds).'.format(k, end-start)+'\n')
     print('HF energy:'.ljust(20)+'%20.16f Eh\n'%(hf_energy))
     print('Converged energy:'.ljust(20)+'%20.16f Eh\n'%(Energy(x, b, n_occs, n_noccs, n_orbitals, f, g, hf_energy)+hf_energy))



def Energy(x, gradient, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    x = x.T
    E = (np.array(gradient).dot(x)+.5*x.T.dot(Hessian_Action(x, n_occs, n_noccs, n_orbitals, f, g, hf_energy)))
    return E

if __name__ == '__main__':
    start = timer()

    print('\n'+'- '*24+'\n')
    print('  '*10+'STARLIGHT'+' '*10+'\n')
    print(' '*4+'An H.R. Grimsley & N.J. Mayhall Algorithm\n')
    print('- '*24+'\n')
    psi4.core.set_output_file('output.dat', False)
    geometry = '''
    H 0 0 0
    Cl 0 0 1
    symmetry c1
    '''
    
    molecule = psi4.geometry(geometry)
    #Psi4 Calculations
    psi4.set_options({'basis': 'cc-pvdz', 'scf_type': 'pk'})

    hf_energy, wfn = psi4.energy('scf', return_wfn = True, scf_type = 'pk')

    #h = 1-electron integral 2-tensor (<i|h|a>) 
    molecule = psi4.geometry(geometry)
    #Psi4 Calculations
    psi4.set_options({'basis': 'cc-pvdz', 'scf_type': 'pk'})
    
    hf_energy, wfn = psi4.energy('scf', return_wfn = True, scf_type = 'pk')
    
    #h = 1-electron integral 2-tensor (<i|h|a>)
    #g = 2-electron integral 4-tensor (<ij|ab>, not <ij||ab>)
    #f = Fock matrix 2-tensor (<i|h|i>+\sum_a(<ii|aa>))
    #All tensors in MO, spinorbital representation
    
    nr = molecule.nuclear_repulsion_energy()
    mints = psi4.core.MintsHelper(wfn.basisset())
    C = wfn.Ca()
    h = mints.ao_kinetic()
    f = wfn.Fa()
    h.add(mints.ao_potential())
    h.transform(C)
    f.transform(C)
    h = np.kron(np.array(h),(np.array([[1,0],[0,1]])))
    g = np.asarray(mints.mo_spin_eri(C, C))
    g = g.swapaxes(1, 2)
    f = np.array(f)
    f = np.kron(f, (np.array([[1,0],[0,1]])))
    n_orbitals = 2*wfn.nmo()
    n_occs = 2*wfn.doccpi()[0]
    n_noccs = n_orbitals-n_occs
    occs = [i for i in range(0, n_occs)]
    noccs = [i for i in range(n_occs, n_orbitals)]
    n_ops = 0
    for i in range(0, n_occs):
         for j in range(i+1, n_occs):
             for a in range(0, n_noccs):
                 for b in range(a+1, n_noccs):
                         n_ops += 1
 
    for i in range(0, n_occs):
        for a in range(0, n_noccs):
            n_ops += 1

    t_vec = [0 for i in range(1, n_ops+1)]
    CG_Solver(t_vec, n_occs, n_noccs, n_orbitals, f, g, hf_energy)

