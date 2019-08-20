import numpy as np
from opt_einsum import contract
import psi4
import copy
import random
from sys import argv
from timeit import default_timer as timer

def get_gradient(n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    gradient = []
    #doubles
    for a in range(n_occs, n_orbitals):
        for b in range(a+1, n_orbitals):
            for i in range(0, n_occs):
               for j in range(i+1, n_occs):
                    gradient.append(2**.5*(g[i,j,b,a]))
    #singles
    for a in range(n_occs, n_orbitals):
        for i in range(0, n_occs):
            gradient.append(0)

    gradient = np.array(gradient)
    return gradient

def Hessian_Action(dx, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    t_vec = copy.copy(dx)
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

    Hdx = np.zeros((n_noccs,n_occs))
    Hdx2 = np.zeros((n_noccs,n_noccs,n_occs,n_occs))
    
    #Singles/Singles
    #BHA
    Hdx -= contract('ji,bi->bj',foo,dx[0])
    Hdx += contract('ba,aj->bj',fvv,dx[0]) 

    Hdx += contract('ibaj,ai->bj',govvo,dx[0])
    

    #Remove this term for CEPA(0)?
    #HAB
    Hdx += contract('ijab,ai->bj',goovv,dx[0])

    #Singles/Doubles
    #BHA
    Hdx2 += contract('cbaj,ak->bckj',gvvvo,dx[0])
    Hdx2 -= contract('cbak,aj->bckj',gvvvo,dx[0])

    Hdx2 -= contract('ibkj,ci->bckj',govoo,dx[0])
    Hdx2 += contract('ickj,bi->bckj',govoo,dx[0])
    
    #Doubles/Singles
    #BHA
    Hdx -= .5*contract('ibac,acij->bj',govvv,dx[1])

    Hdx += .5*contract('ikaj,abik->bj',goovo,dx[1])

    #Doubles/Doubles
    #BHA
    Hdx2 += contract('ab,acjk->bcjk',fvv,dx[1]) 
    Hdx2 -= contract('ac,abjk->bcjk',fvv,dx[1]) 

    Hdx2 -= contract('ij,bcik->bcjk',foo,dx[1]) 
    Hdx2 += contract('ik,bcij->bcjk',foo,dx[1])

    t2_vec = t_vec[:(len(t_vec)-n_occs*n_noccs)]
    t2_vec.shape = (int(n_noccs*(n_noccs-1)/2), int(n_occs*(n_occs-1)/2))
    g2_vec = Populate_VV([gvvvv, f], n_occs, n_noccs, n_orbitals, f, g, hf_energy)
    g2_vec.shape = (int(n_noccs*(n_noccs-1)/2), int(n_noccs*(n_noccs-1)/2))
    con = contract('ab,bc->ac',g2_vec, t2_vec)
    con = con.flatten()
    con = np.concatenate((con,np.zeros((n_occs*n_noccs))))
    con = Populate_Tensor(con, n_occs, n_noccs, n_orbitals, f, g, hf_energy)

    Hdx2 += con[1]
    #Hdx2 += .5*contract('bcad,adjk->bcjk',gvvvv,dx[1])

    Hdx2 += .5*contract('iljk,adil->adjk',goooo,dx[1])
    Hdx2 -= contract('icjd,bdik->cbkj',govov,dx[1])
    Hdx2 += contract('ickd,bdij->cbkj',govov,dx[1])
    Hdx2 += contract('ibjd,cdik->cbkj',govov,dx[1])
    Hdx2 -= contract('ibkd,cdij->cbkj',govov,dx[1])
    return(Populate_Vector([Hdx2, Hdx], n_occs, n_noccs, n_orbitals, f, g, hf_energy))


def Populate_Tensor(t_vec, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
    t1 = np.array(t_vec[(len(t_vec)-n_occs*n_noccs):])
    t1.shape = (n_noccs,n_occs)
    t2 = np.zeros((n_noccs, n_noccs, n_occs, n_occs))
    t2 = np.zeros((n_noccs, n_noccs, n_occs, n_occs))
    t2 = np.zeros((n_noccs, n_noccs, n_occs, n_occs))
    t2 = np.zeros((n_noccs, n_noccs, n_occs, n_occs))
    O = (np.triu_indices(n_occs, 1))
    V = (np.triu_indices(n_noccs, 1))
    v2 = np.array(t_vec[:(len(t_vec)-n_occs*n_noccs)])
    v2.shape = (len(V[0]),len(O[0]))
    t2[V[0][:,None],V[1][:,None],O[0],O[1]] = v2
    t2[V[0][:,None],V[1][:,None],O[1],O[0]] = -v2
    t2[V[1][:,None],V[0][:,None],O[0],O[1]] = -v2
    t2[V[1][:,None],V[0][:,None],O[1],O[0]] = v2
    return(t1,t2)

def Populate_Vector(t_tensor, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
   O = (np.triu_indices(n_occs, 1))
   V = (np.triu_indices(n_noccs, 1))
   t_vector = (t_tensor[0][V][(slice(None),) + O])
   t_vector = t_vector.flatten()
   t1 = t_tensor[1]
   t1 = t1.flatten()
   t_vector = np.concatenate((t_vector,t1))
   return t_vector

def Populate_VV(t_tensor, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
   O = (np.triu_indices(n_noccs, 1))
   V = (np.triu_indices(n_noccs, 1))
   t_vector = (t_tensor[0][V][(slice(None),) + O])
   t_vector = t_vector.flatten()
   t1 = t_tensor[1]
   t1 = t1.flatten()
   return t_vector

def CG_Solver(t_vec, n_occs, n_noccs, n_orbitals, f, g, hf_energy):
     k = 0
     b = -get_gradient(n_occs, n_noccs, n_orbitals, f, g, hf_energy)
     x = t_vec
     Ax0 = Hessian_Action(x, n_occs, n_noccs, n_orbitals, f, g, hf_energy)
     r = Ax0-b
     p = -r
     r_k_norm = np.dot(r,r)
     print('%5s|%20.16s'%(('Iter.', 'Residual Norm')))
     while r_k_norm > 1e-16:
         Ap = (Hessian_Action(p, n_occs, n_noccs, n_orbitals, f, g, hf_energy))
         alpha = r_k_norm/np.dot(p,Ap)
         x += alpha * p
         r += alpha * Ap
         r_kplus1_norm = np.dot(r,r)
         beta = r_kplus1_norm/r_k_norm
         r_k_norm = r_kplus1_norm
         p = beta * p - r
         print('-'*30)
         k += 1
         print('{}'.format(k).ljust(5)+'|'+'%20.16f'%(r_k_norm)) 


     end = timer()
     print('\nConverged in {} iterations  ({} seconds).'.format(k, end-start)+'\n')
     print('HF energy:'.ljust(20)+'%20.16f Eh\n'%(hf_energy))
     print('Converged energy:'.ljust(20)+'%20.16f Eh\n'%(Energy(x, -b, n_occs, n_noccs, n_orbitals, f, g, hf_energy)+hf_energy))

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
    
    hf_energy0, wfn = psi4.energy('scf', return_wfn = True, scf_type = 'pk')
    
    #h = 1-electron integral 2-tensor (<i|h|a>)
    #g = 2-electron integral 4-tensor (<ij|ab>, not <ij||ab>)
    #f = Fock matrix 2-tensor (<i|h|i>+\sum_a(<ii|aa>))
    #All tensors in MO, spinorbital representation
    
    nr = molecule.nuclear_repulsion_energy()
    mints = psi4.core.MintsHelper(wfn.basisset())
    C = wfn.Ca()
    h0 = mints.ao_kinetic()
    f0 = wfn.Fa()
    h0.add(mints.ao_potential())
    h0.transform(C)
    f0.transform(C)
    h0 = np.kron(np.array(h0),(np.array([[1,0],[0,1]])))
    g0 = np.asarray(mints.mo_spin_eri(C, C))
    g0 = g0.swapaxes(1, 2)
    g1 = copy.copy(g0)
    g1 = g1.swapaxes(2, 3)
    g0 -= g1
    f0 = np.array(f0)
    f0 = np.kron(f0, (np.array([[1,0],[0,1]])))
    n_orbitals0 = 2*wfn.nmo()
    n_occs0 = 2*wfn.doccpi()[0]
    n_noccs0 = n_orbitals0-n_occs0
    n_ops = 0

    for i in range(0, n_occs0):
         for j in range(i+1, n_occs0):
             for a in range(0, n_noccs0):
                 for b in range(a+1, n_noccs0):
                     n_ops += 1
                      
    for i in range(0, n_occs0):
        for a in range(0, n_noccs0):
            n_ops += 1

    t_vec0 = np.zeros(n_ops)
    CG_Solver(copy.copy(t_vec0), copy.copy(n_occs0), copy.copy(n_noccs0), copy.copy(n_orbitals0), copy.copy(f0), copy.copy(g0), copy.copy(hf_energy0))


