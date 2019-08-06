import psi4
import numpy as np
import copy
import scipy
print('\n'+'- '*25+'\n')
print('  '*10+'STARLIGHT'+' '*10+'\n')
print(' '*4+'An H.R. Grimsley & N.J. Mayhall Algorithm\n')
print('- '*25+'\n')
psi4.core.set_output_file('output.dat', False)

molecule = psi4.geometry("""
H 0 0 0
H 0 0 1.5
symmetry c1
""")

#Psi4 Calculations
psi4.set_options({'basis': 'sto-3g', 'scf_type': 'pk'})
print('Performing Hartree-Fock Orbital Optimization...\n')
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
#f = np.diag(np.linalg.eig(np.array(wfn.Fa()))[0])
f = np.array(f)
f = np.kron(f, (np.array([[1,0],[0,1]])))
print('HF energy: %20.16f Eh\n'%(hf_energy))
n_orbitals = 2*wfn.nmo()
n_occs = 2*wfn.doccpi()[0]
occs = [i for i in range(0, n_occs)]
noccs = [i for i in range(n_occs, n_orbitals)]
#HF Sanity Check
assert(abs(hf_energy-.5*np.einsum('abab', g[:n_occs,:n_occs,:n_occs,:n_occs])-.5*np.einsum('abba', g[:n_occs,:n_occs,:n_occs,:n_occs])+np.einsum('aa', h[:n_occs,:n_occs])<.000001))

         
#compute gradient
print('Computing gradient ...\n')
gradient = []

#singles
for i in range(0, n_occs):
    for a in range(n_occs, n_orbitals):
        gradient.append(0)

#doubles
for i in range(0, n_occs):
    for j in range(i+1, n_occs):
        for a in range(n_occs, n_orbitals):
            for b in range(a+1, n_orbitals):
                gradient.append(2**.5*(g[i][j][b][a]-g[i][j][a][b]))

gradient = np.array([gradient]).T


#Hessian Function
#Input: dx - trial vector
#Output: Hdx - Hessian action on dx

def Vec_to_Tensor(vector):
    t1 = np.zeros((n_occs, n_orbitals-n_occs)) 
    t2 = np.zeros((n_occs, n_orbitals-n_occs, n_occs, n_orbitals-n_occs)) 
    for i in range(0, len(occs)):
        for a in range(0, len(noccs)):
            if i%2==a%2:
                t1[i][a] = vector[str(i)+','+str(a)]
            for j in range(i, len(occs)):
                for b in range(a, len(noccs)):

                    if (j==i and b==a) or ((a+b)%2 != (j+i)%2):
                        continue                
                    t2[a][b][i][j] = vector[str(i)+','+str(a)+','+str(j)+','+str(b)]
                    t2[b][a][j][i] = vector[str(i)+','+str(a)+','+str(j)+','+str(b)]
                    t2[a][b][j][i] = -vector[str(i)+','+str(a)+','+str(j)+','+str(b)]
                    t2[b][a][i][j] = -vector[str(i)+','+str(a)+','+str(j)+','+str(b)]
    return [t1, t2]

def Tensor_to_Vec(tensor):
    t1 = tensor[0]
    t2 = tensor[1]
    vector = {}    
    for i in range(0, len(occs)):
        for a in range(0, len(noccs)):
            if i%2==a%2:
                vector[str(i)+','+str(a)] = t1[i][a]
            for j in range(i, len(occs)):
                for b in range(a, len(noccs)):

                    if (j==i and b==a) or ((a+b)%2 != (j+i)%2):
                        continue
                    vector[str(i)+','+str(a)+','+str(j)+','+str(b)] = t2[a][b][i][j]
    return vector


def Hessian_Action(dx):
    print('Computing Hessian ...\n')
    Hdx = 0
    Hdx2 = 0

    
    #Singles/Singles
    Hdx += np.einsum('ijab,jb->ia',g[:n_occs,:n_occs,n_occs:,n_occs:],dx[0])
    Hdx -= np.einsum('ijba,jb->ia',g[:n_occs,:n_occs,n_occs:,n_occs:],dx[0])
    Hdx -= np.einsum('ij,ja->ia',f[:n_occs,:n_occs],dx[0]) 
    Hdx += np.einsum('ab,jb->ja',f[n_occs:,n_occs:],dx[0]) 
    Hdx -= np.einsum('ibja,jb->ia',g[:n_occs,n_occs:,:n_occs,n_occs:],dx[0])
    Hdx += np.einsum('ibaj,jb->ia',g[:n_occs,n_occs:,n_occs:,:n_occs],dx[0])
 
    #Singles/Doubles
    Hdx += np.einsum('icjk,jbkc->ib',g[:n_occs,n_occs:,:n_occs,:n_occs],dx[1])
    Hdx -= np.einsum('ickj,jbkc->ib',g[:n_occs,n_occs:,:n_occs,:n_occs],dx[1])     
    Hdx -= np.einsum('bcak,jbkc->ja',g[:n_occs,:n_occs,n_occs:,:n_occs],dx[1])
    Hdx += np.einsum('bcka,jbkc->ja',g[:n_occs,:n_occs,:n_occs,n_occs:],dx[1])     
    
    #Doubles/Singles
    Hdx2 += np.einsum('ikaj,jb->iakb',g[:n_occs,:n_occs,n_occs:,:n_occs],dx[0]) 
    Hdx2 -= np.einsum('ikja,jb->iakb',g[:n_occs,:n_occs,:n_occs,n_occs:],dx[0]) 
    Hdx2 -= np.einsum('ibac,jb->iajc',g[:n_occs,n_occs:,n_occs:,n_occs:],dx[0]) 
    Hdx2 += np.einsum('ibca,jb->iajc',g[:n_occs,n_occs:,n_occs:,n_occs:],dx[0]) 

    #Doubles/Doubles
    gvvvv = g[n_occs:,n_occs:,n_occs:,n_occs:]
    goooo = g[:n_occs,:n_occs,:n_occs,:n_occs]
    goovv = g[:n_occs,:n_occs,n_occs:,n_occs:]
    gvvoo = g[n_occs:,n_occs:,:n_occs,:n_occs]
    gvovo = g[n_occs:,:n_occs,n_occs:,:n_occs]
    gvoov = g[n_occs:,:n_occs,:n_occs,n_occs:]
    govov = g[:n_occs,n_occs:,:n_occs,n_occs:]
    govvo = g[:n_occs,n_occs:,n_occs:,:n_occs]
    foo = f[:n_occs,:n_occs]
    fvv = f[n_occs:,n_occs:]

    Hdx2 += np.einsum('ab,acjk->bcjk',fvv,dx[1]) 
    Hdx2 -= np.einsum('ac,abjk->bcjk',fvv,dx[1]) 

    Hdx2 -= np.einsum('ij,bcik->bcjk',foo,dx[1]) 
    Hdx2 += np.einsum('ik,bcij->bcjk',foo,dx[1])
 
    Hdx2 += .5*np.einsum('bcad,adjk->bcjk',gvvvv,dx[1])
    Hdx2 -= .5*np.einsum('bcda,adjk->bcjk',gvvvv,dx[1])
 
    Hdx2 += .5*np.einsum('iljk,adil->adjk',goooo,dx[1])
    Hdx2 -= .5*np.einsum('ilkj,adil->adjk',goooo,dx[1])
    
    Hdx2 -= np.einsum('icjd,bdik->cbkj',govov,dx[1])
    Hdx2 += np.einsum('icdj,bdik->cbkj',govvo,dx[1])
    Hdx2 += np.einsum('ickd,bdij->cbkj',govov,dx[1])
    Hdx2 -= np.einsum('icdk,bdij->cbkj',govvo,dx[1])
    Hdx2 += np.einsum('ibjd,cdik->cbkj',govov,dx[1])
    Hdx2 -= np.einsum('ibdj,cdik->cbkj',govvo,dx[1])
    Hdx2 -= np.einsum('ibkd,cdij->cbkj',govov,dx[1])
    Hdx2 += np.einsum('ibdk,cdij->cbkj',govvo,dx[1])
    

    return(Hdx, Hdx2)

def Projective_Engine():
    trial_vec = {}
    for i in range(0, len(occs)):
        for a in range(0, len(noccs)):
            if a%2==i%2:
                trial_vec[str(i)+','+str(a)] = 1
            for j in range(i, len(occs)):
                for b in range(a, len(noccs)):
                    if (a+b)%2!=(i+j)%2 or (i==j and b==a):
                        continue
                    trial_vec[str(i)+','+str(a)+','+str(j)+','+str(b)] = 1

    for i in range(0, len(occs)):
        for a in range(0, len(noccs)):
            if a%2==i%2:
                trial_vec[str(i)+','+str(a)] = 1/np.sqrt(len(trial_vec))
            for j in range(i, len(occs)):
                for b in range(a, len(noccs)):
                    if (a+b)%2!=(i+j)%2 or (i==j and b==a):
                        continue
                    trial_vec[str(i)+','+str(a)+','+str(j)+','+str(b)] = 1/np.sqrt(len(trial_vec))
    trial_vec = Vec_to_Tensor(paramdict)

    for i in range(0, 100):
        maxi = max([np.amax(abs(trial_vec[0]))**2,np.amax(abs(trial_vec[1]))**2])
        trial_vec = Hessian_Action([trial_vec[0]/maxi, trial_vec[1]/maxi])
    trial_vec = Tensor_to_Vec(trial_vec)

      
        
      
if __name__ == '__main__':

    paramdict = {'0,0,1,1': 2.2, '0,0': 3.2, '1,1': 2.1}
    dx = Vec_to_Tensor(paramdict)
    print(Tensor_to_Vec(Hessian_Action(dx)))
