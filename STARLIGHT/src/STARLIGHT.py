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
print('Gradient vector:\n\n %s\n'%gradient)

#Hessian Function
#Input: dx - trial vector
#Output: Hdx - Hessian action on dx
def Hessian_Action(dx):
    print('Computing Hessian action on trial vector...\n')

    #singles/singles
    Hdx = .5*np.einsum('ijab,jb->ia',g,dx[0])
    Hdx -= .5*np.einsum('ijba,jb->ia',g,dx[0])
    Hdx += np.einsum('ij->i',f)
    Hdx -= np.einsum('ab->a',f)
    Hdx += np.einsum('ibaj,jb->ia',g,dx[0])
    Hdx -= np.einsum('ibja,jb->ia',g,dx[0])    
    #singles/doubles            
    
    #doubles/singles
    
    #doubles/doubles 
    print("Computed Hessian action: \n\n%s \n"%(Hdx))
if __name__ == '__main__':
    dx = []
    dx.append(np.ones((4,4))) 
    dx.append(np.ones((4,4,4,4)))
    Hessian_Action(dx)
