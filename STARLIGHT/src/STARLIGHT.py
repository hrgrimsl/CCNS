import psi4
import numpy as np
import copy
import scipy

psi4.core.set_output_file('output.dat', False)

h2 = psi4.geometry("""
0 1
H 0 0 0
H 0 0 1.5
H 0 0 3
H 0 0 4.5
symmetry c1
""")


psi4.set_options({'basis': 'sto-3g', 'scf_type': 'pk'})

hf_energy, wfn = psi4.energy('scf', return_wfn = True, scf_type = 'pk')
#wfn = psi4.core.Wavefunction.build(h2, psi4.core.get_global_option('basis'))
print('HF energy: %20.16f Eh'%hf_energy)

mints = psi4.core.MintsHelper(wfn.basisset())
C = wfn.Ca()
h = mints.ao_kinetic()
h.add(mints.ao_potential())
h.transform(wfn.Ca())
h = np.array(h)
g = np.asarray(mints.mo_eri(C, C, C, C))
g = g.swapaxes(1, 2)



n_orbitals = wfn.nmo()
n_occs = wfn.doccpi()[0]

#compute gradient
gradient = []

#singles
for i in range(0, n_occs):
    for a in range(n_occs, n_orbitals):
        gradient.append(0)

#doubles
for i in range(0, n_occs):
    for j in range(i, n_occs):
        for a in range(n_occs, n_orbitals):
            for b in range(a, n_orbitals):
                if a==b and i==j:
                    gradient.append((-(2**.5)*(g[i,j,a,b])))
                elif a!=b and i!=j:
                    gradient.append((-(6**.5)*(g[i,j,a,b]-g[i,j,b,a])))
                    gradient.append((-(2**.5)*(g[i,j,a,b]+g[i,j,b,a])))
                elif i==j and a!=b:
                    gradient.append((-2*g[i,j,a,b]))
                elif i!=j and a==b:
                    gradient.append((-(g[i,j,a,b]+-g[i,j,b,a])))
                else:
                    gradient.append(0)
gradient = np.array(gradient)


hf_1e = 0
for ind in range(0,n_occs):
    hf_1e+=h[ind,ind]

#Compute Hessian
hessian = []
occs = [i for i in range(0, n_occs)]


#single
for i in range(0, n_occs):
    for a in range(n_occs, n_orbitals):
        hessian.append([])
        #single/single
        for j in range(0, n_occs):
            for b in range(n_occs, n_orbitals):
                if a!=b and i!=j:
                    #first operator (6 terms)

                    hab = (2*g[i][j][a][b]-g[i][j][b][a])
                    bha = -g[b][i][a][j]+2*g[b][i][j][a]
                    print(str([i,j,a,b]))

                    print(hab)
                    print(bha)
                    hessian[-1].append(hab+bha)
                    

hessian = hessian
print(hessian)            


