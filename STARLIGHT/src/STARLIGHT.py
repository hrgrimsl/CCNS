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
nr = h2.nuclear_repulsion_energy()
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
        for j in range(i, n_occs):
            for b in range(n_occs, n_orbitals):
                if a!=b and i!=j:
                    #first operator (6 terms)
                    hab = (2*g[i][j][a][b]-g[i][j][b][a])
                    bha = -g[b][i][a][j]+2*g[b][i][j][a]
                    hessian[-1].append(hab+bha)
                elif a==b and i!=j:
                    bha = -g[i][a][j][a]-h[i][j]
                    for n in range(0, n_occs):                       
                        bha -= g[i][n][j][n]
                        bha += g[i][n][n][j]
                    hab = 0 
                    hessian[-1].append(0)
                elif a!=b and i==j:
                    bha = -g[b][i][a][i]-h[a][b]
                    for n in range(0, n_occs):
                        bha -= g[a][n][b][n]
                        bha += g[a][n][n][b]
                    hab = g[i][i][a][b]-g[i][i][b][a]
                    hessian[-1].append(0)
                elif a==b and i==j:
                    hab = -hf_energy+g[i][i][a][a]
                       
                    bha = g[i][a][a][i]+nr
                    indices = [q for q in range(0, n_occs)]

                    indices.append(a)

                    for k in indices:
                        fac = 1
                        if k != a and k != i:
                            fac *= 2
                        bha += fac*h[k][k]    
                        for l in indices:
                            xfac = 2
                            if (k == a and l == i) or (k == i and l == a):
                                xfac = 0
                            elif k == a or l == i or k == i or l == a:
                                xfac = 1
                            fac = 1
                            if k != a and k != i:
                                fac *= 2
                            if l != a and l != i:
                                fac *= 2
                            bha += (.5*fac*g[k][l][k][l]-.5*xfac*g[k][l][l][k])
                    print(bha)
                    print(hab)
                    hessian[-1].append(bha+hab)
hessian = hessian
print(hessian)            


