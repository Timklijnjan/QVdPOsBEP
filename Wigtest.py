from qutip import * 
import numpy as np
import matplotlib.pyplot as plt

psi0 = fock_dm(9,0)

xvec = np.linspace(-5,5,500)
W_fock = wigner(psi0,xvec,xvec)
fig, axes = plt.subplots(1,1,figsize=(12,12))
cont0 = axes.contourf(xvec,xvec,W_fock,100)
lbl0 = axes.set_title("Wigner function")
axes.set_xlabel('x?')
axes.set_ylabel('p?')
cb = fig.colorbar(cont0)
plt.show()