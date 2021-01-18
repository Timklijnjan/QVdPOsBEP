import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def modCss(Domg,V):
    return (V+V**2)/((3+V)*np.sqrt(Domg**2+(3+V)**2))

posdelt = np.linspace(0, 25, 2501)
negdelt = np.linspace(-25,0,2501)
x=np.linspace(-25,25,1000)
y=np.linspace(0,50,1000)
X,Y=np.meshgrid(x,y)
modCssres = modCss(X,Y)
fig, ax = plt.subplots() 
cont0 = ax.contourf(X,Y,modCssres,100)
ax.plot(posdelt, 2*posdelt, 'k--');
ax.plot(negdelt, -2*negdelt,'k--');
ax.set_aspect(1)
ax.set_title("Quantum Arnold tongue");
ax.set_xlabel("$\Delta \omega$"+" (orders of "+"$\gamma$"+")");
ax.set_ylabel("$V$"+" (orders of "+"$\gamma$"+")")
cb =fig.colorbar(cont0)
plt.show()