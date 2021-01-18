import numpy as np
import matplotlib.pyplot as plt

posdelt = np.linspace(0, 25, 2501)
negdelt = np.linspace(-25,0,2501)
partposdelt =np.linspace(0,14.203,2501)
partnegdelt=np.linspace(-14.203,0,2501)

fig, ax = plt.subplots()
ax.plot(posdelt, 2.2725*posdelt, 'r');
ax.plot(negdelt, -2.2725*negdelt,'r');
ax.plot(posdelt, 2*posdelt,'b');
ax.plot(negdelt,-2*negdelt,'b');
ax.plot(partposdelt, 4*partposdelt,'g');
ax.plot(partnegdelt, -4*partnegdelt,'g');
ax.set_xlabel("$\Delta_{i,j}$");
ax.set_ylabel('V');
ax.set_title('Arnold tongue for Synchronisation');
ax.fill_between(posdelt,2*posdelt,56.8125,facecolor='cyan',alpha =1)
ax.fill_between(negdelt,-2*negdelt,56.8125,facecolor='cyan',alpha =1)
ax.fill_between(negdelt, -2.2725*negdelt, 56.8125, facecolor='r',alpha=0.65)
ax.fill_between(posdelt, 2.2725*posdelt, 56.8125, facecolor='r',alpha=0.65)
ax.fill_between(partposdelt,4*partposdelt,56.8125,facecolor='g',alpha =0.6)
ax.fill_between(partnegdelt,-4*partnegdelt,56.8125,facecolor='g',alpha =0.6)
plt.show()
