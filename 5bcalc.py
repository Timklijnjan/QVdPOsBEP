import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import pickle

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

#Loading data
with open("Thetavarsfile0__pi3__pi2","rb") as fp:
    Datlist = pickle.load(fp) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

lentimes = Datlist[0]
lentimes2 = Datlist[1]
dt = Datlist[2]
Dt = Datlist[3]
thetlist = Datlist[4]
nummont = Datlist[5]
thetdata = Datlist[6]


#Retrieving extra time
dtt = 1/dt
extime = math.ceil(Dt/dt)+2

Cphiphaseavgthet=[]
for i in range(len(thetlist)):
    trajdata = thetdata[i]
    Cphiphaseavg=[]
    for j in range(nummont):
        datexp =trajdata[j]
        Cphi = datexp[0][extime:lentimes2]/np.sqrt(datexp[1][extime:lentimes2]*datexp[2][extime:lentimes2])
        Cphiphase = np.angle(Cphi)
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
    Cphiphaseavgthet.append(Cphiphaseavg)


hist1, bin_edges1 =np.histogram(Cphiphaseavgthet[0],range =(-np.pi,np.pi),bins = 100)
hist2, bin_edges2 =np.histogram(Cphiphaseavgthet[1],range =(-np.pi,np.pi),bins = 100)
hist3, bin_edges3 =np.histogram(Cphiphaseavgthet[2],range =(-np.pi,np.pi),bins = 100)
plt.bar(bin_edges1[:-1],hist1/nummont,width = 2*np.pi/100, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/nummont,width = 2*np.pi/100, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/nummont,width = 2*np.pi/100, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],["$-\pi$","$-\pi /2$","$0$","$\pi /2$","$\pi$"])
plt.xlabel("$\Delta \phi_{\psi}$")
plt.ylabel("$P(\Delta \phi_{\psi})$")
plt.show()