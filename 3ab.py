import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

#Loading data
#[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5,(states)],trja1:[...],...],V1:[[]]...]]
Datlist = qload("Vvars5")

lentimes = Datlist[0]
lentimes2 = Datlist[1]
dt = Datlist[2]
Dt = Datlist[3]
Vlist = Datlist[4]
nummont = Datlist[5]
Vdata = Datlist[6]


#Retrieving extra time
dtt = 1/dt
extime = math.ceil(Dt/dt)+2

#defining basic quantum operators
n=2
a1 = tensor(destroy(n),qeye(n))
a1d = a1.dag()
a2 = tensor(qeye(n),destroy(n))
a2d = a2.dag()
e_ops = [a1d*a2,a1d*a1,a2d*a2] # Operators for which expectation values are determined


CphiphaseavgV=[]
SavgV=[]
for i in range(len(Vlist)):
    trajdata = Vdata[i]
    Cphiphaseavg=[]
    Savg =[]
    for j in range(nummont):
        datexp =trajdata[j]
        S=[]
        Cphiphase=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(extime,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            aca = elt1*celt1
            bcb = elt2*celt2
            ccc=elt3*celt3
            dcd=elt4*celt4
            mcmd=ccc+dcd
            S.append(entropy_vn(Qobj([[aca+bcb,elt1*celt3+elt2*celt4],[elt3*celt1+elt4*celt2,mcmd]],dims=[[2],[2]])))
            Cphiphase.append(np.angle(celt3*elt2/np.sqrt(mcmd*(bcb+dcd))))
        Savg.append(timeavg(S,dt,dtt))
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
    #CphiphaseavgV.append(Cphiphaseavg)
    #SavgV.append(Savg)     activate when multiple Vs are used

for l in [200,300,400,500,600,700,800,900]:
    filename = './tempory_state_dumpp_'+str(l)
    trajdata = qload(filename)
    for j in range(100):
        datexp =trajdata[j]
        S=[]
        Cphiphase=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(extime,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            aca = elt1*celt1
            bcb = elt2*celt2
            ccc=elt3*celt3
            dcd=elt4*celt4
            mcmd=ccc+dcd
            S.append(entropy_vn(Qobj([[aca+bcb,elt1*celt3+elt2*celt4],[elt3*celt1+elt4*celt2,mcmd]],dims=[[2],[2]])))
            Cphiphase.append(np.angle(celt3*elt2/np.sqrt(mcmd*(bcb+dcd))))
        Savg.append(timeavg(S,dt,dtt))
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
CphiphaseavgV.append(Cphiphaseavg)
SavgV.append(Savg)

Savg=[]
Cphiphaseavg=[]
for l in [100,200,300,400,500,600,700,800,900,1000]:
    filename = './temporyV50_state_dumpp_'+str(l)
    trajdata = qload(filename)
    for j in range(100):
        datexp =trajdata[j]
        S=[]
        Cphiphase=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(extime,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            aca = elt1*celt1
            bcb = elt2*celt2
            ccc=elt3*celt3
            dcd=elt4*celt4
            mcmd=ccc+dcd
            S.append(entropy_vn(Qobj([[aca+bcb,elt1*celt3+elt2*celt4],[elt3*celt1+elt4*celt2,mcmd]],dims=[[2],[2]])))
            Cphiphase.append(np.angle(celt3*elt2/np.sqrt(mcmd*(bcb+dcd))))
        Savg.append(timeavg(S,dt,dtt))
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
CphiphaseavgV.append(Cphiphaseavg)
SavgV.append(Savg)


hist1, bin_edges1 =np.histogram(SavgV[0],range =(0,0.7),bins = 35)
hist2, bin_edges2 =np.histogram(SavgV[1],range =(0,0.7),bins = 35)
plt.bar(bin_edges1[:-1],hist1/1000,width = 0.70/35, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = 0.70/35, color ='r', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xlabel("$S_{\psi}$")
plt.ylabel("$P(S_{\psi})$")
plt.show()

plt.scatter(SavgV[0],CphiphaseavgV[0],c ='b',alpha=0.5)
plt.scatter(SavgV[1],CphiphaseavgV[1],c ='r',alpha=0.5)
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7],["0.0","","0.2","","0.4","","0.6"])
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],["","$-\pi /2$","0","$\pi /2$",""])
plt.xlabel("$S_{\psi}$")
plt.ylabel("$\Delta\phi_{\psi}$")
plt.show()
