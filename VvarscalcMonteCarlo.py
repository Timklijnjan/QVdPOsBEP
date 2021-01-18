import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import pickle

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

#Loading data
with open("Vvarsfile0_1","rb") as fp:
    Datlist = pickle.load(fp) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

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

CphiphaseavgVvar=[]
CphimodavgVvar=[]
for i in range(len(Vlist)):
    trajdata = Vdata[i]
    Cphiphaseavg=[]
    Cphimodavg=[]
    for j in range(nummont):
        datexp =trajdata[j]
        Cphi = datexp[0][extime:lentimes2]/np.sqrt(datexp[1][extime:lentimes2]*datexp[2][extime:lentimes2])
        Cphiphase = np.angle(Cphi)
        Cphimod = np.abs(Cphi)
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
        Cphimodavg.append(timeavg(Cphimod,dt,dtt))
    CphiphaseavgVvar.append(np.var(np.array(Cphiphaseavg)))
    CphimodavgVvar.append(np.var(np.array(Cphimodavg)))

#2e file
#Loading data
with open("Vvarsfile0_215","rb") as fp:
    Datlist = pickle.load(fp) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

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

for i in range(len(Vlist)):
    trajdata = Vdata[i]
    Cphiphaseavg=[]
    Cphimodavg=[]
    for j in range(nummont):
        datexp =trajdata[j]
        Cphi = datexp[0][extime:lentimes2]/np.sqrt(datexp[1][extime:lentimes2]*datexp[2][extime:lentimes2])
        Cphiphase = np.angle(Cphi)
        Cphimod = np.abs(Cphi)
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
        Cphimodavg.append(timeavg(Cphimod,dt,dtt))
    CphiphaseavgVvar.append(np.var(np.array(Cphiphaseavg)))
    CphimodavgVvar.append(np.var(np.array(Cphimodavg)))

#3e file
#Loading data
with open("Vvarsfile0_464__1__2_15__4_64","rb") as fp:
    Datlist = pickle.load(fp) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

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

for i in range(len(Vlist)):
    trajdata = Vdata[i]
    Cphiphaseavg=[]
    Cphimodavg=[]
    for j in range(nummont):
        datexp =trajdata[j]
        Cphi = datexp[0][extime:lentimes2]/np.sqrt(datexp[1][extime:lentimes2]*datexp[2][extime:lentimes2])
        Cphiphase = np.angle(Cphi)
        Cphimod = np.abs(Cphi)
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
        Cphimodavg.append(timeavg(Cphimod,dt,dtt))
    CphiphaseavgVvar.append(np.var(np.array(Cphiphaseavg)))
    CphimodavgVvar.append(np.var(np.array(Cphimodavg)))

#4e file
#Loading data
with open("Vvarsfile10__21_5__46_4__100","rb") as fp:
    Datlist = pickle.load(fp) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

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

for i in range(len(Vlist)):
    trajdata = Vdata[i]
    Cphiphaseavg=[]
    Cphimodavg=[]
    for j in range(nummont):
        datexp =trajdata[j]
        Cphi = datexp[0][extime:lentimes2]/np.sqrt(datexp[1][extime:lentimes2]*datexp[2][extime:lentimes2])
        Cphiphase = np.angle(Cphi)
        Cphimod = np.abs(Cphi)
        Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
        Cphimodavg.append(timeavg(Cphimod,dt,dtt))
    CphiphaseavgVvar.append(np.var(np.array(Cphiphaseavg)))
    CphimodavgVvar.append(np.var(np.array(Cphimodavg)))

Vvals = [0.1,0.215,0.464,1,2.15,4.64,10,21.5,46.6,100]
fig,ax =plt.subplots()
ax.semilogx(Vvals,CphimodavgVvar,'b--',marker='o')
ax.set_xscale('log')
ax.set_xlabel("V (units of $\gamma$)")
ax.set_ylabel("Var[|$C_{\psi}$|]", color='b')
ax.set_yticks([0.00,0.02,0.04,0.06])
ax2 = ax.twinx()
ax2.semilogx(Vvals,CphiphaseavgVvar,'g--',marker='o')
ax2.set_ylabel("Var[$\Delta\phi$]",color ='g')
ax2.set_yticks([0,1,2,3])
plt.show()