import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import pickle

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

'''#Loading data
Datlist = qload("temporyV20bD30_state_dumpp_1000param")
lentimes = Datlist[0]
lentimes2 = Datlist[1]
dt = Datlist[2]
Dt = Datlist[3]
Vlist = Datlist[4]
nummont = Datlist[5]


#Retrieving extra time
dtt = 1/dt
extime = math.ceil(Dt/dt)+2

Domgvarsphi=[]
DomgvarsS=[]
for Domgint in ["30"]:   #"min30","min23","min17","min10","min3","0","3","10","17","23","30"
    Cphiphaseavg=[]
    Savg = []
    for l in [100,200,300,400,500,600,700,800,900,1000]:
        filename = './temporyV20bD'+Domgint+'_state_dumpp_'+str(l)
        trajdata = qload(filename)
        for j in range(100):
            datexp =trajdata[j]
            Cphiphase=[]
            S =[]
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
                Cphi= celt3*elt2/np.sqrt(mcmd*(bcb+dcd))
                Cphiphase.append(np.angle(Cphi))
                S.append(entropy_vn(Qobj([[aca+bcb,elt1*celt3+elt2*celt4],[elt3*celt1+elt4*celt2,mcmd]],dims=[[2],[2]])))
            Cphiphaseavg.append(timeavg(Cphiphase,dt,dtt))
            Savg.append(timeavg(S,dt,dtt))
    Domgvarsphi.append(np.var(np.array(Cphiphaseavg)))
    DomgvarsS.append(np.var(np.array(Savg)))

qsave(Domgvarsphi,"DomgvarsphiV20bD30")
qsave(DomgvarsS,"DomgvarsSV20bD30")

print('this is varsphi')
print(str(Domgvarsphi[0]))
print('this is varsS')
print(str(DomgvarsS[0]))'''

finDomgvarsphi= qload("DomgvarsphiV20Dmin30_min23_min17_min10_min3_0_3")
Domgvarsphi2 =qload("DomgvarsphiV20D10_17_23_30")
finDomgvarsphi.extend(Domgvarsphi2)
finDomgvarsS = qload("DomgvarsSV20Dmin30_min23_min17_min10_min3_0_3")
DomgvarsS2 = qload("DomgvarsSV20D10_17_23_30")
finDomgvarsS.extend(DomgvarsS2)
correctDomgvarsphi = qload("DomgvarsphiV20bD30")
finDomgvarsphi[len(finDomgvarsphi)-1]=correctDomgvarsphi[0]
correctvarsS = qload("DomgvarsSV20bD30")
finDomgvarsS[len(finDomgvarsS)-1]=correctvarsS[0]



Domgvals =[-30,-23,-17,-10,-3,0,3,10,17,23,30]
fig,ax =plt.subplots()
ax.plot(Domgvals,finDomgvarsphi,'b--',marker='o')
ax.set_xlabel("$\Delta\omega$ (units of $\gamma$)")
ax.set_ylabel("Var[$\Delta\phi_{\psi}$]", color='b')
plt.xticks([-30,-20,-10,0,10,20,30],["-30","-20","-10","0","10","20","30"])
ax2 = ax.twinx()
ax2.plot(Domgvals,finDomgvarsS,'g--',marker='o')
ax2.set_ylabel("Var[$S_{\psi}$]",color ='g')
plt.show()