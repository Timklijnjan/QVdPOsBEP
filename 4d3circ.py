import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import time

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

#Loading data
#[len(times),len(times2),dt,Dt,Vlist,Domgint,nummont]
Datlist = qload("tempory3oscircV20D20_state_dumpp_1000param")

lentimes = Datlist[0]
lentimes2 = Datlist[1]
dt = Datlist[2]
Dt = Datlist[3]
Vlist = Datlist[4]
Domgint=Datlist[5]
nummont = Datlist[6]


#Retrieving extra time
dtt = 1/dt
extime = math.ceil(Dt/dt)+2


"""fileparamlist=['V5D1','V100D1','V20D20']
for fileparam in fileparamlist:
    Savg1 =[]
    Savg2 =[]
    Savg3=[]
    for l in [100,200,300,400,500,600,700,800,900,1000]:
        filename = './tempory3oscirc'+fileparam+'_state_dumpp_' +str(l)
        trajdata = qload(filename)
        for j in range(100):
            start_time=time.time()
            datexp =trajdata[j]
            S1=[]
            S2=[]
            S3=[]
            dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
            for k in range(extime,lentimes2):
                state_obj =dataexp_unroll[k]
                elt1 =state_obj[0][0][0]
                elt2 =state_obj[1][0][0]
                elt3 =state_obj[2][0][0]
                elt4 =state_obj[3][0][0]
                elt5 =state_obj[4][0][0]
                elt6 =state_obj[5][0][0]
                elt7 =state_obj[6][0][0]
                elt8 =state_obj[7][0][0]
                celt1 = np.conjugate(elt1)
                celt2 = np.conjugate(elt2)
                celt3 = np.conjugate(elt3)
                celt4 = np.conjugate(elt4)
                celt5 = np.conjugate(elt5)
                celt6 = np.conjugate(elt6)
                celt7 = np.conjugate(elt7)
                celt8 = np.conjugate(elt8)
                S1.append(entropy_vn(Qobj([[elt1*celt1+elt2*celt2+elt3*celt3+elt4*celt4,elt1*celt5+elt2*celt6+elt3*celt7+elt4*celt8],[elt5*celt1+elt6*celt2+elt7*celt3+elt8*celt4,elt5*celt5+elt6*celt6+elt7*celt7+elt8*celt8]],dims=[[2],[2]])))
                S2.append(entropy_vn(Qobj([[elt1*celt1+elt2*celt2+elt5*celt5+elt6*celt6,elt1*celt3+elt5*celt7+elt2*celt4+elt6*celt8],[elt3*celt1+elt7*celt5+elt4*celt2+elt8*celt6,elt3*celt3+elt4*celt4+elt7*celt7+elt8*celt8]],dims=[[2],[2]])))
                S3.append(entropy_vn(Qobj([[elt1*celt1+elt3*celt3+elt5*celt5+elt7*celt7,elt1*celt2+elt3*celt4+elt5*celt6+elt7*celt8],[elt2*celt1+elt6*celt5+elt4*celt3+elt8*celt7,elt2*celt2+elt4*celt4+elt6*celt6+elt8*celt8]],dims=[[2],[2]])))
            Savg1.append(timeavg(S1,dt,dtt))
            Savg2.append(timeavg(S2,dt,dtt))
            Savg3.append(timeavg(S3,dt,dtt))
            print("--- %2.2f seconds ---" % (time.time() - start_time))
    qsave(Savg1,'3oscircSavg12'+fileparam)
    qsave(Savg2,'3oscircSavg13'+fileparam)
    qsave(Savg3,'3oscircSavg23'+fileparam)"""

SavgV=[qload('3oscircSavg23V5D1'),qload('3oscircSavg23V100D1'),qload('3oscircSavg23V20D20')]

hist1, bin_edges1 =np.histogram(SavgV[0],range =(0,0.7),bins = 40)
hist2, bin_edges2 =np.histogram(SavgV[1],range =(0,0.7),bins = 40)
hist3, bin_edges3 =np.histogram(SavgV[2],range =(0,0.7),bins = 40)
plt.bar(bin_edges1[:-1],hist1/1000,width = 0.7/40, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = 0.7/40, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/1000,width = 0.7/40, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([0.0,0.2,0.4,0.6],["0.0","0.2","0.4","0.6"])
plt.xlabel("$S_{\psi}$")
plt.ylabel("$P(S_{\psi})$")
plt.show()