import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
#import time

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

#Loading data
#[len(times),len(times2),dt,Dt,Vlist,Domgint,nummont]
Datlist = qload("tempory3osV20D20_state_dumpp_1000param")

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
    Cphimodavg12=[]
    Cphimodavg13=[]
    Cphimodavg23=[]
    Cphiphaseavg12=[]
    Cphiphaseavg13=[]
    Cphiphaseavg23=[]

    for l in [100,200,300,400,500,600,700,800,900,1000]:
        filename = './tempory3os'+fileparam+'_state_dumpp_' +str(l)
        trajdata = qload(filename)

        for j in range(100):
            #start_time = time.time()
            datexp =trajdata[j]
            Cphimod12=[]
            Cphimod13=[]
            Cphimod23=[]
            Cphiphase12=[]
            Cphiphase13=[]
            Cphiphase23=[]
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
                aca = elt1*celt1
                bcb = elt2*celt2
                ccc=elt3*celt3
                dcd=elt4*celt4
                ece =elt5*celt5
                fcf= elt6*celt6
                gcg = elt7*celt7
                hch= elt8*celt8
                s5678 = ece+fcf+gcg+hch
                s3478 = ccc+dcd+gcg+hch
                s2468 = bcb+dcd+fcf+hch
                Cphi12= (elt3*celt5+elt4*celt6)/np.sqrt(s5678*s3478)
                Cphi13= (elt2*celt5+elt4*celt7)/np.sqrt(s5678*s2468)
                Cphi23= (elt2*celt3+elt6*celt7)/np.sqrt(s3478*s2468)
                Cphiphase12.append(np.angle(Cphi12))
                Cphimod12.append(np.abs(Cphi12))
                Cphiphase13.append(np.angle(Cphi13))
                Cphimod13.append(np.abs(Cphi13))
                Cphiphase23.append(np.angle(Cphi23))
                Cphimod23.append(np.abs(Cphi23))

            Cphiphaseavg12.append(timeavg(Cphiphase12,dt,dtt))
            Cphiphaseavg13.append(timeavg(Cphiphase13,dt,dtt))
            Cphiphaseavg23.append(timeavg(Cphiphase23,dt,dtt))
            Cphimodavg12.append(timeavg(Cphimod12,dt,dtt))
            Cphimodavg13.append(timeavg(Cphimod13,dt,dtt))
            Cphimodavg23.append(timeavg(Cphimod23,dt,dtt))
            #print("--- %s seconds ---" % (time.time() - start_time))

    qsave(Cphiphaseavg12,'3osCphiphaseavg12'+fileparam)
    qsave(Cphiphaseavg13,'3osCphiphaseavg13'+fileparam)
    qsave(Cphiphaseavg23,'3osCphiphaseavg23'+fileparam)
    qsave(Cphimodavg12,'3osCphimodavg12'+fileparam)
    qsave(Cphimodavg13,'3osCphimodavg13'+fileparam)
    qsave(Cphimodavg23,'3osCphimodavg23'+fileparam)"""


CphiphaseavgV = [qload("3osCphiphaseavg23V100D1"),qload("3osCphiphaseavg23V5D1"),qload("3osCphiphaseavg23V20D20")]
CphimodavgV = [qload("3osCphimodavg23V100D1"),qload("3osCphimodavg23V5D1"),qload("3osCphimodavg23V20D20")]

hist1, bin_edges1 =np.histogram(CphiphaseavgV[0],range =(-2*np.pi/3,2*np.pi/3),bins = 40)
hist2, bin_edges2 =np.histogram(CphiphaseavgV[1],range =(-2*np.pi/3,2*np.pi/3),bins = 40)
hist3, bin_edges3 =np.histogram(CphiphaseavgV[2],range =(-2*np.pi/3,2*np.pi/3),bins = 40)
plt.bar(bin_edges1[:-1],hist1/1000,width = np.pi/30, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = np.pi/30, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/1000,width = np.pi/30, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([-np.pi/2,0,np.pi/2],["$-\pi/2$","","$\pi/2$"])
plt.xlabel("$\Delta\phi_{\psi}$")
plt.ylabel("$P(\Delta\phi_{\psi})$")
plt.show()

hist1, bin_edges1 =np.histogram(CphimodavgV[0],range =(0,1),bins = 40)
hist2, bin_edges2 =np.histogram(CphimodavgV[1],range =(0,1),bins = 40)
hist3, bin_edges3 =np.histogram(CphimodavgV[2],range =(0,1),bins = 40)
plt.bar(bin_edges1[:-1],hist1/1000,width = 1/40, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = 1/40, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/1000,width = 1/40, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([0,0.5,1],["0,0","0.5","1.0"])
plt.xlabel("$|C_{\psi}|$")
plt.ylabel("$P(|C_{\psi}|)$")
plt.show()