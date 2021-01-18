import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import time

def rx1x2(t,Dt,x1,x2,dt,dtt):
    #assuming times is spaced with equal increments
    Dtt =1/Dt
    result = []
    cDt = math.ceil((Dt*0.5)*dtt)
    ismin = math.floor((t[0]-Dt*0.5)*dtt)
    ismax = math.ceil((t[len(t)-1]+Dt*0.5)*dtt)+1
    rint1= []
    rint2= []
    for si in range(ismin,ismax):
        irmin = si - cDt
        irmax = si + cDt+1
        rint1.append(scipy.integrate.simps(Dtt*x1[irmin: irmax],dx=dt))
        rint2.append(scipy.integrate.simps(Dtt*x2[irmin: irmax],dx=dt))
    
    rint1a =np.array(rint1)
    rint2a =np.array(rint2)
    
    for ti in t:
        isminti =math.floor((ti-Dt*0.5)*dtt) #assuming x1 and x2 start at time times[0] and are spaced with dt
        ismaxti = math.ceil((ti+Dt*0.5)*dtt)+1
        
        isminri = isminti-ismin
        ismaxri = ismaxti - ismin

        # first integration for averages of whole exprissions
        sint12 = scipy.integrate.simps((x1[isminti:ismaxti]-rint1a[isminri:ismaxri])*(x2[isminti:ismaxti]-rint2a[isminri:ismaxri]),dx=dt)
        sint11 = scipy.integrate.simps((x1[isminti:ismaxti]-rint1a[isminri:ismaxri])**2,dx=dt)
        sint22 = scipy.integrate.simps((x2[isminti:ismaxti]-rint2a[isminri:ismaxri])**2,dx=dt)
        result.append(sint12/np.sqrt(sint11*sint22))
    return np.array(result)

def timeavg(array, dt, dtt):
    return dtt*scipy.integrate.simps(array,dx=dt)/(len(array)-1)

dsqrt2=1/np.sqrt(2)
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
times = np.linspace(0.0,dt*(lentimes-1),lentimes)
times2 = np.linspace(0.0,times[len(times)-1]+extime*dt,len(times)+extime) #array used to obtain a large enough interval for the pearson indicator integration on all of times
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1] #times for which pearson indicator is calculated


"""fileparamlist=['V5D1','V100D1','V20D20']
for fileparam in fileparamlist:
    ravg12=[]
    ravg13=[]
    ravg23=[]
    for l in [100,200,300,400,500,600,700,800,900,1000]:
        filename = './tempory3os'+fileparam+'_state_dumpp_' +str(l)
        trajdata = qload(filename)
        for j in range(100):
            start_time = time.time()
            datexp =trajdata[j]
            r12=[]
            r13=[]
            r23=[]
            x1=[]
            x2=[]
            x3=[]
            dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
            for k in range(0,lentimes2):
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
                x1.append(dsqrt2*(elt1*celt5+elt2*celt6+elt3*celt7+elt4*celt8+elt5*celt1+elt6*celt2+elt7*celt3+elt8*celt4))
                x2.append(dsqrt2*(elt1*celt3+elt2*celt4+elt3*celt1+elt4*celt2+elt5*celt7+elt6*celt8+elt7*celt5+elt8*celt6))
                x3.append(dsqrt2*(elt1*celt2+elt2*celt1+elt3*celt4+elt4*celt3+elt5*celt6+elt6*celt5+elt7*celt8+elt8*celt7))
            x1=np.array(x1)
            x2=np.array(x2)
            x3=np.array(x3)
            r12=rx1x2(t,Dt,x1,x2,dt,dtt)
            ravg12.append(timeavg(r12,dt,dtt))
            r13=rx1x2(t,Dt,x1,x3,dt,dtt)
            ravg13.append(timeavg(r13,dt,dtt))
            r23=rx1x2(t,Dt,x2,x3,dt,dtt)
            ravg23.append(timeavg(r23,dt,dtt))
            print("--- %2.2f seconds ---" % (time.time() - start_time))
    qsave(ravg12,'3osravg12'+fileparam)
    qsave(ravg13,'3osravg13'+fileparam)
    qsave(ravg23,'3osravg23'+fileparam)"""

rV = [qload('3osravg23V5D1'),qload('3osravg23V100D1'),qload('3osravg23V20D20')]

hist1, bin_edges1 =np.histogram(rV[0],range =(-1.0,1.0),bins = 40)
hist2, bin_edges2 =np.histogram(rV[1],range =(-1.0,1.0),bins = 40)
hist3, bin_edges3 =np.histogram(rV[2],range =(-1.0,1.0),bins = 40)
plt.bar(bin_edges1[:-1],hist1/1000,width = 1/20, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = 1/20, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/1000,width = 1/20, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([-1.0,0.0,1.0],["-1.0","","1.0"])
plt.xlabel("$r_{x_{2},x_{3}}$")
plt.ylabel("$P(r_{x_{2},x_{3}})$")
plt.show()

