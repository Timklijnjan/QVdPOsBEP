import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math

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
times = np.linspace(0.0,dt*(lentimes-1),lentimes)
times2 = np.linspace(0.0,times[len(times)-1]+extime*dt,len(times)+extime) #array used to obtain a large enough interval for the pearson indicator integration on all of times
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1] #times for which pearson indicator is calculated
"""

a=b*vg(w)p

"""
rV=[]
for i in range(len(Vlist)):
    trajdata = Vdata[i]
    ravg =[]
    for j in range(nummont):
        datexp =trajdata[j]
        r=[]
        x1=[]
        x2=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(0,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            x1.append(dsqrt2*(celt1*elt3+celt2*elt4+celt3*elt1+celt4*elt2))
            x2.append(dsqrt2*(celt1*elt2+celt2*elt1+celt3*elt4+celt4*elt3))
        x1=np.array(x1)
        x2=np.array(x2)
        r=rx1x2(t,Dt,x1,x2,dt,dtt)
        ravg.append(timeavg(r,dt,dtt))
    #CphiphaseavgV.append(Cphiphaseavg)
    #SavgV.append(Savg)     activate when multiple Vs are used

for l in [200,300,400,500,600,700,800,900]:
    filename = './tempory_state_dumpp_'+str(l)
    trajdata = qload(filename)
    for j in range(100):
        datexp =trajdata[j]
        r=[]
        x1=[]
        x2=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(0,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            x1.append(dsqrt2*(celt1*elt3+celt2*elt4+celt3*elt1+celt4*elt2))
            x2.append(dsqrt2*(celt1*elt2+celt2*elt1+celt3*elt4+celt4*elt3))
        x1=np.array(x1)
        x2=np.array(x2)
        r=rx1x2(t,Dt,x1,x2,dt,dtt)
        ravg.append(timeavg(r,dt,dtt))
rV.append(ravg)

ravg=[]
for l in [100,200,300,400,500,600,700,800,900,1000]:
    filename = './temporyV100_state_dumpp_'+str(l)
    trajdata = qload(filename)
    for j in range(100):
        datexp =trajdata[j]
        r=[]
        x1=[]
        x2=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(0,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            x1.append(dsqrt2*(celt1*elt3+celt2*elt4+celt3*elt1+celt4*elt2))
            x2.append(dsqrt2*(celt1*elt2+celt2*elt1+celt3*elt4+celt4*elt3))
        x1=np.array(x1)
        x2=np.array(x2)
        r=rx1x2(t,Dt,x1,x2,dt,dtt)
        ravg.append(timeavg(r,dt,dtt))
rV.append(ravg)

ravg=[]
for l in [100,200,300,400,500,600,700,800,900,1000]:
    filename = './temporyV20_state_dumpp_'+str(l)
    trajdata = qload(filename)
    for j in range(100):
        datexp =trajdata[j]
        r=[]
        x1=[]
        x2=[]
        dataexp_unroll = datexp[0] #for single trajectories extra [0] is needed
        for k in range(0,lentimes2):
            state_obj =dataexp_unroll[k]
            elt1 =state_obj[0][0][0]
            elt2 =state_obj[1][0][0]
            elt3 =state_obj[2][0][0]
            elt4 =state_obj[3][0][0]
            celt1 = np.conjugate(elt1)
            celt2 = np.conjugate(elt2)
            celt3 = np.conjugate(elt3)
            celt4 = np.conjugate(elt4)
            x1.append(dsqrt2*(celt1*elt3+celt2*elt4+celt3*elt1+celt4*elt2))
            x2.append(dsqrt2*(celt1*elt2+celt2*elt1+celt3*elt4+celt4*elt3))
        x1=np.array(x1)
        x2=np.array(x2)
        r=rx1x2(t,Dt,x1,x2,dt,dtt)
        ravg.append(timeavg(r,dt,dtt))
rV.append(ravg)

qsave(rV,"SavgV5V100V20_4c")

hist1, bin_edges1 =np.histogram(rV[0],range =(-1.0,1.0),bins = 40)
hist2, bin_edges2 =np.histogram(rV[1],range =(-1.0,1.0),bins = 40)
hist3, bin_edges3 =np.histogram(rV[2],range =(-1.0,1.0),bins = 40)
plt.bar(bin_edges1[:-1],hist1/1000,width = 1/20, color ='r', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges2[:-1],hist2/1000,width = 1/20, color ='b', edgecolor='k',alpha = 0.5)
plt.bar(bin_edges3[:-1],hist3/1000,width = 1/20, color ='g', edgecolor='k',alpha = 0.5)
plt.xlim(min(bin_edges1),max(bin_edges1))
plt.xticks([-1.0,0.0,1.0],["-1.0","","1.0"])
plt.xlabel("$r_{x_{1},x_{2}}$")
plt.ylabel("$P(r_{x_{1},x_{2}})$")
plt.show()

