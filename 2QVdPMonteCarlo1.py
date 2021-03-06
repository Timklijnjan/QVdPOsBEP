import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math

def rx1x2(t,Dt,x1,x2,dt):
    #assuming times is spaced with equal increments
    dtt = 1/dt
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

        # second integration for averages used in in delta<x>
        #for si in s:
            #irmin = si - cDt
            #irmax = si + cDt+1
            #rint1.append(scipy.integrate.simps(Dtt*x1[irmin: irmax],dx=dt))
            #rint2.append(scipy.integrate.simps(Dtt*x2[irmin: irmax],dx=dt))

        # first integration for averages of whole exprissions
        sint12 = scipy.integrate.simps((x1[isminti:ismaxti]-rint1a[isminri:ismaxri])*(x2[isminti:ismaxti]-rint2a[isminri:ismaxri]),dx=dt)
        sint11 = scipy.integrate.simps((x1[isminti:ismaxti]-rint1a[isminri:ismaxri])**2,dx=dt)
        sint22 = scipy.integrate.simps((x2[isminti:ismaxti]-rint2a[isminri:ismaxri])**2,dx=dt)
        result.append(sint12/np.sqrt(sint11*sint22))
    return np.array(result)


n = 2
omg1 = 2*np.pi
gam1u = 0.01
gam2u = gam1u
gam1d = 10000
gam2d = gam1d
theta = 0
V = 10*gam1u
omg =[omg1,omg1+0.1*gam1u]
Domg = omg[1]-omg[0]
Dt = 4*np.pi/omg1
N = (3*gam1u+V)*(3*gam1u*(Domg**2+9*gam1u**2)+(Domg**2+27*gam1u**2)*V+8*gam1u*V**2)
#p = gam1u*(3*gam1u+V)*((omg[1]-omg[0])**2+(3*gam1u+V)**2)/N
times = np.linspace(0.0,100/omg[0], 100000)
dt = times[1]-times[0]
extime = math.ceil(Dt/dt)+2
times2 = np.linspace(0.0,100/omg[0]+extime*dt,len(times)+extime)
pi11 = 1-(gam1u*(5*gam1u+2*V)*(Domg**2+(3*gam1u+V)**2))/N
pi22 = gam1u*(2*gam1u+V)*(Domg**2+(3*gam1u+V)**2)/N
pi33 = pi22
pi44 = gam1u**2*(Domg**2+(3*gam1u+V)**2)/N
pi23 =gam1u*V*(gam1u+V)*(3*gam1u+V-complex(0,Domg))*np.exp(complex(0,-theta))/N
pi32 = np.conjugate(pi23)
pi = np.array([[pi11,0,0,0],[0,pi22,pi23,0],[0,pi32,pi33,0],[0,0,0,pi44]])
pi1 = Qobj([[1],[0],[0],[0]])
p1 = pi11
pi2 = 1/np.sqrt(2)*Qobj([[0],[-pi23/np.abs(pi23)],[1],[0]])
p2 = -np.abs(pi23)+pi22
pi3 = 1/np.sqrt(2)*Qobj([[0],[pi23/np.abs(pi23)],[1],[0]])
p3 = np.abs(pi23)+pi22
pi4 = Qobj([[0],[0],[0],[1]])
p4 = pi44
pn = [pi1,pi2,pi3,pi4]
pick = np.random.choice(4,p=np.array([p1,p2,p3,p4]))
psi0 = pn[pick]
a1 = tensor(destroy(n),qeye(n))
a1d = a1.dag()
a2 = tensor(qeye(n),destroy(n))
a2d = a2.dag()
a=[a1,a2]
ad = [a1d,a2d]
H = 0 #use for ssesolve
for i in range(len(omg)):
    H = H + omg[i]*ad[i]*a[i]
L = [np.sqrt(gam1d)*a1**2,np.sqrt(gam1u)*a1d, np.sqrt(gam2d)*a2**2, np.sqrt(gam2u)*a2d, np.sqrt(V)*(a1-np.exp(complex(0,theta))*a2)]
#Heff = H
#for i in range(len(L)):
#    Heff = Heff -complex(0,1)* 1/2*L[i].dag()*L[i]
e_ops = [a1d*a2,a1d*a1,a2d*a2,(a1+a1d)/np.sqrt(2),(a2+a2d)/np.sqrt(2)]
data = ssesolve(H,psi0,times2,sc_ops=L, e_ops = e_ops, method="homodyne")
Cphi = data.expect[0][1:len(times)]/np.sqrt(data.expect[1][1:len(times)]*data.expect[2][1:len(times)])
Cphimod = np.abs(Cphi)
Cphiphase = np.angle(Cphi)
Cpi = V*(gam1u+V)/((3*gam1u+V)*np.sqrt(Domg**2+(3*gam1u+V)**2))*np.ones(len(times))
Dphipi = (theta-np.arctan(Domg/(3*gam1u+V)))*np.ones(len(times))

x1=data.expect[3]
x2=data.expect[4]
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1]
pears = rx1x2(t,Dt,x1,x2,dt)

fig,ax =plt.subplots()
ax.plot(times[1:],Cphimod,'b')
ax.plot(times[1:],Cpi[1:],'b--')
ax.plot(t,pears,'orange')
plt.xticks([0, 20/omg[0], 40/omg[0],60/omg[0],80/omg[0],100/omg[0]],["0","20","40","60","80","100"])
ax.set_xlabel("$\omega_{1}t$")
ax.set_ylabel("|$C_{\psi (t)}$|", color='b')
ax.set_yticks([-1.00,-0.75,-0.50,-0.25,0.00,0.25,0.50,0.75,1.00])
ax.set_yticklabels(['-1.00','','-0.50','','0.00','','0.50','','1.00'])
ax2 = ax.twinx()
ax2.plot(times[1:],Cphiphase,'g')
ax2.plot(times[1:],Dphipi[1:],'g--')
ax2.set_ylabel("$\Delta\phi (t)$",color ='g')
ax2.set_yticks([-1.00*np.pi,-0.50*np.pi,0.00,0.50*np.pi,1.00*np.pi])
ax2.set_yticklabels(['$-\pi$','$-\pi /2$','$0$','$\pi /2$','$\pi$'])
plt.show()