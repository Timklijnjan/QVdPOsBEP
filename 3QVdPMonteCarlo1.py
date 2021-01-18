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
gamup = 0.01
gamdown = 10000
theta = 0
V = 10*gamup
Domg = 0.1*gamup
omg = [omg1,omg1+Domg,omg1+2*Domg]
Dt = 4*np.pi/omg1
times = np.linspace(0.0,100/omg[0], 100000)
dt = times[1]-times[0]
extime = math.ceil(Dt/dt)+2
times2 = np.linspace(0.0,100/omg[0]+extime*dt,len(times)+extime)
#defining basic quantum operators
a1 = tensor(tensor(destroy(n),qeye(n)),qeye(n))
a2 = tensor(tensor(qeye(n),destroy(n)),qeye(n))
a3 = tensor(tensor(qeye(n),qeye(n)),destroy(n))
aa1 = a1*a1
aa2 = a2*a2
aa3 = a3*a3
c_op0 =np.sqrt(gamup)*a1.dag()
c_op1 = np.sqrt(gamdown)*aa1
c_op2 =np.sqrt(gamup)*a2.dag()
c_op3 =np.sqrt(gamdown)*aa2
c_op4 =np.sqrt(gamup)*a3.dag()
c_op5 =np.sqrt(gamdown)*aa3
a1ma2 = a1-a2
a2ma3 = a2-a3
a1ma3 = a1-a3 
omg1a1da1 =  omg1 * a1.dag()*a1
a2da2 = a2.dag()*a2
a3da3 = a3.dag()*a3

sqrtV = np.sqrt(V)
#Creation of the Lindblad operators
c_op=[c_op0, c_op1, c_op2, c_op3, c_op4, c_op5, sqrtV*a1ma2, sqrtV*a2ma3, sqrtV*a1ma3]

#Creation of the system Hamiltonian
Htot = omg1a1da1 + omg[1]*a2da2 + omg[2]*a3da3
#Heff = H
#for i in range(len(L)):
#    Heff = Heff -complex(0,1)* 1/2*L[i].dag()*L[i]
e_ops = [a1.dag()*a2, a1.dag()*a3, a2.dag()*a3, a1.dag()*a1, a2.dag()*a2, a3.dag()*a3, (a1+a1.dag())/np.sqrt(2), (a2+a2.dag())/np.sqrt(2),(a3+a3.dag())/np.sqrt(2)]
final_state = steadystate(Htot,c_op)
qsave(final_state,'3QVdPallssindi')
ssevals = [expect(a1.dag()*a2,final_state), expect(a1.dag()*a3,final_state), expect(a2.dag()*a3,final_state), expect(a1.dag()*a1,final_state), expect(a2.dag()*a2,final_state),expect(a3.dag()*a3,final_state) ]
Cpi12 = ssevals[0]/np.sqrt(ssevals[3]*ssevals[4])
Cpimod12 = np.abs(Cpi12)*np.ones(len(times))
Cpiphi12 = np.angle(Cpi12)*np.ones(len(times))
Cpi13 = ssevals[1]/np.sqrt(ssevals[3]*ssevals[5])
Cpimod13 = np.abs(Cpi13)*np.ones(len(times))
Cpiphi13 = np.angle(Cpi13)*np.ones(len(times))
Cpi23 = ssevals[2]/np.sqrt(ssevals[4]*ssevals[5])
Cpimod23 = np.abs(Cpi23)*np.ones(len(times))
Cpiphi23 = np.angle(Cpi23)*np.ones(len(times))
pnarr,pinarr = np.linalg.eig(final_state)
pnarr = np.real(pnarr)
pinlist =[]
for j in range(len(pinarr)):
    pinlist.append(Qobj(pinarr[:,j]))
pick = np.random.choice(len(pnarr),p=pnarr) #picking the initial state vector, this starts the Monte Carlo part
psi0 = pinlist[pick]
data = ssesolve(Htot,psi0,times2,sc_ops=c_op, e_ops = e_ops, method="homodyne")
qsave(data.expect,'3QVdPallexpvalsindi')
Cphi12 = data.expect[0][1:len(times)]/np.sqrt(data.expect[3][1:len(times)]*data.expect[4][1:len(times)])
Cphimod12 = np.abs(Cphi12)
Cphiphase12 = np.angle(Cphi12)
Cphi13 = data.expect[1][1:len(times)]/np.sqrt(data.expect[3][1:len(times)]*data.expect[5][1:len(times)])
Cphimod13 = np.abs(Cphi13)
Cphiphase13 = np.angle(Cphi13)
Cphi23 = data.expect[2][1:len(times)]/np.sqrt(data.expect[4][1:len(times)]*data.expect[5][1:len(times)])
Cphimod23 = np.abs(Cphi23)
Cphiphase23 = np.angle(Cphi23)


x1=data.expect[6]
x2=data.expect[7]
x3=data.expect[8]
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1]
pears12 = rx1x2(t,Dt,x1,x2,dt)
pears13 = rx1x2(t,Dt,x1,x3,dt)
pears23 = rx1x2(t,Dt,x2,x3,dt)


fig,ax =plt.subplots()
ax.plot(times[1:],Cphimod12,'b')
ax.plot(times[1:],Cpimod12[1:],'b--')
ax.plot(t,pears12,'orange')
plt.xticks([0, 20/omg[0], 40/omg[0],60/omg[0],80/omg[0],100/omg[0]],["0","20","40","60","80","100"])
ax.set_xlabel("$\omega_{1}t$")
ax.set_ylabel("|$C_{\psi (t)}$|", color='b')
ax.set_yticks([-1.00,-0.75,-0.50,-0.25,0.00,0.25,0.50,0.75,1.00])
ax.set_yticklabels(['-1.00','','-0.50','','0.00','','0.50','','1.00'])
ax2 = ax.twinx()
ax2.plot(times[1:],Cphiphase12,'g')
ax2.plot(times[1:],Cpiphi12[1:],'g--')
ax2.set_ylabel("$\Delta\phi (t)$",color ='g')
ax2.set_yticks([-1.00*np.pi,-0.50*np.pi,0.00,0.50*np.pi,1.00*np.pi])
ax2.set_yticklabels(['$-\pi$','$-\pi /2$','$0$','$\pi /2$','$\pi$'])
plt.show()

fig,ax =plt.subplots()
ax.plot(times[1:],Cphimod13,'b')
ax.plot(times[1:],Cpimod13[1:],'b--')
ax.plot(t,pears13,'orange')
plt.xticks([0, 20/omg[0], 40/omg[0],60/omg[0],80/omg[0],100/omg[0]],["0","20","40","60","80","100"])
ax.set_xlabel("$\omega_{1}t$")
ax.set_ylabel("|$C_{\psi (t)}$|", color='b')
ax.set_yticks([-1.00,-0.75,-0.50,-0.25,0.00,0.25,0.50,0.75,1.00])
ax.set_yticklabels(['-1.00','','-0.50','','0.00','','0.50','','1.00'])
ax2 = ax.twinx()
ax2.plot(times[1:],Cphiphase13,'g')
ax2.plot(times[1:],Cpiphi13[1:],'g--')
ax2.set_ylabel("$\Delta\phi (t)$",color ='g')
ax2.set_yticks([-1.00*np.pi,-0.50*np.pi,0.00,0.50*np.pi,1.00*np.pi])
ax2.set_yticklabels(['$-\pi$','$-\pi /2$','$0$','$\pi /2$','$\pi$'])
plt.show()

fig,ax =plt.subplots()
ax.plot(times[1:],Cphimod23,'b')
ax.plot(times[1:],Cpimod23[1:],'b--')
ax.plot(t,pears23,'orange')
plt.xticks([0, 20/omg[0], 40/omg[0],60/omg[0],80/omg[0],100/omg[0]],["0","20","40","60","80","100"])
ax.set_xlabel("$\omega_{1}t$")
ax.set_ylabel("|$C_{\psi (t)}$|", color='b')
ax.set_yticks([-1.00,-0.75,-0.50,-0.25,0.00,0.25,0.50,0.75,1.00])
ax.set_yticklabels(['-1.00','','-0.50','','0.00','','0.50','','1.00'])
ax2 = ax.twinx()
ax2.plot(times[1:],Cphiphase23,'g')
ax2.plot(times[1:],Cpiphi23[1:],'g--')
ax2.set_ylabel("$\Delta\phi (t)$",color ='g')
ax2.set_yticks([-1.00*np.pi,-0.50*np.pi,0.00,0.50*np.pi,1.00*np.pi])
ax2.set_yticklabels(['$-\pi$','$-\pi /2$','$0$','$\pi /2$','$\pi$'])
plt.show()