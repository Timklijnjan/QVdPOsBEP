import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import pickle


#Defining system parameters
n = 2 #number of states per oscillator
omg1 = 8*np.pi
gamup = 0.01
gamdown = 10000
#theta = 0, implementation of complex powers of e is required when a shifted coupling is used
Vlist = [100*gamup,5*gamup,20*gamup] #0.1*gam1u, 0.215*gam1u,0.464*gam1u,gam1u,2.15*gam1u,4.64*gam1u,
Vstr = ['100','5','20'] # To label files
#omg =[omg1-0*gam1u, omg1+3*gam1u, omg1+10*gam1u, omg1+17*gam1u]
Domgint = [gamup,gamup,20*gamup]
Domgstr = ['1','1','20']
Dt = 80*np.pi/omg1 #window used for pearson indicator

#Setting timesteps used for evaluation
times = np.linspace(0.0,1000/omg1, 10000)
dt = times[1]-times[0]
dtt = 1/dt
extime = math.ceil(Dt/dt)+2
times2 = np.linspace(0.0,times[len(times)-1]+extime*dt,len(times)+extime) #array used to obtain a large enough interval for the pearson indicator integration on all of times
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1] #times for which pearson indicator is calculated

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

#thetvals =[]
for i in range(len(Vlist)):
    V = Vlist[i]
    Domg = Domgint[i]
    omg = [omg1,omg1+Domg,omg1+2*Domg]
    sqrtV = np.sqrt(V)
    #Creation of the Lindblad operators
    c_op=[c_op0, c_op1, c_op2, c_op3, c_op4, c_op5, sqrtV*a1ma2, sqrtV*a2ma3, sqrtV*a1ma3]

    #Creation of the system Hamiltonian
    Htot = omg1a1da1 + omg[1]*a2da2 + omg[2]*a3da3
    
    #Calculation of the steady state to obtain samples of the initial wave function from
    #final_state = steadystate(Htot,c_op)
    finalinitialstatename = './finalstateV'+Vstr[i]+'D'+Domgstr[i]
    final_state = qload(finalinitialstatename) #qsave(final_state, finalinitialstatename), swap this and activate 2 lines above to generate steady states (done independently since steadystate is broken 4/5 of the time)
    #Setting up initial wave functions sampled from the steady state with their respective probabilities
    pnarr,pinarr = np.linalg.eig(final_state)
    pnarr = np.real(pnarr)
    pinlist =[]
    for j in range(len(pinarr)):
        pinlist.append(Qobj(pinarr[:,j]))


    #expvals =[]
    nummont = 1000 #number of trajectories used for Monte Carlo calculations
    temp = []
    
    for m in range(nummont):
        pick = np.random.choice(len(pnarr),p=pnarr) #picking the initial state vector, this starts the Monte Carlo part
        psi0 = pinlist[pick]
        data = ssesolve(Htot,psi0,times2,sc_ops=c_op, method="homodyne")
        #expvals.append(data.states)
        temp.append(data.states)
        tn=m+1

        if (tn % 100)==0:            
            filename = './tempory3osV'+Vstr[i]+'D'+Domgstr[i] +'_state_dumpp_'+str(tn)
            qsave(temp,filename)
            temp = []

    #thetvals.append(expvals)
Datlist =[len(times),len(times2),dt,Dt,Vlist,Domgint,nummont]
qsave(Datlist,filename+"param")
#Datlist.append(thetvals) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

#qsave(Datlist,"Vvars5")'''

