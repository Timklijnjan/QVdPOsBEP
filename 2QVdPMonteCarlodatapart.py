import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import scipy as scipy
import math
import pickle

#Defining system parameters
n = 2 #number of states per oscillator
omg1 = 8*np.pi
gam1u = 0.01
gam2u = gam1u
gam1d = 10000
gam2d = gam1d
theta = 0
Vlist = [5*gam1u] #0.1*gam1u, 0.215*gam1u,0.464*gam1u,gam1u,2.15*gam1u,4.64*gam1u,
omg =[omg1,omg1+1*gam1u]
Domg = omg[1]-omg[0]
Dt = 80*np.pi/omg1 #window used for pearson indicator

#Setting timesteps used for evaluation
times = np.linspace(0.0,1000/omg[0], 20000)
dt = times[1]-times[0]
dtt = 1/dt
extime = math.ceil(Dt/dt)+2
times2 = np.linspace(0.0,times[len(times)-1]+extime*dt,len(times)+extime) #array used to obtain a large enough interval for the pearson indicator integration on all of times
t= times2[extime:math.floor((times2[len(times2)-1]-Dt)/dt)-1] #times for which pearson indicator is calculated

thetvals =[]
for V in Vlist:
    #Setting up initial wave functions sampled from the steady state with their respective probabilities
    N = (3*gam1u+V)*(3*gam1u*(Domg**2+9*gam1u**2)+(Domg**2+27*gam1u**2)*V+8*gam1u*V**2)
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

    #defining basic quantum operators
    a1 = tensor(destroy(n),qeye(n))
    a1d = a1.dag()
    a2 = tensor(qeye(n),destroy(n))
    a2d = a2.dag()
    a=[a1,a2]
    ad = [a1d,a2d]
    H = 0 #Creation of the system Hamiltonian
    for i in range(len(omg)):
        H = H + omg[i]*ad[i]*a[i]
    L = [np.sqrt(gam1d)*a1**2,np.sqrt(gam1u)*a1d, np.sqrt(gam2d)*a2**2, np.sqrt(gam2u)*a2d, np.sqrt(V)*(a1-np.exp(complex(0,theta))*a2)] #Lindblad operators
    e_ops = [a1d*a2,a1d*a1,a2d*a2,(a1+a1d)/np.sqrt(2),(a2+a2d)/np.sqrt(2)] # Operators for which expectation values are determined

    
    expvals =[]
    nummont = 1 #number of trajectories used for Monte Carlo calculations
    for m in range(nummont):
        pick = np.random.choice(4,p=np.array([p1,p2,p3,p4])) #picking the initial state vector, this starts the Monte Carlo part
        psi0 = pn[pick]
        data = ssesolve(H,psi0,times2,sc_ops=L, e_ops = e_ops, method="homodyne")
        expv1=data.expect[0]
        expv2=data.expect[1]
        expv3=data.expect[2]
        expv4=data.expect[3]
        expv5=data.expect[4]
        expvals.append([expv1,expv2,expv3,expv4,expv5])
    thetvals.append(expvals)
Datlist =[len(times),len(times2),dt,Dt,Vlist,nummont]
Datlist.append(thetvals) #[len(times),len(times2),dt,Dt,[5*gamu1],nummont,[V0:[trja0:[expv1,expv2,expv3,expv4,expv5],trja1:[...],...],V1:[[]]...]]

qsave(Datlist,"Vvars5")

