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
V = 20*gamup #[100*gamup, 5*gamup,20*gamup] #0.1*gam1u, 0.215*gam1u,0.464*gam1u,gam1u,2.15*gam1u,4.64*gam1u,
Vstr = '20'# ['100','5','20'] # To label files
#omg =[omg1-0*gam1u, omg1+3*gam1u, omg1+10*gam1u, omg1+17*gam1u]
Domg = 20*gamup #[gamup,gamup,20*gamup]
Domgstr = '20'#['1','1','20']
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

omg = [omg1,omg1+Domg,omg1+2*Domg]
sqrtV = np.sqrt(V)
#Creation of the Lindblad operators
c_op=[c_op0, c_op1, c_op2, c_op3, c_op4, c_op5, sqrtV*a1ma2, sqrtV*a2ma3]

#Creation of the system Hamiltonian
Htot = omg1a1da1 + omg[1]*a2da2 + omg[2]*a3da3

final_state = steadystate(Htot,c_op)
qsave(final_state,'./finalstatecircV'+Vstr+'D'+Domgstr)