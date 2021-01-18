from qutip import * 
import numpy as np
import matplotlib.pyplot as plt

n = 5 #number of states
omg = 1 #frequency of oscillator
gamup = 0.1 # nonlinear pumping rate
gamdown = 100 # nonlinear damping rate
a = destroy(n)
aa = a*a
Htot = omg*a.dag()*a
fopu = 0
for i in range(n):
    fopu = fopu + fock(n,i)
wave = 1/np.sqrt(n)*fopu
psi0 = fock(n,0) #wave*wave.dag()
c_op=[np.sqrt(gamup)*a.dag(), np.sqrt(gamdown)*aa]
times = np.linspace(0.0, 50, 10000)
result = mesolve(Htot,psi0,times,c_op)
final_state = steadystate(Htot,c_op)

pr0 =[]
pr1 =[]
""" pr2 = []
pr3 =[]
pr4=[]
pr5=[]
pr6=[]
pr7=[]
pr8=[] """
for t in range(0,len(times)):
    pr0.append(result.states[t][0,0])
    pr1.append(result.states[t][1,1])
    """ pr2.append(result.states[t][2,2])
    pr3.append(result.states[t][3,3])
    pr4.append(result.states[t][4,4])
    pr5.append(result.states[t][5,5])
    pr6.append(result.states[t][6,6])
    pr7.append(result.states[t][7,7])
    pr8.append(result.states[t][8,8]) """

""" print(1-(final_state[0,0]+final_state[1,1]))
print(1-min(np.array(pr0)+np.array(pr1)+np.array(pr2)+np.array(pr3)+np.array(pr4)+np.array(pr5)+np.array(pr6)+np.array(pr7)+np.array(pr8))) # """
fig, ax = plt.subplots()
ax.plot(times, pr0);
ax.plot(times, pr1);
""" ax.plot(times, pr2);
ax.plot(times, pr3);
ax.plot(times, pr4);
ax.plot(times, pr5);
ax.plot(times, pr6);
ax.plot(times, pr7);
ax.plot(times, pr8); """
ax.plot(times, [final_state[0,0]]*np.ones(len(times)),'r');
ax.plot(times, [final_state[1,1]]*np.ones(len(times)),'r');
ax.set_xlabel('Time');
ax.set_ylabel('Pure state probability');
ax.legend(("$P_0$","$P_1$","Steady state")) #,"$P_2$","$P_3$","$P_4$","$P_5$","$P_6$","$P_7$","$P_8$"
plt.show(fig)

""" xvec = np.linspace(-5,5,500)
W_fock = wigner(result.states[len(times)-1],xvec,xvec)
fig, axes = plt.subplots(1,1,figsize=(12,12))
cont0 = axes.contourf(xvec,xvec,W_fock,100)
lbl0 = axes.set_title("Wigner function")
axes.set_xlabel('$x/x_{zpf}$')
axes.set_ylabel('$p/p_{zpf}$')
cb = fig.colorbar(cont0)
plt.show() """