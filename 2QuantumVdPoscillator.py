from qutip import * 
import numpy as np
import matplotlib.pyplot as plt

n1 = 4 #number of states 1
n2 = 4 #number of states 2
omg1 = 1010 #frequency of oscillator1
omg2 = 1000 #frequency of oscillator2
gamup1 = 0.4 # nonlinear pumping rate1
gamup2 = 0.1 # nonlinear pumping rate2
gamdown1 = 1000 # nonlinear damping rate1
gamdown2 = 1000 # nonlinear damping rate2
a1 = tensor(destroy(n1),qeye(n2))
a2 = tensor(qeye(n1),destroy(n2))
aa1 = a1*a1
aa2 = a2*a2
theta = 0.5*np.pi
V =0.1
Htot = omg1*a1.dag()*a1 +omg2*a2.dag()*a2
fopu1 = 0
for i in range(n1):
    fopu1 = fopu1 + fock(n1,i)
wave1 = 1/np.sqrt(n1)*fopu1
fopu2 = 0
for i in range(n2):
    fopu2 = fopu2 + fock(n2,i)
wave2 = 1/np.sqrt(n2)*fopu2
psi0 = tensor(wave1*wave1.dag(),wave2*wave2.dag())
c_op=[np.sqrt(gamup1)*a1.dag(), np.sqrt(gamdown1)*aa1, np.sqrt(gamup2)*a2.dag(), np.sqrt(gamdown2)*aa2, np.sqrt(V)*(a1-np.exp(complex(0,theta))*a2)]
times = np.linspace(0.0, 50.0, 10000)
result = mesolve(Htot,psi0,times,c_op)
final_state = steadystate(Htot,c_op)

pr00 =[]
pr01 =[]
pr02 = []
pr03 = []
pr10 =[]
pr11=[]
pr12=[]
pr13=[]
pr20=[]
pr21=[]
pr22=[]
pr23=[]
pr30=[]
pr31=[]
pr32=[]
pr33=[]

for t in range(0,len(times)):
    pr00.append(result.states[t][0,0])
    pr01.append(result.states[t][1,1])
    pr02.append(result.states[t][2,2])
    pr03.append(result.states[t][3,3])
    pr10.append(result.states[t][n2,n2])
    pr11.append(result.states[t][n2+1,n2+1])
    pr12.append(result.states[t][n2+2,n2+2])
    pr13.append(result.states[t][n2+3,n2+3])
    pr20.append(result.states[t][2*n2,2*n2])
    pr21.append(result.states[t][2*n2+1,2*n2+1])
    pr22.append(result.states[t][2*n2+2,2*n2+2])
    pr23.append(result.states[t][2*n2+3,2*n2+3])
    pr30.append(result.states[t][3*n2,3*n2])
    pr31.append(result.states[t][3*n2+1,3*n2+1])
    pr32.append(result.states[t][3*n2+2,3*n2+2])
    pr33.append(result.states[t][3*n2+3,3*n2+3])


#print(1-(final_state[0,0]+final_state[1,1]))
print(1-min(np.array(pr00[len(times)-2:len(times)-1])+np.array(pr01[len(times)-2:len(times)-1])+np.array(pr02[len(times)-2:len(times)-1])+np.array(pr10[len(times)-2:len(times)-1])+np.array(pr11[len(times)-2:len(times)-1])+np.array(pr20[len(times)-2:len(times)-1]))) #+np.array(pr4)+np.array(pr5)+np.array(pr6)+np.array(pr7)+np.array(pr8)"""
fig, ax = plt.subplots()
ax.plot(times, pr00);
ax.plot(times, pr01);
ax.plot(times, pr02);
ax.plot(times, pr03);
ax.plot(times, pr10);
ax.plot(times, pr11);
ax.plot(times, pr12);
ax.plot(times, pr13);
ax.plot(times, pr20);
ax.plot(times, pr21);
ax.plot(times, pr22);
ax.plot(times, pr23);
ax.plot(times, pr30);
ax.plot(times, pr31);
ax.plot(times, pr32);
ax.plot(times, pr33);

ax.plot(times, [final_state[0,0]]*np.ones(len(times)),'r');
ax.plot(times, [final_state[1,1]]*np.ones(len(times)),'r');
ax.set_xlabel('Time');
ax.set_ylabel('Pure state probability');
ax.legend(("$P_{00}$","$P_{01}$","$P_{02}$","$P_{03}$","$P_{10}$","$P_{11}$","$P_{12}$","$P_{13}$","$P_{20}$","$P_{21}$","$P_{22}$","$P_{23}$","$P_{30}$","$P_{31}$","$P_{32}$","$P_{33}$","Steady state"))
plt.show(fig) 

Findens = result.states[len(times)-1]
Fintr1 = Findens.ptrace(0)
Fintr2 =Findens.ptrace(1)

xvec = np.linspace(-5,5,500)
W_fock = wigner(Findens,xvec,xvec)
fig, axes = plt.subplots(1,1,figsize=(12,12))
cont0 = axes.contourf(xvec,xvec,W_fock,100)
lbl0 = axes.set_title("Wigner function full system")
axes.set_xlabel('x?')
axes.set_ylabel('p?')
cb = fig.colorbar(cont0)
plt.show()

xvec = np.linspace(-5,5,500)
W_fock = wigner(Fintr1,xvec,xvec)
fig, axes = plt.subplots(1,1,figsize=(12,12))
cont0 = axes.contourf(xvec,xvec,W_fock,100)
lbl0 = axes.set_title("Wigner function system1")
axes.set_xlabel('x?')
axes.set_ylabel('p?')
cb = fig.colorbar(cont0)
plt.show()

xvec = np.linspace(-5,5,500)
W_fock = wigner(Fintr2,xvec,xvec)
fig, axes = plt.subplots(1,1,figsize=(12,12))
cont0 = axes.contourf(xvec,xvec,W_fock,100)
lbl0 = axes.set_title("Wigner function system2")
axes.set_xlabel('x?')
axes.set_ylabel('p?')
cb = fig.colorbar(cont0)
plt.show()