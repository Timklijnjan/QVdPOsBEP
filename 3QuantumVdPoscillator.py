from qutip import * 
import numpy as np
import matplotlib.pyplot as plt

n1 = 2 #number of states 1
n2 = 2 #number of states 2
n3 = 2
gamup = 0.01 # nonlinear pumping rate1
omg1 = 8*np.pi #frequency of oscillator1
omg2 = omg1+5*gamup #frequency of oscillator2
omg3 = omg1+10*gamup
gamdown = 10000 # nonlinear damping rate1
a1 = tensor(tensor(destroy(n1),qeye(n2)),qeye(n3))
a2 = tensor(tensor(qeye(n1),destroy(n2)),qeye(n3))
a3 = tensor(tensor(qeye(n1),qeye(n2)),destroy(n3))
aa1 = a1*a1
aa2 = a2*a2
aa3 = a3*a3
theta = 0
V =50
Htot = omg1*a1.dag()*a1 +omg2*a2.dag()*a2+omg3*a3.dag()*a3
fopu1 = 0
for i in range(n1):
    fopu1 = fopu1 + fock(n1,i)
wave1 = 1/np.sqrt(n1)*fopu1
fopu2 = 0
for i in range(n2):
    fopu2 = fopu2 + fock(n2,i)
wave2 = 1/np.sqrt(n2)*fopu2
fopu3 = 0
for i in range(n3):
    fopu3 = fopu3 + fock(n3,i)
wave3 = 1/np.sqrt(n3)*fopu3
psi0 = tensor(tensor(wave1*wave1.dag(),wave2*wave2.dag()),wave3*wave3.dag())
c_op=[np.sqrt(gamup)*a1.dag(), np.sqrt(gamdown)*aa1, np.sqrt(gamup)*a2.dag(), np.sqrt(gamdown)*aa2, np.sqrt(gamup)*a3.dag(),np.sqrt(gamdown)*aa3, np.sqrt(V)*(a1-np.exp(complex(0,theta))*a2),np.sqrt(V)*(a2-np.exp(complex(0,theta))*a3),np.sqrt(V)*(a1-np.exp(complex(0,theta))*a3)]
times = np.linspace(0.0, 500.0, 10000)
#result = mesolve(Htot,psi0,times,c_op)
final_state = steadystate(Htot,c_op)
#qsave(result,'3QVdPtoss')

pr000 =[]
pr001 =[]
pr010 = []
pr011 =[]
pr100=[]
pr101=[]
pr110=[]
pr111=[]




'''for t in range(0,len(times)):
    pr000.append(result.states[t][0,0])
    pr001.append(result.states[t][1,1])
    pr010.append(result.states[t][2,2])
    pr011.append(result.states[t][3,3])
    pr100.append(result.states[t][4,4])
    pr101.append(result.states[t][5,5])
    pr110.append(result.states[t][6,6])
    pr111.append(result.states[t][7,7])



#print(1-(final_state[0,0]+final_state[1,1]))
#print(1-min(np.array(pr00[len(times)-2:len(times)-1])+np.array(pr01[len(times)-2:len(times)-1])+np.array(pr02[len(times)-2:len(times)-1])+np.array(pr10[len(times)-2:len(times)-1])+np.array(pr11[len(times)-2:len(times)-1])+np.array(pr20[len(times)-2:len(times)-1]))) #+np.array(pr4)+np.array(pr5)+np.array(pr6)+np.array(pr7)+np.array(pr8)"""
fig, ax = plt.subplots()
ax.plot(times, pr000);
ax.plot(times, pr001);
ax.plot(times, pr010);
ax.plot(times, pr011);
ax.plot(times, pr100);
ax.plot(times, pr101);
ax.plot(times, pr110);
ax.plot(times, pr111);'''
'''ax.plot(times, pr20);
ax.plot(times, pr21);
ax.plot(times, pr22);
ax.plot(times, pr23);
ax.plot(times, pr30);
ax.plot(times, pr31);
ax.plot(times, pr32);
ax.plot(times, pr33);'''

labels =['pr000','pr001','pr010','pr011','pr100','pr101','pr110','pr111']
for j in range(n1*n2*n3):
    ax.plot(times,[final_state[j,j]]*np.ones(len(times)))
    labels.append("P"+str(j))
ax.set_xlabel('Time');
ax.set_ylabel('Pure state probabilityss');
ax.legend(labels)
plt.show(fig) 

'''Findens = result.states[len(times)-1]
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
plt.show()'''