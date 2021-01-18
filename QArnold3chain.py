from qutip import * 
import numpy as np
import matplotlib.pyplot as plt

posdelt = np.linspace(0, 12.5, 2501)
negdelt = np.linspace(-12.5,0,2501)
posdeltfull = np.linspace(0, 25, 2501)
negdeltfull = np.linspace(-25,0,2501)

n1 = 2 #number of states 1
n2 = 2 #number of states 2
n3 = 2
gamup = 0.01 # nonlinear pumping rate1
gamdown = 10000 # nonlinear damping rate1
a1 = tensor(tensor(destroy(n1),qeye(n2)),qeye(n3))
a2 = tensor(tensor(qeye(n1),destroy(n2)),qeye(n3))
a3 = tensor(tensor(qeye(n1),qeye(n2)),destroy(n3))
aa1 = a1*a1
aa2 = a2*a2
aa3 = a3*a3
omg1 = 8*np.pi
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
Domgsp = np.linspace(-25,25,500)
Vspplot = np.linspace(0,50,1000)
Vsp = Vspplot*gamup
VClistlist =[]
for V in Vsp:
    VClist =[]
    sqrtV = np.sqrt(V)
    c_op=[c_op0, c_op1, c_op2, c_op3, c_op4, c_op5, sqrtV*a1ma2, sqrtV*a2ma3]
    for Domg in Domgsp:
        omg2 = omg1 +Domg*gamup
        omg3 = omg2 +Domg*gamup
        Htot = omg1a1da1 + omg2*a2da2 + omg3*a3da3
        final_state = steadystate(Htot,c_op)
        pi77ppi88 = final_state[6,6]+final_state[7,7]
        C = (final_state[2,4]+final_state[3,5])/np.sqrt((final_state[4,4]+final_state[5,5]+pi77ppi88)*(final_state[2,2]+final_state[3,3]+pi77ppi88))
        VClist.append(np.abs(C))
    VClistlist.append(VClist)
Cvals = np.array(VClistlist)

qsave(Cvals,'Cmodvalschain')

fig, ax = plt.subplots() 
cont0 = ax.contourf(Domgsp,Vspplot,Cvals,levels=100)
ax.plot(posdelt, 4*posdelt, 'k--');
ax.plot(negdelt, -4*negdelt,'k--');
ax.plot(posdeltfull,1.3333*posdeltfull,'k-');
ax.plot(negdeltfull,-1.3333*negdeltfull,'k-');
ax.set_aspect(1)
ax.set_title("Quantum Arnold tongue");
ax.set_xlabel("$\Delta\omega_{2,1}$"+" (orders of "+"$\gamma$"+")");
ax.set_ylabel("$V$"+" (orders of "+"$\gamma$"+")")
cb =fig.colorbar(cont0)
plt.show()