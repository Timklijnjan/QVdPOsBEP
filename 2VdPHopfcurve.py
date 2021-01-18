import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import newton
import matplotlib.pyplot as plt

def somg(mu1, sig, omg):
    return sig**2/((sig-omg**2)*mu1)

def frapre(mu2, mu1,omg1s,omg2s): 
    return (mu1*omg2s+mu2*omg1s)/(mu1+mu2)

def sighop(mu1, mu2, omg1s, omg2s):
    x = 1/(mu1+mu2)
    z = x**2
    a = 1*np.ones(10000)
    b = -x*(mu1**2*mu2+mu1*mu2**2)
    c = z*(mu1*mu2*(omg1s-omg2s)**2)+x*(mu1**2*mu2*omg2s+mu1*mu2**2*omg1s)
    Det = b**2-4*a*c
    up = []
    down = []
    for i in range(len(Det)):
        if Det[i] >= 0:
            up.append((-b[i]+np.sqrt(Det[i]))/(2*a[i]))
            down.append((-b[i]-np.sqrt(Det[i]))/(2*a[i]))
        else:
            up.append(complex((-b[i])/(2*a[i]), np.sqrt(-Det[i])/(2*a[i])))
            down.append(complex((-b[i])/(2*a[i]), -np.sqrt(-Det[i])/(2*a[i])))
    return [up, down]

def calceig(mu1, mu2, sig, omg1s, omg2s):
    A = np.array([[0,1,0,0],[sig-omg1s,mu1,-sig,0],[0,0,0,1],[-sig,0,sig-omg2s,mu2]])
    return np.linalg.eig(A)[0]

def poly(lamb, mu1, mu2, sig, omg1s, omg2s):
    return lamb**4 -(mu1+mu2)*lamb**3+(mu1*mu2-2*sig+omg1s +omg2s)*lamb**2 +(mu1*sig-mu1*omg2s-mu2*omg1s+mu2*sig)*lamb-sig*(omg1s+omg2s)+omg1s*omg2s

def deriv(lamb, mu1, mu2, sig, omg1s, omg2s):
    return 4*lamb**3 -3*(mu1+mu2)*lamb**2+2*(mu1*mu2-2*sig+omg1s+omg2s)*lamb+mu1*sig-mu1*omg2s-mu2*omg1s+mu2*sig

sig = 0.4
omg =1
omg1 = 1
omg2 = 1
omg1s = omg1**2*np.ones(10000)
#omg2s = omg2**2
omg2s = omg2**2*np.ones(10000)
#fromg = np.linspace(0.1, 50, 1000)
#mu1 = np.linspace(-50,50,10000)
#mu2 = 1*np.ones(10000)
mu1 = np.linspace(0.001,0.25, 2000)
mu1n = np.linspace(-0.25, -0.001, 2000)

mu2 = somg(mu1,sig,omg)
mu2n =somg(mu1n,sig,omg)

fig, ax = plt.subplots()
ax.plot(mu1, mu2,'b');
ax.plot(mu1n,mu2n,'b');
ax.set_xlabel('$\mu_1$');
ax.set_ylabel('$\mu_2$');
plt.show(fig)

"""mu2c = fsolve(domg, np.ones(len(mu1))*0.1, (mu1, sig*np.ones(len(mu1)),omg1s*np.ones(len(mu1)),omg2s*np.ones(len(mu1))))
mu2cn = fsolve(domg, np.ones(len(mu1n))*0.1, (mu1n, sig*np.ones(len(mu1n)),omg1s*np.ones(len(mu1n)),omg2s*np.ones(len(mu1n))))

fig, ax = plt.subplots()
ax.plot(mu1, mu2c,'b');
ax.plot(mu1n,mu2cn,'b');
ax.set_xlabel('$\mu_1$');
ax.set_ylabel('$\mu_2$');
plt.show(fig)"""

"""sigsol = sighop(mu1, mu2, omg1s, omg2s)
feas = frapre(mu2,mu1,omg1s,omg2s)
upsol = sigsol[0]
downsol = sigsol[1]
mu1f = []
upsolf =[]
downsolf=[]
for j in range(len(upsol)):
    if isinstance(upsol[j],float):
        mu1f.append(mu1[j])
        upsolf.append(upsol[j])
        downsolf.append(downsol[j])

fig, ax = plt.subplots()
ax.plot(mu1f, upsolf,'b');
ax.plot(mu1f, downsolf, 'r');
#ax.plot(mu1,midsol,'y')
ax.plot(mu1,feas,'g');
ax.set_xlabel('$\mu_{1}$');
ax.set_ylabel('$\sigma$');
plt.show(fig)

emu1= mu1f[0]
esig = downsolf[0]
eigval =calceig(emu1,mu2[0],esig, omg1s[0],omg2s[0])
print(str(eigval))

#lamb = complex(0,-9999.8000135)
#print(poly(lamb,mu1,mu2,esig,omg1s,eomg2s))"""