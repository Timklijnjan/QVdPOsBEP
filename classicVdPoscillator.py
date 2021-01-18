import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import solve_ivp

def VdP(t,z):
	x, y = z
	return [y, mu*(1-x**2)*y-omg**2*x]

def VdP2(X, t = 0):
	return np.array([X[1], mu*(1-X[0]**2)*X[1]-omg**2*X[0]])


t0 = 0
t1 = 50
x0 = 1
dx0 = 0
mu = 200
omg = 10

xg = np.linspace(-3,3,60)
yg = np.linspace(-3,3,60)
X1, Y1 = np.meshgrid(xg,yg)
DX1, DY1 = VdP2([X1,Y1])
M = (np.hypot(DX1,DY1))
M[ M==0] = 1
DX1 /= M
DY1 /= M


t = scipy.linspace(t0,t1,10000)


sol = solve_ivp(VdP, [t0, t1], [x0, dx0], t_eval = t)
fig1 = plt.figure(figsize=(8,6))
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(sol.y[0],sol.y[1])

#muval = [1,1,1]
#for mu in muval:
#	sol = solve_ivp(VdP, [t0, t1], [x0, dx0], t_eval = t)
#	plt.plot(sol.y[0],sol.y[1])

plt.xlim([-3,3])
plt.ylim([-3,3])
plt.axes().set_aspect(1)
plt.title('phase diagram')
plt.xlabel('x')
plt.ylabel("$\dot{x}$")
ax1.quiver(X1, Y1, DX1, DY1, M, pivot ='mid')

plt.show()

#for mu in muval:
	#sol = solve_ivp(VdP, [t0, t1], [x0, dx0], t_eval = t)
	#plt.plot(sol.t, sol.y[0])

sol = solve_ivp(VdP, [t0, t1], [x0, dx0], t_eval = t)
plt.plot(sol.t, sol.y[0])
plt.title('x solution')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
