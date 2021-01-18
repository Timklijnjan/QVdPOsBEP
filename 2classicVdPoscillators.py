import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import solve_ivp

def VdP(t,z):
	x1, y1, x2, y2 = z
	return [y1, mu1*(1-x1**2)*y1-omg1**2*x1 -sig*(x2-x1), y2, mu2*(1-x2**2)*y2-omg2**2*x2 -sig*(x1-x2)] 


t0 = 0
t1 = 150
x10 = 1
dx10 = 0
x20 = 1
dx20 = 0

mu1 = -0.2
mu2 = 0.5
sig = 0.1
omg1 = 1
omg2 = 1

t = scipy.linspace(t0,t1,1000000)

sol = solve_ivp(VdP, [t0, t1], [x10, dx10, x20, dx20 ], t_eval = t)

plt.plot(sol.y[0],sol.y[1], ls = '--')
plt.plot(sol.y[2],sol.y[3], ls = ':')
plt.legend(["$x_1$","$x_2$"])
plt.title("phase plot of oscillator 1 and 2")
plt.xlabel('x')
plt.ylabel("$\dot{x}$")
plt.axes().set_aspect(1)
plt.show()


"""plt.plot(sol.t, sol.y[0])
plt.plot(sol.t, sol.y[2])

plt.legend(["$x_1$","$x_2$"])
plt.title('x solution')
plt.xlabel('t')
plt.ylabel('x')
plt.show()

plt.plot(sol.t, sol.y[0]-sol.y[2])
plt.show()"""