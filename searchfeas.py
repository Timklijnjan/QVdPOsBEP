import numpy as np 

def Detomg2s(mu1, mu2, omg1s):
    return (8*mu1*mu2*omg1s-4*mu1**2*mu2**2-4*mu1**3*mu2)**2+16*mu1*mu2*(mu1**4*mu2**2+2*mu1**3*mu2**3+mu1**2*mu2**4-4*mu1*mu2*omg1s**2-4*mu1**2*mu2**2*omg1s-4*mu1*mu2**3*omg1s)

def omg2scalc(mu1, mu2, omg1s, Det):
    return[(-(8*mu1*mu2*omg1s-4*mu1**2*mu2**2-4*mu1**3*mu2)+Det)/(2*-4*mu1*mu2), (-(8*mu1*mu2*omg1s-4*mu1**2*mu2**2-4*mu1**3*mu2)-Det)/(2*-4*mu1*mu2)]
mu1 = 0.2
mu2= 0.3
omg1s = np.linspace(0.01, 10, 100)

Detomg2ssol = Detomg2s(mu1,mu2,omg1s)
print(max(Detomg2ssol))
indexmax = np.argmax(max(Detomg2ssol))
print(omg2scalc(mu1,mu2,0.01,Detomg2ssol[0]))