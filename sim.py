import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m_a = 0.1
m_d = 0.5
r = 0.05
l = 0.1
I_a = (1/12) * m_a * l**2
I_d = (1/2)  * m_d * r**2
b = 0.1
c = 0.01
T=0
g=9.8

def T_law():
    return 0
def dX_dt(t, X):
    T =T_law()
    
    num=b*X[3]\
        -c*X[2]\
        -(m_a+2*m_d)*g*l*np.cos(X[0])\
        -T
    den = (m_a+4*m_d)*l**2 + I_a
    ddth = num/den
    ddph = (T-I_d*ddth-b*X[3])/I_d
    
    return [X[2],X[3], ddth, ddph]

sol = solve_ivp(dX_dt, [0,10], [0,0,0,0],max_step=0.01)

fig = plt.figure()
ax = plt.axes()
ax.plot(sol.t,np.rad2deg(np.transpose(sol.y[0:2,:])))
plt.show()
