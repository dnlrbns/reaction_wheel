import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m1= 0.1
m2 = 0.1
r = 0.04
l1 = 0.1
l2 = 0.2
I1 = (1/12) * m1 * l1**2
I2 = (1/2)  * m2 * r**2
b = 0.01
c = 0.001
T=0
g=9.8
ref = -np.pi/2

def dZ_dt(t, Z):
    dZ = [None, None, None, None, None, None]
    dZ[0] = Z[1]
    n1 = (m1*g*l1 - m2*g*l2)*np.cos(Z[0]) \
        - b*Z[1] \
        + c*Z[3] \
        - T
    d1 = m1*l1**2 + m2*l2**2 + I1
    dZ[1] = n1/d1
    dZ[2] = Z[3]
    n2 = T \
        - c*Z[3] \
        - dZ[1]*I2
    d2 = I2
    dZ[3] = n2/d2
    dZ[4] = ref - Z[0]
    dZ[5] = -Z[1]
    return dZ

sol = solve_ivp(dZ_dt, [0,4], [-np.pi/2,0,0,20*np.pi,ref+np.pi/2,0], max_step=.05)

try:
    plt.close('all')
except:
    print('no figures to close')

fig,[ax0,ax1,ax2,ax3]= plt.subplots(4,1,sharex=True)
l0, = ax0.plot(sol.t,np.rad2deg(sol.y[[0],:].T))
ax0.grid()
l1, = ax1.plot(sol.t,sol.y[[1],:].T/2/np.pi)
ax1.grid()
l2, = ax2.plot(sol.t,sol.y[[3],:].T/2/np.pi)
ax2.grid()
l3, = ax3.plot(sol.t,np.rad2deg(sol.y[[4],:].T))
ax3.grid()
fig.legend()
plt.show()
