import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m_a = 0.1
m_d = 0.5
r = 0.05
l1 = 0.1
l2 = 0.2
I_a = (1/12) * m_a * l1**2
I_d = (1/2)  * m_d * r**2
b = 0.001
c = 0.01
T=0
g=9.8
ref = 0

def dZ_dt(t, Z):
    dZ = [None, None, None, None, None, None]
    dZ[0] = Z[1]
    n1 = (m_a*g*l1 - m_d*g*l2)*np.cos(Z[0]) - b*Z[1] + c*Z[3] - T
    d1 = m_a*l1 + m_d*l2 + I_a
    dZ[1] = n1/d1
    dZ[2] = Z[3]
    n2 = T - c*Z[3] - dZ[1]*I_d + m_d*g*l2*np.cos(Z[0])
    d2 = I_d
    dZ[3] = n2/d2
    dZ[4] = ref - Z[0]
    dZ[5] = -Z[1]
    dZ = [None, None, None, None, None, None]
    return dZ

sol = solve_ivp(dZ_dt, [0,10], [-np.pi/2,0,0,0,0,0])

plt.plot(sol.t,np.rad2deg(np.transpose(sol.y)))
plt.show()
