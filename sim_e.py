import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m_a = 0.1
m_d = 0.2
r = 0.05
l = 0.2
I_a = (1/12) * m_a * l**2
I_d = (1/2)  * m_d * r**2
b = 1
c = 0.01
T=0
g=9.8
rpm =30*0.1047

ref= np.deg2rad(90)

def T_law(th,ei):
    Kp =-10
    Ki = -15
    e =ref -th
    return Kp*e + Ki*ei
def dX_dt(t, X):
    if X[3]>rpm:
        X[3] = rpm
    elif X[3]<-rpm:
        X[3] = -rpm
    else:
        pass
    T =T_law(X[0],X[4])
    
    num=b*X[3]\
        -c*X[2]\
        -(m_a+2*m_d)*g*l*np.cos(X[0])\
        -T
    den = (m_a+4*m_d)*l**2 + I_a
    ddth = num/den
    
    ddph = (T-I_d*ddth-b*X[3])/I_d
    
    e = ref - X[0]
    
    return [X[2],X[3], ddth, ddph,e]

sol = solve_ivp(dX_dt, [0,10], [np.deg2rad(93),0,0,0,0],max_step=0.01)

try:
    plt.close('all')
except:
    print('no figures to close')
    
fig, (ax1, ax2) = plt.subplots(2,1)

ax1.plot(sol.t,np.rad2deg(sol.y[0,:]))
ax1.grid()
ax2.plot(sol.t,sol.y[3,:]/0.1047)
ax2.grid()
plt.show()
