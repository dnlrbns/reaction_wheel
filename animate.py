import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib import animation

m_a = 0.1
m_d = 0.5
r = 0.05
l = 0.1
I_a = (1/12) * m_a * l**2
I_d = (1/2)  * m_d * r**2
b = 0.001
c = 0.01
T=0
g=9.8

def dX_dt(t, X):
    ddth = (b*X[3]-c*X[2]-(m_a+2*m_d)*g*l*np.cos(X[0])-T)/((m_a+4*m_d)*l**2 + I_a)
    ddph = (T-I_d*ddth-b*X[3])/I_d
    
    return [X[2],X[3], ddth, ddph]

fps = 30
time = 10
dt = 1/fps
steps = fps*time
th_list = np.empty([steps,],dtype=np.float)
th_list.fill(np.nan)
states = [-np.pi/2+np.deg2rad(1),0,0,0]
for i in range(steps):
    sol = solve_ivp(dX_dt, [i*dt,(i+1)*dt], states)
    th_list[i] = sol.y[0][-1]
    states = [sol.y[0][-1],sol.y[1][-1],sol.y[2][-1],sol.y[3][-1]]
    print(states)
print(th_list)
plt.plot(range(steps),th_list)
#plt.plot(sol.t,np.rad2deg(np.transpose(sol.y)),'*')
plt.show()
