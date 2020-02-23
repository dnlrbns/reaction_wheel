import numpy as np
import matplotlib.pyplot as plt
import control

m1= 0.1
m2 = 0.1
r = 0.04
l1 = 0.1
l2 = 0.2
I1 = (1/12) * m1 * l1**2
I2 = (1/2)  * m2 * r**2
b = 0.01
c = 0.01
T=0
g=9.8
ref = np.pi/2

J = m1*l1**2 + m2*l2**2 + I1

A = [\
    [-b/J, (g/J)*(m1*l1+m2*l2), c/J, 0],\
    [1, 0, 0, 0],\
    [b/J, (-g/J)*(m1*l1+m2*l2), -c/J, 0],\
    [0, 0, 1, 0]]
B = [[-1/J],[0],[1/(J+I2)],[0]]

C = [0,1,0,0]
D = 0
sys = control.StateSpace(A,B,C,D)
Q = np.eye(4)
R = np.eye(1)

print(Q)
print(R)
print(sys)
K = control.lqr(A,B,Q,R)
print(K)

# K, S, E = lqr(A,B, Q, R, [N])

# plt.plot(control.impulse_response(sys))
# plt.show()