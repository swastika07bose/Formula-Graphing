# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 21:01:14 2023

@author: USER
"""


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Assign variables
a1=5.0
a2=0.1
b1=3.0
b2=2.0
c=0.85
d1=0.4
d2=0.01
alpha=2.0
zai=0.0
N1init=0.75
N2init=0.3
N3init=2.0


#define the function
def dN1dt(t,N1,N2,N3):
    return (N1*(1-N1)-((a1*N1*N2)/(1+b1*N1)))
def dN2dt(t,N1,N2,N3):
    return (((a1*N1*N2)/(1+b1*N1))-((a2*N2*N3)/(1+alpha*zai+b2*N2))-d1*N2)
def dN3dt(t,N1,N2,N3):
    return (((a2*(N2+c*zai)*N3)/(1+(alpha*zai)+b2*N2))-d2*N3)

def rk(dt,tend):
    t=np.arange(0,tend+dt,dt)
    N1=np.zeros_like(t)
    N2=np.zeros_like(t)
    N3=np.zeros_like(t)
    N1[0]=N1init
    N2[0]=N2init
    N3[0]=N3init
    for i in range(0,t.shape[0]-1):
        h1=dt*dN1dt(t[i],N1[i],N2[i],N3[i])
        j1=dt*dN2dt(t[i],N1[i],N2[i],N3[i])
        k1=dt*dN3dt(t[i],N1[i],N2[i],N3[i])
        h2=dt*dN1dt(t[i]+dt/2,N1[i]+h1/2,N2[i]+j1/2,N3[i]+k1/2)
        j2=dt*dN2dt(t[i]+dt/2,N1[i]+h1/2,N2[i]+j1/2,N3[i]+k1/2)
        k2=dt*dN3dt(t[i]+dt/2,N1[i]+h1/2,N2[i]+j1/2,N3[i]+k1/2)
        h3=dt*dN1dt(t[i]+dt/2,N1[i]+h2/2,N2[i]+j2/2,N3[i]+k2/2)
        j3=dt*dN2dt(t[i]+dt/2,N1[i]+h2/2,N2[i]+j2/2,N3[i]+k2/2)
        k3=dt*dN3dt(t[i]+dt/2,N1[i]+h2/2,N2[i]+j2/2,N3[i]+k2/2)
        h4=dt*dN1dt(t[i]+dt,N1[i]+h3,N2[i]+j3,N3[i]+k3)
        j4=dt*dN2dt(t[i]+dt,N1[i]+h3,N2[i]+j3,N3[i]+k3)
        k4=dt*dN3dt(t[i]+dt,N1[i]+h3,N2[i]+j3,N3[i]+k3)
        N1[i+1]=N1[i]+(h1+2.0*h2+2.0*h3+h4)/6.0
        N2[i+1]=N2[i]+(j1+2.0*j2+2.0*j3+j4)/6.0
        N3[i+1]=N3[i]+(k1+2.0*k2+2.0*k3+k4)/6.0
    return N1,N2,N3,t

N1,N2,N3,t=rk(0.01,100)
fig=plt.figure(figsize=(20,6))
plt.plot(t,N1,label='N1')
plt.plot(t,N2,label='N2')
plt.plot(t,N3,label='N3')
plt.title("Figure 1")


fig2 = plt.figure()
ax = Axes3D(fig2)
ax.plot(N1, N2, N3, 'b', lw=0.5)


ax.plot(N1, N2, N3, lw=0.5)
ax.set_xlabel('N2')
ax.set_ylabel('N1')
ax.set_zlabel('N3')
plt.tick_params(labelsize=15)
ax.set_title('Food_Chain', fontsize=15)

plt.show()


fig1, (ax1) = plt.subplots(1,1, figsize = (6, 6))
ax1.plot(N1, N2)
ax1.set_title("N1 & N2")
plt.show()

fig3, (ax1) = plt.subplots(1,1, figsize = (6,6))
ax1.plot(N1, N3)
ax1.set_title("N1 & N3")
plt.show()

fig3, (ax1) = plt.subplots(1,1, figsize = (6,6))
ax1.plot(N2, N3)
ax1.set_title("N2 & N3")
plt.show()

plt.plot(t, N1, t, N2, t, N3)
plt.show()