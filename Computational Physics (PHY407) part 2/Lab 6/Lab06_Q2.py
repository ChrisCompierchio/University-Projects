#!/usr/bin/env python
# coding: utf-8

# Sky Kapoor and Chris Compierchio
# 
# This program will simulate the motion of an arbitrary number of particles under the influence of the Leonard-Jones potential

# In[1]:


#import libraries
import numpy as np
import matplotlib.pyplot as plt


# PART A

# In[2]:


#intitialize N, dt, T, U, K
N = 16
dt = 0.01
T = np.arange(0, .2, dt)
U = np.zeros(len(T))

#create arrays for the xs, ys, rs, velocities and accelerations
xs = [[0]*N for i in range(T.size)]
ys = [[0]*N for i in range(T.size)]
rxs = [[0]*N for i in range(T.size)]
rys = [[0]*N for i in range(T.size)]
vxs = [[0]*N for i in range(T.size)]
vys = [[0]*N for i in range(T.size)]
axs = [[0]*N for i in range(T.size)]
ays = [[0]*N for i in range(T.size)]

#Define constants
m = 1
epsilon = 1
sigma = 1

#set initial positions of particles
for i in range(N):
    if i == 0 or i == 4 or i == 8 or i == 12:
        xs[0][i] = 0.5
    if i == 1 or i == 5 or i == 9 or i == 13:
        xs[0][i] = 1.5
    if i == 2 or i == 6 or i == 10 or i == 14:
        xs[0][i] = 2.5
    if i == 3 or i == 7 or i == 11 or i == 15:
        xs[0][i] = 3.5
for j in range(N):
    if j == 0 or j == 1 or j == 2 or j == 3:
        ys[0][j] = 0.5
    if j == 4 or j == 5 or j == 6 or j == 7:
        ys[0][j] = 1.5
    if j == 8 or j == 9 or j == 10 or j == 11:
        ys[0][j] = 2.5
    if j == 12 or j == 13 or j == 14 or j == 15:
        ys[0][j] = 3.5

#plot initial positions
plt.scatter(xs[:][0], ys[:][0], label = "particle", color="grey")
plt.xlim(-10,16)
plt.ylim(-10,16)

#define a function that will calculate the acceleration of the particles
def acc(x, y):
    u=0
    for i in range(N):
        for j in range(N):
            if i != j:
                x = xs[0][i] - xs[0][j]
                y = ys[0][i] - ys[0][j]
                axs[0][i] += 4*((12/(np.sqrt(x**2 + y**2))**13)-(6/(np.sqrt(x**2 + y**2))**7))*x/np.sqrt(x**2 + y**2)
                ays[0][i] += 4*((12/(np.sqrt(x**2 + y**2))**13)-(6/(np.sqrt(x**2 + y**2))**7))*y/np.sqrt(x**2 + y**2)
    return axs, ays, u
    
#loop through the particles
for i in range (T.size-1):
    for j in range(N):
        
        #For particle 1, update its x, y, accelerations, and velocities
        xs[i+1][0] = xs[i][0] + vxs[i][0]*dt + acc(xs, ys)[0][i][0]*.5*dt**2
        ys[i+1][0] = ys[i][0] + vys[i][0]*dt + acc(xs, ys)[1][i][0]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][0]
        ay = acc(xs[i], ys[i])[1][i+1][0]
        
        vxs[i+1][0] = vxs[i][0] + 0.5*(acc(xs[i],ys[i])[0][i][0]+ax)*dt
        vys[i+1][0] = vys[i][0] + 0.5*(acc(xs[i],ys[i])[1][i][0]+ay)*dt
        
        #repeat for particles 2-16
        xs[i+1][1] = xs[i][1] + vxs[i][1]*dt + acc(xs, ys)[0][i][1]*.5*dt**2
        ys[i+1][1] = ys[i][1] + vys[i][1]*dt + acc(xs, ys)[1][i][1]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][1]
        ay = acc(xs[i], ys[i])[1][i+1][1]
        
        vxs[i+1][1] = vxs[i][1] + 0.5*(acc(xs[i],ys[i])[0][i][1]+ax)*dt
        vys[i+1][1] = vys[i][1] + 0.5*(acc(xs[i],ys[i])[1][i][1]+ay)*dt
        
        
        xs[i+1][2] = xs[i][2] + vxs[i][2]*dt + acc(xs, ys)[0][i][2]*.5*dt**2
        ys[i+1][2] = ys[i][2] + vys[i][2]*dt + acc(xs, ys)[1][i][2]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][2]
        ay = acc(xs[i], ys[i])[1][i+1][2]
        
        
        vxs[i+1][2] = vxs[i][2] + 0.5*(acc(xs[i],ys[i])[0][i][2]+ax)*dt
        vys[i+1][2] = vys[i][2] + 0.5*(acc(xs[i],ys[i])[1][i][2]+ay)*dt
        
        
        xs[i+1][3] = xs[i][3] + vxs[i][3]*dt + acc(xs, ys)[0][i][3]*.5*dt**2
        ys[i+1][3] = ys[i][3] + vys[i][3]*dt + acc(xs, ys)[1][i][3]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][3]
        ay = acc(xs[i], ys[i])[1][i+1][3]
        
        vxs[i+1][3] = vxs[i][3] + 0.5*(acc(xs[i],ys[i])[0][i][3]+ax)*dt
        vys[i+1][3] = vys[i][3] + 0.5*(acc(xs[i],ys[i])[1][i][3]+ay)*dt
       
    
        xs[i+1][4] = xs[i][4] + vxs[i][4]*dt + acc(xs, ys)[0][i][4]*.5*dt**2
        ys[i+1][4] = ys[i][4] + vys[i][4]*dt + acc(xs, ys)[1][i][4]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][4]
        ay = acc(xs[i], ys[i])[1][i+1][4]
        
        
        vxs[i+1][4] = vxs[i][4] + 0.5*(acc(xs[i],ys[i])[0][i][4]+ax)*dt
        vys[i+1][4] = vys[i][4] + 0.5*(acc(xs[i],ys[i])[1][i][4]+ay)*dt
        
        
        xs[i+1][5] = xs[i][5] + vxs[i][5]*dt + acc(xs, ys)[0][i][5]*.5*dt**2
        ys[i+1][5] = ys[i][5] + vys[i][5]*dt + acc(xs, ys)[1][i][5]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][5]
        ay = acc(xs[i], ys[i])[1][i+1][5]
        
        
        vxs[i+1][5] = vxs[i][5] + 0.5*(acc(xs[i],ys[i])[0][i][5]+ax)*dt
        vys[i+1][5] = vys[i][5] + 0.5*(acc(xs[i],ys[i])[1][i][5]+ay)*dt
        
        
        xs[i+1][6] = xs[i][6] + vxs[i][6]*dt + acc(xs, ys)[0][i][6]*.5*dt**2
        ys[i+1][6] = ys[i][6] + vys[i][6]*dt + acc(xs, ys)[1][i][6]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][6]
        ay = acc(xs[i], ys[i])[1][i+1][6]
       
        
        vxs[i+1][6] = vxs[i][6] + 0.5*(acc(xs[i],ys[i])[0][i][6]+ax)*dt
        vys[i+1][6] = vys[i][6] + 0.5*(acc(xs[i],ys[i])[1][i][6]+ay)*dt
        
        
        xs[i+1][7] = xs[i][7] + vxs[i][7]*dt + acc(xs, ys)[0][i][7]*.5*dt**2
        ys[i+1][7] = ys[i][7] + vys[i][7]*dt + acc(xs, ys)[1][i][7]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][7]
        ay = acc(xs[i], ys[i])[1][i+1][7]
        
        
        vxs[i+1][7] = vxs[i][7] + 0.5*(acc(xs[i],ys[i])[0][i][7]+ax)*dt
        vys[i+1][7] = vys[i][7] + 0.5*(acc(xs[i],ys[i])[1][i][7]+ay)*dt
        
        
        xs[i+1][8] = xs[i][8] + vxs[i][8]*dt + acc(xs, ys)[0][i][8]*.5*dt**2
        ys[i+1][8] = ys[i][8] + vys[i][8]*dt + acc(xs, ys)[1][i][8]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][8]
        ay = acc(xs[i], ys[i])[1][i+1][8]
       
        
        vxs[i+1][8] = vxs[i][8] + 0.5*(acc(xs[i],ys[i])[0][i][8]+ax)*dt
        vys[i+1][8] = vys[i][8] + 0.5*(acc(xs[i],ys[i])[1][i][8]+ay)*dt
        
        
        xs[i+1][9] = xs[i][9] + vxs[i][9]*dt + acc(xs, ys)[0][i][9]*.5*dt**2
        ys[i+1][9] = ys[i][9] + vys[i][9]*dt + acc(xs, ys)[1][i][9]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][9]
        ay = acc(xs[i], ys[i])[1][i+1][9]
        
        
        vxs[i+1][9] = vxs[i][9] + 0.5*(acc(xs[i],ys[i])[0][i][9]+ax)*dt
        vys[i+1][9] = vys[i][9] + 0.5*(acc(xs[i],ys[i])[1][i][9]+ay)*dt
        
        
        xs[i+1][10] = xs[i][10] + vxs[i][10]*dt + acc(xs, ys)[0][i][10]*.5*dt**2
        ys[i+1][10] = ys[i][10] + vys[i][10]*dt + acc(xs, ys)[1][i][10]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][10]
        ay = acc(xs[i], ys[i])[1][i+1][10]
        
        
        vxs[i+1][10] = vxs[i][10] + 0.5*(acc(xs[i],ys[i])[0][i][10]+ax)*dt
        vys[i+1][10] = vys[i][10] + 0.5*(acc(xs[i],ys[i])[1][i][10]+ay)*dt
        
        
        xs[i+1][11] = xs[i][11] + vxs[i][11]*dt + acc(xs, ys)[0][i][11]*.5*dt**2
        ys[i+1][11] = ys[i][11] + vys[i][11]*dt + acc(xs, ys)[1][i][11]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][11]
        ay = acc(xs[i], ys[i])[1][i+1][11]
        
        vxs[i+1][11] = vxs[i][11] + 0.5*(acc(xs[i],ys[i])[0][i][11]+ax)*dt
        vys[i+1][11] = vys[i][11] + 0.5*(acc(xs[i],ys[i])[1][i][11]+ay)*dt
        
        
        xs[i+1][12] = xs[i][12] + vxs[i][12]*dt + acc(xs, ys)[0][i][12]*.5*dt**2
        ys[i+1][12] = ys[i][12] + vys[i][12]*dt + acc(xs, ys)[1][i][12]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][12]
        ay = acc(xs[i], ys[i])[1][i+1][12]
        
        vxs[i+1][12] = vxs[i][12] + 0.5*(acc(xs[i],ys[i])[0][i][12]+ax)*dt
        vys[i+1][12] = vys[i][12] + 0.5*(acc(xs[i],ys[i])[1][i][12]+ay)*dt
        
        
        xs[i+1][13] = xs[i][13] + vxs[i][13]*dt + acc(xs, ys)[0][i][13]*.5*dt**2
        ys[i+1][13] = ys[i][13] + vys[i][13]*dt + acc(xs, ys)[1][i][13]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][13]
        ay = acc(xs[i], ys[i])[1][i+1][13]
        
        vxs[i+1][13] = vxs[i][13] + 0.5*(acc(xs[i],ys[i])[0][i][13]+ax)*dt
        vys[i+1][13] = vys[i][13] + 0.5*(acc(xs[i],ys[i])[1][i][13]+ay)*dt
        
        
        xs[i+1][14] = xs[i][14] + vxs[i][14]*dt + acc(xs, ys)[0][i][14]*.5*dt**2
        ys[i+1][14] = ys[i][14] + vys[i][14]*dt + acc(xs, ys)[1][i][14]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][14]
        ay = acc(xs[i], ys[i])[1][i+1][14]
        
        vxs[i+1][14] = vxs[i][14] + 0.5*(acc(xs[i],ys[i])[0][i][14]+ax)*dt
        vys[i+1][14] = vys[i][14] + 0.5*(acc(xs[i],ys[i])[1][i][14]+ay)*dt
        
        
        xs[i+1][15] = xs[i][15] + vxs[i][15]*dt + acc(xs, ys)[0][i][15]*.5*dt**2
        ys[i+1][15] = ys[i][15] + vys[i][15]*dt + acc(xs, ys)[1][i][15]*.5*dt**2
        
        ax = acc(xs[i], ys[i])[0][i+1][15]
        ay = acc(xs[i], ys[i])[1][i+1][15]
        
        vxs[i+1][15] = vxs[i][15] + 0.5*(acc(xs[i],ys[i])[0][i][15]+ax)*dt
        vys[i+1][15] = vys[i][15] + 0.5*(acc(xs[i],ys[i])[1][i][15]+ay)*dt
    
#plot trajectories
plt.plot(xs, ys, ".")
plt.title("Particle Trajectories", fontsize = 16)
plt.xlabel("X", fontsize = 14)
plt.ylabel("Y", fontsize = 14)
plt.show()


# PART B

# In[3]:


#Initialize Kinetic and potential energies
K=[[0]*N for i in range(T.size)]
U = [[0]*N for i in range(T.size)]
#Calculate Kinetic and potential energies
for i in range (1,len(T)):
    for j in range(N):
        K[i][j] += .5*vxs[i][j]**2+0.5*vys[i][j]**2
        U[i][j] += 0.5*(4*((1/(np.sqrt(xs[i][j]**2 + ys[i][j]**2))**12)-(1/(np.sqrt(xs[i][j]**2 + ys[i][j]**2))**6)))




# In[67]:


#for i in range(N):
    #print(len(T[:16]), len(np.add(U, K)[1][i]))
for i in range (2,N):    
    plt.plot(T[i], np.add(U, K)[i][2], ".")
    plt.ylim(10000, 1.3e4)
    plt.title("Total Energies for Each Particle", fontsize = 16)
    plt.ylabel("Energy (J)", fontsize = 14)
    plt.xlabel("Time(s)", fontsize = 14)


# PART C

# In[ ]:




