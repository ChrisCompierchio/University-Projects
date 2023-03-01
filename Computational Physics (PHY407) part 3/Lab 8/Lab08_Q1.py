#!/usr/bin/env python
# coding: utf-8

# Sky Kapoor and Christopher Compierchio
# 
# This program will calculate the electrostatic potential at each grid point in an electrostatic capacitor.

# In[1]:


#import libraries
import numpy as np
import matplotlib.pyplot as plt


# PART A

# In[2]:


# Code based on example 9.1 from Newman (with adjustment as needed)

#Constants
M = 100 # 100 x 100 points for a 10cm x 10cm box
V1, V2 = 1.0, -1.0
precision = 1e-6

# Create arrays to hold potential values
phi = np.zeros([M+1, M+1], float)

# Using mm to better fit with the 100 x 100 pts
phi[20:80, 20] = V1 # First metal plate, 20mm from the left
phi[20:80, 80] = V2 # Second metal plate, 80mm from the left

phiprime = np.empty([M+1, M+1], float)


# In[3]:


# Main loop
delta = 1.0
while delta > precision:
    phiprime = phi.copy()
    # Calculate new values of the potential
    for i in range(1, M):
        for j in range(1, M):
            if i == 0 or i == M or j == 0 or j == M:
                phi[i, j] = phi[i, j]
            elif (j == 20 or j == 80) and (i >= 20 and i <= 80):
                phi[i, j] = phi[i, j]
            else:
                phi[i, j] = (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])/4
    # Calculte maximum difference from old values
    delta = np. max(np.abs(phi-phiprime))


# In[9]:


#plot results
plt.contour(phi)
plt.title("Contour Plot using Gauss-Siedel Method", fontsize = 16)
plt.xlabel('X', fontsize = 14)
plt.ylabel('Y', fontsize = 14)
plt.show()


# PART B

# In[5]:


#set omega to 0.1
omega = 0.1

# Main loop
delta = 1.0
while delta > precision:
    phiprime = phi.copy()
    # Calculate new values of the potential
    for i in range(1, M):
        for j in range(1, M):
            if i == 0 or i == M or j == 0 or j == M:
                phi[i, j] = phi[i, j]
            elif (j == 20 or j == 80) and (i >= 20 and i <= 80):
                phi[i, j] = phi[i, j]
            else:
                phi[i, j] = (1+omega)*(phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])/4 - omega*phi[i, j]
    # Calculte maximum difference from old values
    delta = np. max(np.abs(phi-phiprime))


# In[10]:


#plot results
plt.contour(phi)
plt.title("Contour Plot using Gauss-Siedel Method, $\omega$ = 0.1", fontsize = 16)
plt.xlabel('X', fontsize = 14)
plt.ylabel('Y', fontsize = 14)
plt.show()


# In[7]:


#set omega to 0.5
omega = 0.5

# Main loop
delta = 1.0
while delta > precision:
    phiprime = phi.copy()
    # Calculate new values of the potential
    for i in range(1, M):
        for j in range(1, M):
            if i == 0 or i == M or j == 0 or j == M:
                phi[i, j] = phi[i, j]
            elif (j == 20 or j == 80) and (i >= 20 and i <= 80):
                phi[i, j] = phi[i, j]
            else:
                phi[i, j] = (1+omega)*(phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])/4 - omega*phi[i, j]
    # Calculte maximum difference from old values
    delta = np. max(np.abs(phi-phiprime))


# In[11]:


#plot results
plt.contour(phi)
plt.title("Contour Plot using Gauss-Siedel Method, $\omega$ = 0.5", fontsize = 16)
plt.xlabel('X', fontsize =14)
plt.ylabel('Y', fontsize = 14)
plt.show()


# In[ ]:




