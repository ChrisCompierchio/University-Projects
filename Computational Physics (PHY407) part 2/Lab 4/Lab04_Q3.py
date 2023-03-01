#!/usr/bin/env python
# coding: utf-8

# Sky Kapoor and Chris Compierchio
# 
# This Program will complee textbook exercises 6.10, 6.11, and 6.13

# In[2]:


#Import Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as con


# PART A

# In[15]:


#Define x and c
x = 1
c = 2
#Calculate x for an accuracy of 10e-6
for i in range(100):
    x = 1 - np.exp(-c*x)
    print(x)
    if round(x, 6) == 0.796812:
        break
    


# In[9]:


#Repeat but plot for C = 0-3
x = 1
c = np.arange(0,3,0.01)
for i in range(100):
    x = 1 - np.exp(-c*x)
plt.plot(c,x)
plt.title("c as a function of x", fontsize = "16")
plt.xlabel("c", fontsize = "16")
plt.ylabel("x", fontsize = "16")
plt.show()


# PART B

# In[15]:


#Repeat but print the iterations it takes to get to an accuracy of 10e-6
x = 1
c = 2
for i in range(100):
    x = 1 - np.exp(-c*x)
    print(i)
    if round(x, 6) == 0.796812:
        break


# In[32]:


#Repeat using overrelaxation
x2 = 0.5
omega = 0.69
for j in range(10):
    x2 = (1 - np.exp(-c*x2))*(1 + omega) - omega*x2
    print(j,x2)


# PART C

# In[49]:


#Define f
def f(x):    
    return 5*np.exp(-x) + x - 5

#Define Constants
eps = 10e-6
x1, x2 = 1, 100
num = 0

#Perform binary search algorithm described in textbook page 264 
while np.abs(x1-x2) > eps:
    midpoint = .5*(x1+x2)

    if (f(midpoint) > 0 and f(x1) > 0) or (f(midpoint) < 0 and f(x1) < 0):
        x1 = midpoint
    else:
        x2 = midpoint
    
    num+=1
    
    x = 0.5*(x1+x2)
#print results 
print(x, num)


# In[52]:


#Calculate b
b = con.h*con.c/(con.k*x)

#Calculate and print the surface temperature of the sun
surface_temp = b/502e-9
surface_temp


# In[ ]:




