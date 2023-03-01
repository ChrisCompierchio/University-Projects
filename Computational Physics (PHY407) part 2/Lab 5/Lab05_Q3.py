#!/usr/bin/env python
# coding: utf-8

# Sky Kapoor and Chris Compierchio
# 
# This program will analyze a contour plot of an SLP anomaly using the fast Fourier tranform

# In[5]:


#import libraries
import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw
from scipy import constants as c


# PART A

# In[8]:


#Load data from text files
SLP = np.loadtxt('SLP.txt')
longitude = np.loadtxt('lon.txt')
times = np.loadtxt('times.txt')

#Plot contour of data
plt.contourf(longitude, times, SLP)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar()
plt.show()


# In[52]:


#Get wavenumbers
wavenumbers = np.concatenate((range(0,int(len(SLP)/2)),range(-int(len(SLP)/2),0) ))

#Plot contours for the wavenumbers m=3 and 5
plt.contourf(longitude, times[:3], SLP[wavenumbers[:3]])
plt.title("SLP anomaly (hPa) for m = 3", fontsize = 16)
plt.xlabel('longitude(degrees)', fontsize = 14)
plt.ylabel('days since Jan. 1 2015', fontsize = 14)
plt.colorbar()
plt.show()

plt.contourf(longitude, times[:5], SLP[wavenumbers[:5]])
plt.title("SLP anomaly (hPa) for m = 5", fontsize = 16)
plt.xlabel('longitude(degrees)', fontsize = 14)
plt.ylabel('days since Jan. 1 2015', fontsize = 14)
plt.colorbar()
plt.show()


# In[ ]:




