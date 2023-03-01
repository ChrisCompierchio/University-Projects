#!/usr/bin/env python
# coding: utf-8

# Sky Kapoor and Chris Compierchio
# 
# This program will calulate the wavefunction for the ground state and the first two excited states of an asymmetric well.

# In[2]:


#Import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as con


# PART B

# In[24]:


#Define given constants
L = 5e-10
a = 1.6022e-18
M = 9.1094e-31

#Define H_mn function
def Hmn(m, n):
    if m == n:
        H = .5*a + ((np.pi**2)*(con.hbar**2)*(m**2))/(2*M*(L**2))
        return H
    elif (m + n) % 2 != 0:
        H = -(8*a*m*n)/((np.pi**2)*(m**2 - n**2)**2)
        return H
    else:
        H = 0
        return H
    
#Define function to integrate sing simpsons rule as in Lab 2
def simp(a, b, psi):
    h = (b-a)/psi.size
   
    m = (psi[0] + psi[-1])

    sum1 = 0

    for k in range(1,100,2):
        sum1 += psi[k] 

    for k in range(2,100, 2):
        sum2 = psi[k]
    
    simp = (h/3)*(m+4*sum1+2*sum2)
    
    return simp


# PART C

# In[25]:


#Define maximum number of indices for the H matrix
mmax, nmax = 10, 10
#Define H
H = np.empty((10, 10))
#Calculate each entry of H
for m in range(1, mmax+1):
    for n in range(1, nmax+1):
        H[m-1, n-1] = Hmn(m, n)
#Calculate the eigenvalues of H
evs = np.linalg.eigvalsh(H)
#Print results in eV
print(evs*6.241509e18)


# PART D

# In[26]:


#Repeat for a 100x100 matrix
mmax, nmax = 100, 100
H = np.empty((100, 100))
for m in range(1, mmax+1):
    for n in range(1, nmax+1):
        H[m-1, n-1] = Hmn(m, n)
evs = np.linalg.eigvalsh(H)
print(evs*6.241509e18)


# PART E

# In[47]:


#Calcualte the eigenvalues and eugenvectors for H
evals, evecs = np.linalg.eigh(H)

#Define position points for the particle
x = np.linspace(0,L,100)

#Loop for n up to 3
for n in range(3):
    #Define psi
    psi = 0
    
    #Calculate psi using the eigenvectors of H
    psi += evecs[n][0]*np.sin((n+1)*np.pi*x/L)
    
    #Integrate psi
    A = simp(0, L, np.abs(psi)**2)

    #Normalize psi
    psi /= np.sqrt(A)

    #Plot results
    plt.plot(x, np.abs(psi)**2, label = "n = "+ str(n))
    plt.legend(loc = "upper left")
    plt.title("WaveFunctions for n = 0,1,2", fontsize = 16)
    plt.xlabel("x",fontsize = 16)
    plt.ylabel("Psi(x)",fontsize = 16)
plt.show()


# In[ ]:




