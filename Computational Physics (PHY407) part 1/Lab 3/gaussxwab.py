#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from gaussxw import gaussxw
def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

