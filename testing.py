# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:22:44 2021

@author: 20183816
"""
import numpy as np
import time
import random
import math

def dot(u: list, v: list):
    time0 = time.time()
    #a = u@v
    time1 = time.time()
    b = u.T.dot(v)
    time2= time.time()
    c = sum([i*j for (i,j) in zip(u,v)])
    time3 = time.time()
    d =  sum(map(lambda x, y: x*y, u, v))
    time4 = time.time()
    print(#'@: '+time1-time0,a, 
          'numpy: '+str(time2-time1), b, '\n'+
          'zip: '+str(time3-time2), c, '\n'+
          'lambda: '+str(time4-time3), d)
    return d

def sqrt(x):
    time0 = time.time()
    a = x**(1/2)
    time1 = time.time()
    b = np.sqrt(x)
    time2 = time.time()
    c = math.sqrt(x)
    time3 = time.time()
    print('python: '+str(time1-time0), a, '\n'+
          'numpy: '+str(time2-time1), b, '\n'+ 
          'math: '+str(time3-time2), c)

u = np.random.rand(1000000,1)
v = np.random.rand(1000000,1)
d = dot(u,v)

sqrt(14234523454)