#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 11:39:07 2018

@author: djhwueng
"""


# https://en.wikipedia.org/wiki/Euler–Maruyama_method

"""The following Python code implements the Euler-Maruyama method 
and uses it to solve the Ornstein-Uhlenbeck process"""

import numpy as np
import matplotlib.pyplot as plt
num_sims = 5 #display five runs

t_init = 3 
t_end = 7
N = 1000 # copute 1000 grid points
dt = float(t_end-t_init) / N
y_init = 0

c_theta = 0.7
c_mu = 1.5
c_sigma = 0.06

def mu(y,t):
    """implement the Ornstein-Uhlenbeck mu,""" 
    ### =\theta(\mu-Y_t)
    return c_theta*(c_mu-y)

def sigma(y,t):
    """implement the Ornstein-Uhlenbeck sigma""" 
    ### = \sigma
    return c_sigma


def sigma_cir(y,t):
    """implement the CIR sigma"""
    return c_sigma*np.sqrt(y)

def dW(delta_t):
    """sample a random number at each call"""
    return np.random.normal(loc =0.0,scale=np.sqrt(delta_t))

ts = np.arange(t_init,t_end,dt)
ys = np.zeros(N)

ys[0] = y_init

for _ in range(num_sims):
    for i in range(1,ts.size):
        t = (i-1)*dt
        y = ys[i-1]
        ys[i] = y + mu(y,t)*dt +  sigma_cir(y,t)*dW(dt)
    plt.plot(ts,ys)
plt.show()

    
    

