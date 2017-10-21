#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:14:30 2017

@author: andrew

Plot results of a main run
"""

import numpy as np
import matplotlib.pyplot as plt



data = np.genfromtxt('data/dipole-dipole6.L_5beta_1.0.csv', delimiter=',')

X = data[:,:1]
Y = data[:,1:2]
U = data[:,2:3]
V = data[:,3:4]
W = data[:,4:5]

#================
# 2D VECTOR FIELD
#================
plt.figure()
ax = plt.gca()
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
plt.draw()
plt.show()




