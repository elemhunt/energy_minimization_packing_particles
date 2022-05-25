#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 12:36:52 2020

@author: elemhunt

This is a code to calculate and compare the minimiazation methods of steepest desecnt, 
non-linear conjugent gradient, and Broyden-Fletcher-Goldfarb-Shanno for a randomly created set of particles.
"""

import random
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Circle
from energies_blind import energy, d_energy
import conjugate_gradient_functions as CGF 
from scipy.optimize import minimize



R = 1   # radius of the circle
L = 20*R    # Length of the window
N = 10    # number of cirlces
K =1.   # penalty number
mg =9.8   # gravity constant
r = N  # multiplcation factor 
tol = N*0.01*mg   # calulated tolerance


#%% Creating necessary plots for each graph

fig, ax = plt.subplots()
ax.set_xlim([0 ,L])
ax.set_ylim([0 ,L])
ax.set_title('Steepest Descent Minimization')

fig, ax2 = plt.subplots()
ax2.set_xlim([0 ,L])
ax2.set_ylim([0 ,L])
ax2.set_title('Nonlinear Conjugate Gradient Minimization')
#
fig, ax3 = plt.subplots()
ax3.set_xlim([0 ,L])
ax3.set_ylim([0 ,L])
ax3.set_title('Scipy.opt-Broyden-Fletcher-Goldfarb-Shanno Minimization')

#%% Creating the inital guess for x and y of the purposed circles
num=[]
# List of N random numbers from 20. 2 numbers per cirlce. The x and the y
for i in range(2*N):
    num.append(random.uniform(1,L))

# Create an array of the random numbers 
X=np.array(num)


#%% Calculations

while (K<500):
    #Steepest Descent
    Min_SD = CGF.steepest_descent(d_energy,energy,X,tol,K)
    
    #Nonlinear Conjugate Gradient
    Min_NCG = CGF.nonlinear_conjugate_gradient(d_energy,energy,X,tol,K)
    
    #Broyden-Fletcher-Goldfarb-Shanno
    SciP = minimize(energy,X,K, method='BFGS', jac=d_energy,tol=tol,options={'disp': True})
    Min_Sci = SciP.x
    
    K=r*K
   
#%% Creating the centers to all the calculated circles
    
for q in range(N):

#list of Xvalues && list of Yvalues
    
    x1 = Min_SD[0:2*N:2]
    y1 = Min_SD[1:2*N:2]
    x2 = Min_NCG[0:2*N:2]
    y2 = Min_NCG[1:2*N:2]
    x3 = Min_Sci[0:2*N:2]
    y3 = Min_Sci[1:2*N:2]

# list of cirlces and array of center points
    
Circles_1 = []
Circles_2 = []
Circles_3 = []

centers_1 = np.column_stack((x1,y1))
centers_2 = np.column_stack((x2,y2))
centers_3 = np.column_stack((x3,y3))

#%% Creating the cirlces from the calculaed centers and Plotting circles
 
# Steepest Descent   
for l in centers_1:    
    circle1 = Circle(l, R, color='blue')
    Circles_1.append(circle1)                
for m in range(len(Circles_1)):
    ax.add_artist(Circles_1[m])

# NonLinear CG
for n in centers_2:    
    circle2 = Circle(n, R, color='red')
    Circles_2.append(circle2)               
for o in range(len(Circles_2)):
    ax2.add_artist(Circles_2[o])

# Broyden-Fletcher-Goldfarb-Shanno  
for p in centers_3:    
    circle3 = Circle(p, R, color='green')
    Circles_3.append(circle3)               
for q in range(len(Circles_3)):
    ax3.add_artist(Circles_3[q])
plt.show()