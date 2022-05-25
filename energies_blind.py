# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 13:26:41 2015

@author: javierseguradoescudero
"""

D=2.
R=D/2.
mg =9.8
L=10*D
# Note : k are the x and k+1 are the y
import numpy as np

def V_floor(x):
    ## complete this function
    if (x < R):
        return (0.5)*mg*(x-R)**2
    elif (R <=x):
        return 0
  
def dV_floor(x):
    ## complete this function
     if (x < R):
        return mg*(x-R)
     elif (R <= x):
        return 0
def V(r):
    ## complete this function

    if (r < 2*R):
        return (0.5)*mg*(r-2*R)**2
    elif (2*R <= r):
        return 0

def dV(r):
    ## complete this function
    penal = 0.0
    if (r < 2*R):
        return mg*(r-2*R)
    elif (2*R <= r):
        return penal

def V_wall(x):
    ## complete this function
    penal = 0.0
    if (x<R):
        return (0.5)*mg*(x-R)**2
    elif (R <=x<= (L-R)):
        return penal
    elif (x>(L-R)):
        return (0.5)*(mg)*(x+R-L)**2
def dV_wall(x):
    ## complete this function
     penal = 0.0
     if (x<R):
        return mg*(x-R)
     elif (R <=x<= (L-R)):
        return penal
     elif (x>(L-R)):
        return (mg)*(x+R-L)
def energy(x,penal):
    n=len(x)
    energy=0.
    for k in range(0,int(n/2)):
        xk=np.array([x[2*k],x[2*k+1]])
        for i in range(k+1,int(n/2)):
            xi=np.array([x[2*i],x[2*i+1]])
            xki=xi-xk
            rki=np.linalg.norm(xki)
            energy+=penal*V(rki)

        energy+=mg*xk[1]
        # add floor and wall contributions to energy
        energy+= penal*(V_wall(xk[0]))
        energy+= penal*(V_floor(xk[1]))
    return energy



def d_energy(x,penal):   
    import numpy as np
    n=len(x)
    denergy=np.zeros(n)
   
    for k in range(0, int(n/2)):
        xk=np.array([x[2*k],x[2*k+1]]) 
        for i in range(0,int(n/2)):  
            if (k!=i):
                xi=np.array([x[2*i],x[2*i+1]])                           
                xki=xi-xk                  
                rki=np.linalg.norm(xki)    
                ddV=dV(rki)
                if(ddV!=0):
                    ddV=ddV/rki                                            
                    denergy[2*k]+=-penal*ddV*xki[0]
                    denergy[2*k+1]+=-penal*ddV*xki[1]

        denergy[2*k+1]+=mg   
        # add floor and wall contributions to dericative of energy
        denergy[2*k+1]+= penal*(dV_floor(xk[1]))
        denergy[2*k]+= penal*(dV_wall(xk[0]))

    return denergy

    