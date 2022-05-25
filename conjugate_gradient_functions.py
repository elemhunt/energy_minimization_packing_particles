# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:38:22 2015

@author: javierseguradoescudero
"""
import numpy as np

   
def steepest_descent(dfun,fun,x0,tol,k):
    x=np.copy(x0)
    res=-1.*dfun(x,k)
    res_scalar=np.linalg.norm(res)
    iter=0.
    alpha=1
    while(res_scalar>tol):
         iter=iter+1
         p=res
         V_old=fun(x,k)
         alpha=line_search_golden(x,p,V_old,fun,alpha,k)
         x = x+ alpha*p
         res=-1.*dfun(x,k)
         res_scalar=np.linalg.norm(res)
         print ('Steepest, iter, energy,res',iter,fun(x,k),res_scalar)
    return x
    
    
def nonlinear_conjugate_gradient(dfun,fun,x0,tol,k):
    x=np.copy(x0)
    res=-1.*dfun(x,k) 
    p=np.copy(res)
    V=fun(x,k) 
    res_scalar=np.linalg.norm(res)
    iter=0
    alpha=1
    while(res_scalar>tol):

        iter=iter+1
        p_old=np.copy(p)      
        res_old=np.copy(res)
        V_old=V
        alpha=line_search_golden(x,p,V_old,fun,alpha,k)
        x = x+alpha*p_old    
        res=-1.*dfun(x,k) 
        V=fun(x,k) 
        res_scalar=np.linalg.norm(res)
        
        print ('NLCG, iter,energy,res',iter,V,res_scalar)
        if(res_scalar<tol ):
            break
# Several options to choose beta:        
#        beta=np.dot(res,res)/np.dot(res_old,res_old)   
        beta=(np.dot(res,res)-np.dot(res,res_old))/np.dot(res_old,res_old) 
        p=res+beta*p_old  
    return x


def line_search_golden(x,p,Vold,fun,alphaold,k):
    N=len(x)/2.
    r=(np.sqrt(5)-1.)/2.
    p_norm=np.max(abs(p))
    alpha_2=(1./p_norm)*N/10.
    alpha_2=min(alpha_2,1E12)
    alpha_1=max(alphaold*.01,alpha_2*1E-6)
    
    alpha_min=alpha_1
    Vmin=Vold
    alpha=alpha_1
    its=0
    
    while (alpha<alpha_2): # increasing alpha
        its+=1
        V=fun(x+alpha*p,k)
        if(V<Vmin):
            alpha_min=alpha
        alpha/=r    
#    if abs(alpha_min-alpha_1)/alpha_1<1E-10: # decrasing alpha
#        print('alpha_1 original',alpha_1) 
#        alpha_2=alpha_1
#        alpha_1=alpha_1/1000.
#        alpha=alpha_1
#        r=0.9
#        while (alpha<alpha_2):
#            its+=1
#            V=fun(x+alpha*p,k)
#            if(V<Vmin):
#                alpha_min=alpha
#            alpha/=r   
#            print ('new alphas in LS',alpha_min,V,Vold,Vmin)
    return alpha_min

  