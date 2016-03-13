# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 19:08:06 2015

@author: diego
"""
import numpy as np
import pyLegendre_anybf as pyL

def leg_interp(var,nbfs,par):
    if par == "EVEN": parity = 0
    if par == "ODD": parity = 1
    if par == "OE": parity = 2
    global leg
    xarr = np.linspace(10,80,nbfs)
    #data = np.genfromtxt(self.models[0]+"/temp_surf"+str(key))
    interp_var = np.empty((len(var[:,0]),100))
    for i in range(len(var[:,0])):
        if parity==0:
            newdata = 1.*var[i,:]
            leg = pyL.legendre(newdata,nbfs)
        elif parity==1:
            newdata = np.cos(np.deg2rad(xarr))*var[i,:]
            leg = pyL.legendre_odd(newdata,nbfs)
        elif parity==2:
            newdata = var[i,:]
            leg = pyL.legendre_odd(newdata,nbfs)
        
        interp_var[i] = leg[:,1]
        
    return interp_var

"""
import seaborn as sns
sns.set(style="white",rc={"figure.figsize": (6, 6)})
new_r = leg_interp(r,8,"EVEN")
new_zp = leg_interp(xi_r,8,"EVEN")
levels = np.linspace(np.min(new_zp),np.max(new_zp), 40)
theta = np.linspace(0,np.deg2rad(90),100)
newt,n_r = np.meshgrid(theta,new_r[:,0])

CS = plt.contourf((new_r*np.cos(newt)), (new_r*np.sin(newt)),new_zp, 100, cmap=plt.cm.jet,vmax=np.max(new_zp), vmin=np.min(new_zp))
CSl = plt.contour((new_r*np.cos(newt)), (new_r*np.sin(newt)), new_zp, 20, colors="k")
#CS = plt.contourf((new_r*np.cos(newt)), (new_r*np.sin(newt)), new_zp, cmap=plt.cm.Spectral,levels=levels)
#plt.axes().set_aspect('equal')
plt.xlabel("Radius [R$_{\odot}$]")
plt.ylabel("Radius [R$_{\odot}$]")
"""