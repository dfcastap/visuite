# -*- coding: utf-8 -*-
"""
Created on Sat May 14 17:10:46 2016

@author: castaned
"""
import numpy as np
import main_modules as mmod

def xi(i,r1,theta,dr,phi,tstp):
    global dth,dr_arr,alambda
    dx1 = np.sin(np.deg2rad(i))
    #    dy1 = 0.
    dz1 = np.cos(np.deg2rad(i))
    x1 = r1*np.sin(np.deg2rad(theta))*np.cos(np.deg2rad(phi))
    y1 = r1*np.sin(np.deg2rad(theta))*np.sin(np.deg2rad(phi))
    z1 = r1*np.cos(np.deg2rad(theta))
    pert = 0.1*r1
    alambda = dr/(r1*np.deg2rad(tstp))
    r3 = pert*np.tan(alambda)
    r2 = np.sqrt((r1+pert)**2+r3**2)
    psi = np.rad2deg(np.arcsin(r3/r2))
    th2 = theta-psi
    x2 = r2*np.sin(np.deg2rad(th2))*np.cos(np.deg2rad(phi))
    y2 = r2*np.sin(np.deg2rad(th2))*np.sin(np.deg2rad(phi))    
    z2 = r2*np.cos(np.deg2rad(th2))
    dx2 = x2-x1
    dy2 = y2-y1
    dz2 = z2-z1
    v2 = np.sqrt(dx2**2+dy2**2+dz2**2)
    cosang = (dx1*dx2+dz1*dz2)/v2
    return cosang
    
def modelread(model_name,dr_arr,par):
    model = np.genfromtxt(model_name)
    if len(model[:,0])==10:
        temp1 = model[::-1]
        temp2 =  np.tile(0.,(20,5))
        for i in range(20):
            if i<10:
                temp2[i] = model[i]
            else:
                temp2[i,0] = i*(9.)+4.5
                temp2[i,1:5] = temp1[i-10,1:5]
        model=1.*temp2
    
    if par=="ODD":
        model = np.genfromtxt(model_name)
    
    dr_arr[0] = 0.
    dr_arr[1:20] = (np.diff(model[:,1]))
    dr_arr[20] = 0.
    #interp = modelinterpolate(model,dr_arr,nzth,nzphi)
    return model
    
def  modelinterpolate(model,dr_arr,nzth,nzphi):
    xi = model[:,0]
    #xi_dr = np.linspace(13.5,166.5,len(dr_arr))
    xi_dr = np.zeros(len(dr_arr))
    xi_dr[0] = 0.
    xi_dr[1:-1] = np.linspace(9,171,len(dr_arr)-2)
    xi_dr[-1] = 180.
    
    y = []
    x1 = np.linspace(0,180,nzth+1)
    midx = x1[0:nzth]+np.diff(x1)/2.
    midx = np.linspace(0,180,nzth)
    #phi = np.linspace(0,360,nzphi)
    # spline order: 1 linear, 2 quadratic, 3 cubic ... 
    #order = 1
    y.append(midx)
    # do inter/extrapolation
    for i in range(1,len(model[0,:])):
        #s = InterpolatedUnivariateSpline(xi, model[:,i], k=order)
        s = np.interp(midx,xi,model[:,i])
        #y.append(s(midx))
        y.append(s)
        
    #s = InterpolatedUnivariateSpline(xi_dr,dr_arr, k=order)
    y.append(np.interp(midx,xi_dr,dr_arr))
    y = np.array(y)
    y = np.transpose(y)
    return y

def find_cosxi(model_name,incl,par,fine_model=False):
    ### Stellar grid definition #####
    nzphi = 400	# Phi zones
    nzth = 200	# Theta zones
    dth = 180./nzth
    dphi = 360./nzphi


#    dr_arr = np.empty(21,dtype=float)
    phi = np.linspace(0,360.-360./nzphi,nzphi)
#    model = modelread(model_name,dr_arr,par)
#    interpmodel = modelinterpolate(model,dr_arr,nzth,nzphi)
#    stp = abs(model[0,0]-model[1,0])
    
    if fine_model==True:
        interpmodel_fine = np.genfromtxt(model_name)
        interpmodel = np.zeros((nzth,len(interpmodel_fine[0,:])+1))
        stp = interpmodel_fine[1,0]-interpmodel_fine[0,0]
        
        midp = np.empty(201)
        midp[0] = 0.
        midp[1:-1] = np.arange(stp/2., interpmodel_fine[-1,0], stp)
        midp[-1] = 190.
        dr_fine = np.empty(201)
        dr_fine[0] = 0. 
        dr_fine[1:-1] = np.diff(interpmodel_fine[:,1])
        dr_fine[-1] = 0.
        dr_fine = np.interp(interpmodel_fine[:,0],midp,dr_fine)
              
        
        interpmodel[:,0:-1] = interpmodel_fine
        interpmodel[:,-1] = dr_fine
        
        
        
    cossquiggle = []
    darea = []
    #ROTORC model dtheta:
    model_dth = stp
    
    #Geometric correction factor:
    correct = np.sqrt(1.+interpmodel[:,-1]**2/(interpmodel[:,1]*np.deg2rad(model_dth))**2)

    #Cossquiggle and darea calculation for every point in the grid:
    for angle in phi:
        xitemp = xi(incl[0],interpmodel[:,1],interpmodel[:,0],interpmodel[:,-1],angle,model_dth)

        cossquiggle.append(xitemp)
        radius = interpmodel[:,1]*6.9598e10
        darea.append(correct*np.sin(np.deg2rad(interpmodel[:,0]))*dth*dphi*(np.pi*radius)**2/(180.**2))
        
        
    
    #Convert lists to numpy.array:    
    cossquiggle = np.array(cossquiggle)
    
    return cossquiggle

#homedir = "/home/castaned/Documents/"
#static_m = homedir+"ROTORCmodels/visibilities/"
#
#rzone = 490
#
#vv = 0
#mass = "1p875"
#mde = "0 p5"
#m = mmod.emode(mass,vv,mde)
#where =static_m+mass+'Msun/V'+m.eq_vel+"/MODE_"+m.parity+"_"+str(m.index)+"/"
#modeln = "fine_model_MODE_"+str(m.index)+"_r"+str(rzone)
#incl = [0]
#
#cosxi = find_cosxi(where+modeln,incl,m.parity,fine_model=True)
