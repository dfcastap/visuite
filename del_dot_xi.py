# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 00:10:11 2015

@author: Diego
"""

import numpy as np

#import pylab as plt

import pandas as pd
#import pyLegendre_anybf as pyL
#import scipy.interpolate as interp
#from scipy.misc import derivative

par = "temp"
model = ["2p5"]
vel = 1 #index, not velocity!
# NRO Mode filename:
modefname="MODE_temp_18"
norm_f = True
scale = 0.1
depth = 10 #radial zones down from the surface

def calcdeldotxi(par,model,vel,modeloc,modefname):
    G = 6.67259e-8
    Msun = 1.99e33
    Rsun = 6.958e10

    PIGR = 4.*np.pi*G*((Rsun)**2)
    SPIGR = np.sqrt(PIGR)
    freq_unit = np.sqrt(4.*np.pi*G)
    
    global folder, rotorc_f,g3m1
    #--------------------------
    #M1p875
    m1vels = ["0","35","62","83","105","125","146","165","187","207"]
    #M2
    m2vels = ["0","36","63","84","106","127","148","168","190","211"]
    #M2p25
    m3vels = ["0","36","65","87","109","131","152","173","195","217"]
    #M2p5
    m4vels = ["0","37p5","67","89","111","134","156","178","200","222"]
    #M3
    m5vels = ["0"]
    if model[0] == "1p875": vels=m1vels
    if model[0] == "2p25": vels=m3vels
    if model[0] == "2": vels=m2vels
    if model[0] == "2p5": vels=m4vels
    if model[0] == "3": vels=m5vels
    #---------------------------
    
    
    folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+vels[vel]+"/"
    
    rotorc_f = "/home/diego/Documents/From_Bob/Delta_Scuti_2010/"+model[0]+"Msun/"+model[0]+"Msun_V"+vels[vel]+"/"
    
    modefname = modeloc+modefname
    #modefname = folder+modefname
    
    #Set up rotoc 'rs' angles:
    rs_ang = np.empty(12) #This initializes an array with 12 elements
    rs_ang[0] = 0
    rs_ang[1:11] = np.arange(4.5,90,9) #[4.5,...,85.5]
    rs_ang[11] = 90
    
    #NRO angles:
    nro_ang = np.linspace(10,80,8) #[10,20,30,40,50,60,70,80]
    
    #Open and read MODE file:
    f = open(modefname).read().split("\n")[26:] # This also skips the header of the MODE file
    tempout = [i.split() for i in f] # Get individual elements in each line
    
    #Get only the Z values:-------------------------
    nro_out=[]
    a=0
    for i in range(len(f)):
        b=i
        if len(f[i])==0:
            nro_out.append(np.array(tempout[a:b-1],dtype=float)) # Get the 8 patches of 'Z' values
            a = b+2
    #-----------------------------------------------
    
    
    # for each NRO angle nro_out[i][:,0-4], each columns are:
    # ZR -> [:,0] , ZT -> [:,1], ZP -> [:,2], ZG -> [:,3], r -> [:,4]
    # Next: define variables with names that make more sense
    ZR = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    ZT = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    ZP = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    ZG = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    r = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    
    # Extract and save relevant MODE quantities in arrays with appropriate names
    for i in range(len(nro_ang)):
        ZR[:,i] = nro_out[i][:,0]
        ZT[:,i] = nro_out[i][:,1]
        ZP[:,i] = nro_out[i][:,2]
        ZG[:,i] = nro_out[i][:,3]
        r[:,i] = nro_out[i][:,4]
    
    # read model properties:
    RS_read = np.genfromtxt(folder+"RS_Dmod_"+model[0]+"M_V"+vels[vel])
    VS_read = np.genfromtxt(folder+"VS_Dmod_"+model[0]+"M_V"+vels[vel])
    GR_read = np.genfromtxt(folder+"GR_Dmod_"+model[0]+"M_V"+vels[vel])
    GT_read = np.genfromtxt(folder+"GT_Dmod_"+model[0]+"M_V"+vels[vel])
    G3m1_read = np.genfromtxt(folder+"GAM3m1_Dmod_"+model[0]+"M_V"+vels[vel])
    Gamma1_read = np.genfromtxt(folder+"GAMMA1_Dmod_"+model[0]+"M_V"+vels[vel])
    
    #### Do the r interpolation:
    idx = np.empty((len(r),2)) # this will store the i and i+1 position to use in the radial interpolations
    wgt = np.empty(len(r)) # this will save the interpolation factor to use with all the rotorc variables
    j=0
    for i in range(len(r[:,0])):
        if i>0:
            j = idx[i-1,0]
        while j < len(RS_read[:,11])-1:
            if r[i,0]>=RS_read[j,11] and r[i,0]<=RS_read[j+1,11]:
                idx[i,:] = [j,j+1] # store the appropriate indexes
                wgt[i] = 1.-(r[i,0]-RS_read[j,11])/(RS_read[j+1,11]-RS_read[j,11]) # store the interpolation factor
                j = 1e99 #setting j to this value breaks the while loop
            j+=1
    # Interpolate model properties to NRO angles: (12) 'rs' angles -> (8) nro angles
    rs_nro_ang = np.empty((len(RS_read),len(nro_ang)))
    vs_nro_ang = np.empty((len(RS_read),len(nro_ang)))
    gr_nro_ang = np.empty((len(RS_read),len(nro_ang)))
    gt_nro_ang = np.empty((len(RS_read),len(nro_ang)))
    g3m1_nro_ang = np.empty((len(G3m1_read),len(nro_ang)))
    gamma1_nro_ang = np.empty((len(Gamma1_read),len(nro_ang)))
    
    for i in range(len(RS_read)):
        # The f(new x) = np.interp(new x, x, f(x)) function does a linear interpolation
        rs_nro_ang[i] = np.interp(nro_ang,rs_ang,RS_read[i])
    
        vs_nro_ang[i] = np.interp(nro_ang,rs_ang,VS_read[i])
    
        gr_nro_ang[i] = np.interp(nro_ang,rs_ang,GR_read[i])
    
        gt_nro_ang[i] = np.interp(nro_ang,rs_ang,GT_read[i])
        
        g3m1_nro_ang[i] = np.interp(nro_ang,rs_ang,G3m1_read[i])
        
        gamma1_nro_ang[i] = np.interp(nro_ang,rs_ang,Gamma1_read[i])
    
    
    # Interpolate from RS to NRO r:
    RS = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    VS = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    GR = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    GT = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    g3m1 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    gamma1 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    
    
    for i in range(len(nro_ang)):
        for j in range(len(RS[:,0])):
            # For each value of r, the appropriate RS, VS, GR, GT are used using the indexes found in line 77
            RS[j,i] = RS_read[idx[j,0],11]*wgt[j]+RS_read[idx[j,1],11]*(1.-wgt[j])
            VS[j,i] = VS_read[idx[j,0],11]*wgt[j]+VS_read[idx[j,1],11]*(1.-wgt[j])
            GR[j,i] = GR_read[idx[j,0],11]*wgt[j]+GR_read[idx[j,1],11]*(1.-wgt[j])
            GT[j,i] = GT_read[idx[j,0],11]*wgt[j]+GT_read[idx[j,1],11]*(1.-wgt[j])
            g3m1[j,i] = G3m1_read[idx[j,0],11]*wgt[j]+G3m1_read[idx[j,1],11]*(1.-wgt[j])
            gamma1[j,i] = Gamma1_read[idx[j,0],11]*wgt[j]+Gamma1_read[idx[j,1],11]*(1.-wgt[j])
    
    
    ################ calculate DEL-DOT-XI with EQ 10 (clement98)
    deldotxi10 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    # Calculate xi_r and xi_t:
    xi_r = ZR[:,:]*RS[:,:]
    # xi_t also needs to be multiplied by sin(theta)cos(theta)
    xi_t = np.transpose([ZT[:,i]*RS[:,i]*np.sin(np.deg2rad(nro_ang[i]))*np.cos(np.deg2rad(nro_ang[i])) for i in range(len(nro_ang))])
    
    for i in range(len(nro_ang)):
        # Calculation of xi dot g to be used in eq 10
        xi_dot_g = xi_r[:,i]*GR[:,i]+xi_t[:,i]*GT[:,i]
        # Calculation of deldotxi:
        deldotxi10[:,i] = (1./(-1.*VS[:,i]))*(ZP[:,i]+ZG[:,i]+xi_dot_g)
        
    #print (1./(-1.*VS[-1,-1])),ZP[-1,-1],ZG[-1,-1],xi_dot_g[-1]
    
    """
    ############### calculate DEL-DOT-XI with EQ 14 (clement98)
    deldotxi14 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    
    
    #dxi_t = np.zeros(RS.shape)
    
    # Define xi_t*sin(theta) for the eq14 calculation
    xi_tsint = np.transpose([xi_t[:,i]*(np.sin(np.deg2rad(nro_ang[i]))) for i in range(len(nro_ang))])
    # Also initialize the container of its derivative
    dxi_tsint = np.zeros(RS.shape)
    #### The next block is a very generic script I took from the NRO_viewer to get the legendre function
    #### from a set of points. In this case I use xi_t :------------------------------------------------
    for i in range(len(ZR[:,0])):
        # For sticking with the generic code I rename some variables:
        data1 = xi_tsint[-1,:]
        # container1 will hold the theta values and the legendre polynomial values every 0.9 degrees
        container1 = pyL.legendre(data1,8) 
        
        # I also use a generic function to characterize the legendre solution in order to evaluate
        # derivative at the NRO points. At the end it's just a linear interpolation but done in 1 line
        f = interp.interp1d(container1[:,0],container1[:,1])
    
        df = np.empty(len(nro_ang))
        for j in range(len(df)):
            # Another generic function that takes the legendre solution of xi_t and finds the slope at the
            # selected location with respect to theta
            df[j] = (derivative(f,nro_ang[j]))
        
        # Store the values on a matrix with a name that makes more sense
        dxi_tsint[i,:]=df[:]
    ############################################################--------------------------------
    
    
    
    # Also initialize the container of the derivative of r^2*xi_r:
    dr2xi_r = np.zeros(xi_r.shape)
    
    
    #### Similar code to calculate the derivative to the one used for xi_t on line 141
    #### This thime there's no need to find a Legendre solution.
    for i in range(len(nro_ang)):
        
        # Definition of the function that I will derive:
        f = interp.interp1d(RS[:,i],xi_r[:,i]*((RS[:,i])**2))
        df = np.zeros(len(RS[:,0]))
        for j in range(len(df)):
            try:
                # Find the derivative at the RS points:
                df[j] = derivative(f,RS[j,i],dx=1e-6)
            except:
                a=1
        
        # Store the derivatives
        dr2xi_r[:,i]=df[:]
    ############################################################--------------------------------
    
    # The r-component of EQ14 can be calculated now:
    r_comp2 = (1./RS[:,:]**2)*(dr2xi_r[:,:])
    
    # The theta-component calculation follows:
    for i in range(len(nro_ang)):
        theta = np.deg2rad(nro_ang[i])    
        t_comp = (1./(RS[:,i]*np.sin(theta)))*(dxi_tsint[:,i])    
        
        # Final calculation of deldotxi for EQ14:    
        deldotxi14[:,i]= r_comp2[:,i]+t_comp
        
    
    #np.savetxt(modefname+"_deldotxi_eq10",deldotxi10)
    #np.savetxt(modefname+"_deldotxi_eq14",deldotxi14)
    
    # Import gammas and other info (see below) for dT/T calculation. this comes from pulset_non_ad
    supportf = np.genfromtxt(rotorc_f+"visibility_file")
    supportf[:,1] = supportf[:,1]-2
    supportf_df = pd.DataFrame(supportf)
    supportf_df.columns = ['i','j','gamma1','g3m1','T','P','R','rho','g','v']
    g3m1_p = []
    gamma1_p = []
    r_p = []
    for i in range(10):
        g3m1_p.append(np.array(supportf_df[supportf_df['j']==i]['g3m1']))
        gamma1_p.append(np.array(supportf_df[supportf_df['j']==i]['gamma1']))
        r_p.append(np.array(supportf_df[supportf_df['j']==i]['R']))
    g3m1_p = np.transpose(np.array(g3m1_p))
    gamma1_p = np.transpose(np.array(gamma1_p))
    r_p = np.transpose(np.array(r_p))
    #need to go from rotorc angles and radii to NRO's

    g3m1_prset = np.empty((len(RS[:,0]),len(g3m1_p)))
    gamma1_prset =  np.empty((len(RS[:,0]),len(g3m1_p)))
    #r_pulset = np.empty(RS.shape)
    for i in range(len(g3m1_p)):
        g3m1_prset[:,i] = np.interp(RS[:,i],r_p[i],g3m1_p[i])
        gamma1_prset[:,i] = np.interp(RS[:,i],r_p[i],gamma1_p[i])

    g3m1_pulset = np.empty(RS.shape)
    gamma1_pulset = np.empty(RS.shape)
    #r_p_pulset = np.empty(RS.shape)
    
    
    for i in range(len(g3m1_prset[:,0])):
        g3m1_pulset[i,:] = np.interp(nro_ang,rs_ang[1:-1],g3m1_prset[i])
        gamma1_pulset[i,:] = np.interp(nro_ang,rs_ang[1:-1],gamma1_prset[i])
        #r_p_nro_ang[i,:] = np.interp(nro_ang,rs_ang[1:-1],r_p[i])
        
    
    #--------------------
    """ 
        
    #dt_t = -(g3m1_pulset)*deldotxi10
    dt_t = -(g3m1)*deldotxi10
    
    
    return xi_r,xi_t,dt_t,ZG*PIGR/Rsun,RS

def norm_and_scale(xi_r,xi_t,dt_t,ZG,norm_f,scale,depth):
    global folder, rotorc_f
    
    xi_r_n = np.empty(xi_r[-(depth+1):-1,:].shape)
    xi_t_n = np.empty(xi_r[-(depth+1):-1,:].shape)
    dt_t_n = np.empty(xi_r[-(depth+1):-1,:].shape)
    ZG_n = np.empty(xi_r[-(depth+1):-1,:].shape)
    
    if norm_f:    
        for i in np.arange(-(depth),0,1):
            i_max=np.argmax(np.abs(xi_r[i,:]))
            xi_r_n[i,:] =  xi_r[i,:]/xi_r[i,i_max]
            xi_t_n[i,:] =  xi_t[i,:]/xi_r[i,i_max]
            dt_t_n[i,:] =  dt_t[i,:]/xi_r[i,i_max]
            ZG_n[i,:] =  ZG[i,:]/xi_r[i,i_max]
        
        xi_r_n *= scale
        xi_t_n *= scale
        dt_t_n *= scale
        ZG_n *= scale
    else:
        xi_r_n = xi_r * scale
        xi_t_n = xi_t * scale
        dt_t_n = dt_t * scale
        ZG_n = ZG * scale
    
    return xi_r_n,xi_t_n,dt_t_n,ZG_n
    
def to_rotorc(xi_r_n,xi_t_n,dt_t_n,ZG_n):
    nro_ang = np.linspace(10,80,8)
    rot_ang = np.linspace(4.5,85.5,10)
    
    xi_r_rot = np.empty((len(xi_r_n[:,0]),len(rot_ang)))
    xi_t_rot = np.empty((len(xi_r_n[:,0]),len(rot_ang)))
    dt_t_rot = np.empty((len(xi_r_n[:,0]),len(rot_ang)))
    ZG_rot = np.empty((len(xi_r_n[:,0]),len(rot_ang)))
    
    for i in range(len(xi_r_n[:,0])):
        #container_r = pyL.legendre(xi_r_n[-1,:],8)
        #container_t = pyL.legendre(xi_t_n[-1,:],8)
        #container_dt = pyL.legendre(dt_t_n[-1,:],8)
    
        xi_r_rot[i] = np.interp(rot_ang,nro_ang+xi_t_n[i,:],xi_r_n[i,:])
        xi_t_rot[i] = np.interp(rot_ang,nro_ang+xi_t_n[i,:],xi_t_n[i,:])
        dt_t_rot[i] = np.interp(rot_ang,nro_ang+xi_t_n[i,:],dt_t_n[i,:])
        ZG_rot[i] = np.interp(rot_ang,nro_ang+xi_t_n[i,:],ZG_n[i,:])
        #xi_r_rot[i] = np.interp(rot_ang,container_r[:,0],xi_r_n[i,:])
        #xi_t_rot[i] = np.interp(rot_ang,container_r[:,0],xi_t_n[i,:])
        #dt_t_rot[i] = np.interp(rot_ang,container_r[:,0],dt_t_n[i,:])
        
    return xi_r_rot,xi_t_rot,dt_t_rot,ZG_rot
    

#r,xi_r,xi_t,dt_t,r = calcdeldotxi(par,model,vel,"_",modefname)
        
#xi_r_n,xi_t_n,dt_t_n = norm_and_scale(xi_r,xi_t,dt_t,norm_f,scale,depth)

#xi_r_rot,xi_t_rot,dt_t_rot = to_rotorc(xi_r_n,xi_t_n,dt_t_n)

#container_t = pyL.legendre(dt_t[-1,:],8)
#container_r = pyL.legendre(xi_r[-1,:],8)
#bob_deldotxi_mode13 = np.genfromtxt("BOB_June8_2015/M2p5_V0_mode13_surf_perturbations",skip_header=2)
