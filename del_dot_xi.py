# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 00:10:11 2015

@author: Diego
"""

import numpy as np
import pylab as plt
import scipy.misc as dtv
from scipy.special import sph_harm
import os

vis_path = os.getcwd()
global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/diego/Documents/"
else:
    homedir = "/home/castaned/"

#import pandas as pd

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
l = 3

def sphm(theta):
    global l
    tmp = sph_harm(0, l, 0, theta).real
    return (tmp/np.max(tmp))

def calcdeldotxi(par,model,vel,modeloc,modefname):
    G = 6.67259e-8
    Msun = 1.99e33
    Rsun = 6.958e10

    PIGR = 4.*np.pi*G*((Rsun)**2)
    SPIGR = np.sqrt(PIGR)
    freq_unit = np.sqrt(4.*np.pi*G)
    
    global folder, rotorc_f,g3m1,ZT,ZP
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
    
    
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+vels[vel]+"/"
    
    rotorc_f = homedir+"From_Bob/Delta_Scuti_2010/"+model[0]+"Msun/"+model[0]+"Msun_V"+vels[vel]+"/"
    
    modefname = modeloc+modefname
    #modefname = modefname #viscalc local test
    
    #Set up rotoc 'rs' angles:
    rs_ang = np.empty(12) #This initializes an array with 12 elements
    rs_ang[0] = 0
    rs_ang[1:11] = np.arange(4.5,90,9) #[4.5,...,85.5]
    rs_ang[11] = 90
    
    #NRO angles:
    nro_ang = np.linspace(10,80,8) #[10,20,30,40,50,60,70,80]
    
    #Open and read MODE file:
    f = open(modefname).read().split("\n")
    sigma = float(f[6].split()[0].strip())
    f = f[26:] # This also skips the header of the MODE file
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
    
    global zpcito, zgcito
    zpcito = 1.*ZP
    zgcito = 1.*ZG
    
    # read model properties:
    RS_read = np.genfromtxt(folder+"RS_Dmod_"+model[0]+"M_V"+vels[vel])
    VS_read = np.genfromtxt(folder+"VS_Dmod_"+model[0]+"M_V"+vels[vel])
    GR_read = np.genfromtxt(folder+"GR_Dmod_"+model[0]+"M_V"+vels[vel])
    GT_read = np.genfromtxt(folder+"GT_Dmod_"+model[0]+"M_V"+vels[vel])
    G3m1_read = np.genfromtxt(folder+"GAM3m1_Dmod_"+model[0]+"M_V"+vels[vel])
    Gamma1_read = np.genfromtxt(folder+"GAMMA1_Dmod_"+model[0]+"M_V"+vels[vel])
    
    #### Do the r interpolation:
    idx = np.zeros((len(r),2),dtype=int) # this will store the i and i+1 position to use in the radial interpolations
    wgt = np.empty(len(r)) # this will save the interpolation factor to use with all the rotorc variables
    j=0
    for i in range(len(r[:,0])):
        if i>0:
            j = idx[i-1,0]
        while j < len(RS_read[:,11])-1:
            if r[i,0]>=RS_read[j,11] and r[i,0]<=RS_read[j+1,11]:
                idx[i,:] = [j,j+1] # store the appropriate indexes
                wgt[i] = 1.-(r[i,0]-RS_read[j,11])/(RS_read[j+1,11]-RS_read[j,11]) # store the interpolation factor
                break
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
            RS[j,i] = rs_nro_ang[idx[j,0],i]*wgt[j]+rs_nro_ang[idx[j,1],i]*(1.-wgt[j])
            VS[j,i] = vs_nro_ang[idx[j,0],i]*wgt[j]+vs_nro_ang[idx[j,1],i]*(1.-wgt[j])
            GR[j,i] = gr_nro_ang[idx[j,0],i]*wgt[j]+gr_nro_ang[idx[j,1],i]*(1.-wgt[j])
            GT[j,i] = gt_nro_ang[idx[j,0],i]*wgt[j]+gt_nro_ang[idx[j,1],i]*(1.-wgt[j])
            g3m1[j,i] = g3m1_nro_ang[idx[j,0],i]*wgt[j]+g3m1_nro_ang[idx[j,1],i]*(1.-wgt[j])
            gamma1[j,i] = gamma1_nro_ang[idx[j,0],i]*wgt[j]+gamma1_nro_ang[idx[j,1],i]*(1.-wgt[j])
    
    
    ################ calculate DEL-DOT-XI with EQ 10 (clement98)
    deldotxi10 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    # Calculate xi_r and xi_t:
    xi_r = ZR[:,:]*RS[:,:]
    # xi_t also needs to be multiplied by sin(theta)cos(theta)
    xi_t = np.transpose([ZT[:,i]*RS[:,i]*np.sin(np.deg2rad(nro_ang[i]))*np.cos(np.deg2rad(nro_ang[i])) for i in range(len(nro_ang))])
    
    dPhi = 1.*ZG
    
    dP = 1.*ZP
    
    
    xi_r /= RS
    xi_t /= RS
    dPhi /= RS
    dP /= RS
    
    if par=="ODD":
        print "Odd mode..."
        rstar = RS[-1,-1]
        xi_r = rstar*np.transpose([ZR[:,i]*np.cos(np.deg2rad(nro_ang[i])) for i in range(len(nro_ang))])
        
        #xi_t = np.transpose([(ZT[:,i] - ZP[:,i]/(sigma*np.cos(np.deg2rad(nro_ang[i])))**2)*np.sin(np.deg2rad(nro_ang[i]))*np.cos(np.deg2rad(nro_ang[i]))**2 for i in range(len(nro_ang))])
        xi_t = np.transpose([(np.sin(np.deg2rad(nro_ang[i]))*(np.cos(np.deg2rad(nro_ang[i]))**2)) * ZT[:,i] * rstar - np.sin(np.deg2rad(nro_ang[i]))*(ZP[:,i]*RS[:,i]*rstar)/(sigma**2) for i in range(len(nro_ang))])
        #xi_t = np.transpose([ZT[:,i]*((np.sin(np.deg2rad(nro_ang[i]))**2))*RS[:,i] + (ZP[:,i]/(sigma**2))*(RS[:,i]*np.sin(np.deg2rad(nro_ang[i]))-1./(np.cos(np.deg2rad(nro_ang[i]))**2)) for i in range(len(nro_ang))])
        
        
        dPhi = np.transpose([ZG[:,i]*RS[:,i]*rstar*np.cos(np.deg2rad(nro_ang[i])) for i in range(len(nro_ang))])
        #dPhi = ZG*RS
        
        dP = np.transpose([ZP[:,i]*RS[:,i]*rstar*np.cos(np.deg2rad(nro_ang[i])) for i in range(len(nro_ang))])
        #dP = ZP*RS
    
    for i in range(len(nro_ang)):
        # Calculation of xi dot g to be used in eq 10
        xi_dot_g = xi_r[:,i]*GR[:,i]+xi_t[:,i]*GT[:,i]
        # Calculation of deldotxi:
        deldotxi10[:,i] = (1./(-1.*VS[:,i]))*(dP[:,i]+dPhi[:,i]+xi_dot_g)
        
    """
    #EQUATION 13 CALCULATIONS TO COMPARE WITH xi_t
    ####### eq 13:
    global xi_t_e13,sph_vals,dsph
    xi_t_e13 = np.empty(xi_t.shape)
    
    
    for i in range(len(ZR[:,0])):
        data1 = dP[i,:]
        container1 = pyL.legendre(data1,8)
        if par=="ODD":
            container1 = pyL.legendre_odd(data1,8)
        
        #plt.plot(container1[:,1])
        dtp1 = np.deg2rad(container1[1:-1,0])
        dtp0 = np.deg2rad(container1[0:-2,0])
        df = (container1[1:-1,1]-container1[0:-2,1])/(dtp1-dtp0)
        f_df = np.interp(np.deg2rad(nro_ang),np.deg2rad(container1[0:-2,0]),df)
        
        xi_t_e13[i,:] = (1./sigma**2)*(f_df/RS[i,:])
        #xi_t_e13[i,:] = f_df
        
        
    plt.plot(xi_t_e13[-1,:]/np.max(xi_t_e13[-1,:]))

    deldotxi10_e13 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    
    sph_vals = np.empty(nro_ang.shape)
    #for i in range(len(nro_ang)):
        #sph_vals[i] = pyL.newLeg(3,np.cos(np.deg2rad(nro_ang[i])))
    sph_vals = sphm(np.deg2rad(nro_ang))
    
    data1 = sph_vals
    container1 = pyL.legendre(data1,8)
    if par=="ODD":
        container1 = pyL.legendre_odd(data1,8)
    
    
    dtp1 = np.deg2rad(container1[1:-1,0])
    dtp0 = np.deg2rad(container1[0:-2,0])
    df = (container1[1:-1,1]-container1[0:-2,1])/(dtp1-dtp0)
    dsph = np.interp(np.deg2rad(nro_ang),np.deg2rad(container1[0:-2,0]),df)
    plt.plot(dsph/np.max(dsph))
    
    
    for i in range(len(nro_ang)):
        # Calculation of xi dot g to be used in eq 10
        xi_dot_g_e13 = xi_r[:,i]*GR[:,i]+xi_t_e13[:,i]*GT[:,i]
        # Calculation of deldotxi:
        deldotxi10_e13[:,i] = (1./(-1.*VS[:,i]))*(dP[:,i]+dPhi[:,i]+xi_dot_g_e13)
    
    """
    """
    ############### calculate DEL-DOT-XI with EQ 14 (clement98)
    import pyLegendre_anybf as pyL
    deldotxi14 = np.empty((len(nro_out[0][:,0]),len(nro_ang)))
    deldotxi142 = np.empty((len(nro_out[0][:,0])-2,len(nro_ang)))
    
    #dxi_t = np.zeros(RS.shape)
    
    # Define xi_t*sin(theta) for the eq14 calculation
    xi_tsint = np.transpose([xi_t[:,i]*(np.sin(np.deg2rad(nro_ang[i]))) for i in range(len(nro_ang))])
    # Also initialize the container of its derivative
    dxi_tsint = np.zeros(RS.shape)
    #### The next block is a very generic script I took from the NRO_viewer to get the legendre function
    #### from a set of points. In this case I use xi_t :------------------------------------------------
    global fine_ZT,fine_ZP,fine_ang,fine_ang_mid,fine_RS,fine_xi_t,d_fine_xi_t
    for i in range(len(ZR[:,0])):
        # For sticking with the generic code I rename some variables:
        data1 = xi_tsint[i,:]
        container1 = pyL.legendre(data1,8)
        # container1 will hold the theta values and the legendre polynomial values every 0.9 degrees
        
        if par=="EVEN":
            container1 = pyL.legendre(data1,8)
            fine_ZT = pyL.legendre(ZT[i,:],8)
            fine_ZP = pyL.legendre(ZP[i,:],8)
            fine_RS = pyL.legendre(RS[i,:],8)
        else:
            container1 = pyL.legendre_odd(data1,8)
            fine_ZT = pyL.legendre(ZT[i,:],8)
            fine_ZP = pyL.legendre(ZP[i,:],8)
            fine_RS = pyL.legendre(RS[i,:],8)
            
        
        fine_ang = np.deg2rad(fine_ZT[:,0])
        fine_ang_mid = 0.5*(fine_ang[0:-2]+fine_ang[1:-1])
        
        if par=="EVEN":
            fine_xi_t = np.array([fine_ZT[j,1]*fine_RS[j,1]*np.sin(fine_ang[j])*np.cos(fine_ang[j]) for j in range(len(fine_ang))])
        else:
            fine_xi_t = np.array([((np.sin(fine_ang[j])**2)*(np.cos(fine_ang[j])**2)) * fine_ZT[j,1] * rstar - (np.sin(fine_ang[j])**2)*(fine_ZP[j,1]*fine_RS[j,1]*rstar)/(sigma**2) for j in range(len(fine_ang))])
            
        d_fine_xi_t = (fine_xi_t[1:-1] - fine_xi_t[0:-2])/(fine_ang[1:-1] - fine_ang[0:-2])
        d_fine_xi_t = np.interp(np.deg2rad(nro_ang),fine_ang_mid,d_fine_xi_t)
        
        #if i == len(ZR[:,0]) - 1 :
            #plt.plot(np.rad2deg(fine_ang),fine_ZT[:,1])
            #plt.plot(nro_ang,ZT[i,:],"o")
        
        df = (container1[1:-1,1]-container1[0:-2,1])/(container1[1:-1,0]-container1[0:-2,0])
        f_df = (np.interp(nro_ang,container1[0:-2,0],df))
        # I also use a generic function to characterize the legendre solution in order to evaluate
        # derivative at the NRO points. At the end it's just a linear interpolation but done in 1 line
        #f = interp.interp1d(container1[:,0],container1[:,1])
    
        #df = np.empty(len(nro_ang))
        #for j in range(len(df)):
            # Another generic function that takes the legendre solution of xi_t and finds the slope at the
            # selected location with respect to theta
            #df[j] = (derivative(f,nro_ang[j]))
        
        # Store the values on a matrix with a name that makes more sense
        #dxi_tsint[i,:]=df[:]
        #dxi_tsint[i,:]=f_df
        dxi_tsint[i,:]=d_fine_xi_t
    ############################################################--------------------------------
    
    
    # Also initialize the container of the derivative of r^2*xi_r:
    dr2xi_r = np.zeros(xi_r.shape)
    dr2xi_r2 = np.zeros(deldotxi142.shape)
    
    #### Similar code to calculate the derivative to the one used for xi_t on line 141
    #### This thime there's no need to find a Legendre solution.
    for i in range(len(nro_ang)):
        temp_val = xi_r[:,i]*((RS[:,i])**2)
        df = (temp_val[1:-1] - temp_val[0:-2])/(RS[1:-1,i]-RS[0:-2,i])
        f_df = np.interp(RS[:,i],RS[0:-2,i],df)
    
        # Definition of the function that I will derive:
        #f = interp.interp1d(RS[:,i],xi_r[:,i]*((RS[:,i])**2))
        #df = np.zeros(len(RS[:,0]))
        #for j in range(len(df)):
            #try:
                # Find the derivative at the RS points:
                #df[j] = derivative(f,RS[j,i],dx=1e-6)
            #except:
                #a=1
        
        # Store the derivatives
        #dr2xi_r[:,i]=df[:]
        dr2xi_r[:,i]=f_df
        dr2xi_r2[:,i]=df
    ############################################################--------------------------------
    # The r-component of EQ14 can be calculated now:
    RS_mid = 0.5*(RS[1:-1,:] + RS[0:-2,:])
    r_comp = (1./RS[:,:]**2)*(dr2xi_r[:,:])
    r_comp2 = (1./RS_mid**2)*(dr2xi_r2[:,:])
            
    
    global t_comp,t_comp2
    # The theta-component calculation follows:
    for i in range(len(nro_ang)):
        theta = np.deg2rad(nro_ang[i])
        if par == "EVEN":
            t_comp = (1./(RS[:,i]*np.sin(theta)))*(dxi_tsint[:,i])    
        else:
            t_comp = (1./(np.sin(theta)))*(dxi_tsint[:,i])
            
        t_comp2 = 0.5*(t_comp[1:-1] + t_comp[0:-2])
        # Final calculation of deldotxi for EQ14:    
        deldotxi14[:,i]= r_comp[:,i]+t_comp
        deldotxi142[:,i]= r_comp2[:,i]+t_comp2

        
    
    #np.savetxt(modefname+"_deldotxi_eq10",deldotxi10)
    #np.savetxt(modefname+"_deldotxi_eq14",deldotxi14)
    
    for i in range(len(nro_ang)):
        pl1, = plt.plot(RS[:,i],deldotxi10[:,i],"b",label=r"$\nabla\cdot xi$ Eq. 10")
        pl2, = plt.plot(RS[:,i],deldotxi14[:,i],"r",label=r"$\nabla\cdot xi$ Eq. 14")
        plt.legend([pl1,pl2],[r"$\nabla\cdot xi$ Eq. 10",r"$\nabla\cdot xi$ Eq. 14"],loc="best")
        plt.grid(b='on')
        plt.xlim(0,2)
        plt.xlabel(r"R [R$_{\odot}$]")
        plt.ylabel(r"$\nabla\cdot xi$")
        plt.title(r"$\theta=$"+str((i+1)*10))
        #plt.savefig("ddotxi_ang_"+str(i))
        plt.clf()
    
    
    for i in range(len(nro_ang)):
        pl1, = plt.plot(RS[:,i],deldotxi10[:,i],"b",label=r"$\nabla\cdot xi$ Eq. 10")
        pl2, = plt.plot(RS[:,i],deldotxi14[:,i],"r",label=r"$\nabla\cdot xi$ Eq. 14")
        #pl2, = plt.plot(RS[:,i],deldotxi14[:,i],"purple",label="Interp R")
    

    plt.legend([pl1,pl2],[r"$\nabla\cdot xi$ Eq. 10",r"$\nabla\cdot xi$ Eq. 14"],loc="best")
    plt.grid(b='on')
    plt.xlim(0,2)
    plt.xlabel(r"R [R$_{\odot}$]")
    plt.ylabel(r"$\nabla\cdot xi$")
    plt.title("All angles")
    #plt.savefig("ddotxi_allang")

    """
    
    """
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
    #dt_t_e13 = -(g3m1)*deldotxi10_e13
    
    print dt_t[-1,:]
    #print dt_t_e13[-1,:]
    
    return xi_r,xi_t,dt_t,dPhi*PIGR/Rsun,RS,dP,sigma,VS

def norm_and_scale(xi_r,xi_t,dt_t,ZG,norm_f,scale,depth,reese,sig):
    global folder, rotorc_f
    
    t_reese = scale/(sig**2)
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
    
    tmp=np.max(np.abs(xi_r[:,:]))
    if reese:
        """
        for i in np.arange(-(depth),0,1):
            i_max=np.argmax(np.abs(xi_r[i,:]))
            xi_r_n[i,:] =  xi_r[i,:]/xi_r[i,i_max]
            xi_t_n[i,:] =  xi_t[i,:]/xi_r[i,i_max]
            dt_t_n[i,:] =  dt_t[i,:]/xi_r[i,i_max]
            ZG_n[i,:] =  ZG[i,:]/xi_r[i,i_max]
        """
        
        xi_r_n = xi_r[-(depth+1):-1,:]*t_reese*(1./tmp)
        xi_t_n = xi_t[-(depth+1):-1,:]*t_reese*(1./tmp)
        dt_t_n = dt_t[-(depth+1):-1,:]*t_reese*(1./tmp)
        ZG_n = ZG[-(depth+1):-1,:]*t_reese*(1./tmp)
        print "Reese normalization: "+str(t_reese)
        
        
    
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
    

"""
import pyLegendre_anybf as pyL
import seaborn as sns
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':18})  


mode = "MODE14"

#xi_r,xi_t,dt_t,zg,r,zp,sigma = calcdeldotxi("EVEN",["2p5"],0,"",mode)
xi_r,xi_t,dt_t,dPhi,RS,dP,sigma,VS = calcdeldotxi("ODD",["2p5"],0,"",mode)
"""

#f = open("M2p5_V0_"+mode+"_perturbations","w")
#f.write("M2p5, V=0, "+ mode+" sigma = "+str(sigma)+"\n")
#f.write("n_angle, r, xi_r, xi_t, dT/T\n")
#for i in range(len(r[:,0])-10,len(r[:,0])):
#    for j in range(len(xi_r[-1,:])):
        #print "%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j])
        #f.write("%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j]))

#f.close()


#xi_r_n,xi_t_n,dt_t_n = norm_and_scale(xi_r,xi_t,dt_t,norm_f,scale,depth)

#xi_r_rot,xi_t_rot,dt_t_rot = to_rotorc(xi_r_n,xi_t_n,dt_t_n)

#container_t = pyL.legendre(dt_t[-1,:],8)
#container_r = pyL.legendre(xi_r[-1,:],8)
#bob_deldotxi_mode13 = np.genfromtxt("BOB_June8_2015/M2p5_V0_mode13_surf_perturbations",skip_header=2)
