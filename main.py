# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 01:04:06 2015

@author: diego


"""

import numpy as np
import pylab as plt
import glob,subprocess,os
import del_dot_xi as ddxi
import pyNRO_1_mode as pyNRO
import scipy.integrate as scint
import legendre_interp as lint
import seaborn as sns
from timeit import default_timer as timer

from main_modules import find_vel,find_index,find_sigma,find_name,index_by_name,emode,phaser,phaser2,find_dr_ofmax

vis_path = os.getcwd()

global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/castaned/Documents/"
else:
    homedir = "/home/castaned/"


def run_visc_pert(model,vel,mode,par,sigma,reese,force_f,minL,phase=0.,ampl=1.,tlmax=0.,tcut=False):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    global kind
    
    norm_f = True
    ssc = 0.02
    if model[0]=="2p5": ssc = 0.01
    if kind=="g":
        scale = ssc
    else:
        scale = ssc/(sigma**2)
    scale =0.001*drtempmax
    depth = 10 #radial zones down from the surface
    forcenro = False
    #---------------------------------------------------
    
    #clic_run = True
    
    
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    
    rotorc_f = homedir+"From_Bob/Delta_Scuti_2010/"+model[0]+"Msun/"+model[0]+"Msun_V"+v+"/"
    
    bob_bin = homedir+"From_Bob/clotho_disc10_bin/"
    
    static_m = homedir+"ROTORCmodels/visibilities/"
    
    
    
    
    #check if visibility file from pulset exits for the selected model:
    
    #v_file = glob.glob(rotorc_f+"visibility_file")
    
    
    ###### FULL MODE FILE GENERATION:
    
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
    tfreq = temp_freqs[mode-1]
    if force_f==True:
        tfreq = sigma
        print sigma
    #tfreq = 1.59692
    #print pyNRO.run_nro(tfreq,folder,model[0],v,par,mode)
    
    
    
    where =static_m+model[0]+'Msun/V'+v+"/MODE_"+par+"_"+str(mode)
    
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+v):
        os.makedirs(static_m+model[0]+'Msun/'+'V'+v)
        
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+v+"/MODE_"+par+"_"+str(mode)):
        os.makedirs(where)
    
    if os.path.isfile(where+'/MODE_'+par+'_'+str(mode))==False or forcenro==True:
        #check if visibility file from pulset exits for the selected model:
        
        #v_file = glob.glob(rotorc_f+"visibility_file")
        if True:
            os.chdir(rotorc_f)
            r_mod_name = glob.glob("*_ZAMS")
            subprocess.call(["cp",r_mod_name[0],"orcmod_pulset"])
            subprocess.call([bob_bin+"pulsetnonadb.exe"])
            print "Generated visibility_file in "+rotorc_f
            subprocess.call([bob_bin+"pulset_gammas.exe"]) #Generates Cmod with gammas
            subprocess.call(["cp","Cmod",folder+"Dmod_"+model[0]+"M_V"+v])
            print "Generated new Dmod"
            #print(glob.glob("*_ZAMS"))
            os.chdir(vis_path)
        
        print "Dmod with GAMMAS generated!"
        
        nro_stat = 1
        idx_iter = 0
        while nro_stat==1:
        #RUN NRO!
            nro_stat = pyNRO.run_nro(tfreq-idx_iter*1e-4,folder,model[0],v,par,mode,homedir)
            idx_iter += 1

            
        subprocess.call(['mv',folder+'MODE_'+par+'_'+str(mode),where+'/MODE_'+par+'_'+str(mode)])
        print "Mode file generation complete!"
    else:
        print "Mode file found! Not running NRO..."


    
    
    ###### del_DOT_xi>
    #modeloc = static_m+model[0]+'Msun/V'+vels[vel]+"/"
    modeloc = where+"/"
    modefname = 'MODE_'+par+'_'+str(mode)
    #modefname = 'MODE13'
    global s_model
    s_model = np.genfromtxt(glob.glob(static_m+model[0]+'Msun/'+model[0]+'Msun_V'+v+"*")[0])
    
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    global xi_r,xi_t,dt_t,zg,r,zp,cs,xi_dot_g
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    xi_r,xi_t,dt_t,zg,r,zp,sig,cs,xi_dot_g = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname)
            
    xi_r_n,xi_t_n,dt_t_n,zg_n,dr_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,r,zg,norm_f,scale,depth,reese,sig,par)

    if phase!=0:
        #dt_t_n = phaser(dt_t_n,phase,1.,tlmax)
        #xi_r_n = phaser(xi_r_n,0.,1.,tlmax)
        dt_t_n,psi_T,psi_L,mx = phaser2(dt_t_n,phase,np.max(np.abs(dt_t[-depth])),np.max(np.abs(xi_r[-depth])))
        xi_r_n = phaser(xi_r_n,0.,1.,mx+np.pi)
        print "Theta L_max-----------> ",np.rad2deg(mx),"deg"
        np.savetxt(modeloc+"phases.txt",[psi_T,psi_L,np.rad2deg(mx+np.pi)],header="psi_T[deg],psi_L[deg],theta_max_L[deg]")
        
        

    a_r = scint.trapz(dt_t_n[-1,:])
    
    inv_f = 1.
    if a_r<0:
        print "-T area"
        #xi_r_n *= -1.
        #xi_t_n *= -1.
        #dt_t_n *= -1.
        #zg_n *= -1.
        inv_f = -1.
        
    if minL==True:
        inv_f *= -1.

    #xi_r_rot,xi_t_rot,dt_t_rot,zg_rot = ddxi.to_rotorc(xi_r_n,xi_t_n,dt_t_n,zg_n)
    
    if par == "EVEN":
        xi_r_fine = lint.leg_interp(xi_r_n[:,:],8,"EVEN")
        xi_t_fine = lint.leg_interp(xi_t_n[:,:],8,"EVEN")
        dt_t_fine = lint.leg_interp(dt_t_n[:,:],8,"EVEN")
    else:
        xi_r_fine = lint.leg_interp(xi_r_n[:,:],8,"OE")
        xi_t_fine = lint.leg_interp(xi_t_n[:,:],8,"OE")
        dt_t_fine = lint.leg_interp(dt_t_n[:,:],8,"OE")
        
    global rs_n,rs_rot
    rs_n = np.empty(xi_r_fine.shape)
    for d in range(depth):
        rs_n[d,:] = np.interp(np.linspace(0,90,100),np.linspace(10,80,8),r[-d-1,:])
    
#    dr_n = xi_r_fine*rs_n
#    for i in range(len(dr_n[:,0])):
#        imax = np.argmax(np.abs(dr_n[i,:]))
#        dr_n[i,:] = dr_n[i,:]/(dr_n[i,imax]/2.28534026)
        
    
    
    p_angles = np.empty(np.shape(xi_r_n))
    p_angles_fine = np.empty(np.shape(xi_r_fine))
    rot_ang = np.arange(4.5,90.,9)
    xi_r_rot = np.empty((depth,len(rot_ang)))
    xi_t_rot = np.empty((depth,len(rot_ang)))
    dt_t_rot = np.empty((depth,len(rot_ang)))
    rs_rot = np.empty((depth,len(rot_ang)))
    for d in range(depth):
        p_angles[d,:] = (np.linspace(10,80,8))*(1.+xi_t_n[d,:])
        p_angles_fine[d,:] = (np.linspace(0,90,100))*(1.+xi_t_fine[d,:])
        xi_r_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],xi_r_fine[d,:]) 
        xi_t_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],xi_t_fine[d,:])
        dt_t_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],dt_t_fine[d,:])
        rs_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],rs_n[d,:])
        
    
    
#    
#    #plt.plot(p_angles_fine[2,:],xi_t_fine[2,:])
#    if dr_n[0,-1]>0:
#        ivrt = 1.
#    else:
#        ivrt = -1.
#    vr = ivrt*dr_n[0,:]
#    #vr = ivrt*dt_t_fine[0,:]
#    #lblb = find_vel(model,vel)+" km s$^{-1}$ "+find_name(model,vel,par,mode)
#    lblb = find_vel(model,vel)+" km s$^{-1}$ "
#    #lblb = r"$\ell$ = "+find_name(model,vel,par,mode).split()[0]
#    etsy = "-"
#    if find_name(model,vel,par,mode).split()[0] == "2": etsy = "-"
#    if par == "EVEN":
#        #plt.plot(p_angles_fine[2,:],ivrt*lint.leg_interp(vr[:,:],8,"EVEN")[0,:],label=lblb)
#        plt.plot(p_angles_fine[2,:],vr,ls=etsy,label=lblb)
#        
#    else:
#        #plt.plot(p_angles_fine[2,:],ivrt*lint.leg_interp(vr[:,:],8,"OE")[0,:],label=lblb)
#        plt.plot(p_angles_fine[2,:],vr,ls=etsy,label=lblb)
#    #plt.plot(p_angles_fine[2,:],lint.leg_interp(xi_r[-depth::,:],8,"OE")[2,:],label=find_name(model,vel,par,mode))
#    #plt.plot(p_angles[2,:],dt_t[-depth+2,:],"o",mfc="white",mec="k",mew=1)
#    #plt.grid()
#    plt.xlim(0,90)
#    plt.yticks()
#    plt.xlabel("Colatitude [deg]")
#    plt.ylabel(r"$\delta$R [R$\odot$]")
#    #plt.ylabel(r"$\delta$T/T")
#    plt.legend(loc="best")
#    #plt.title("Mode:" + find_name(model,vel,par,mode))
#    #ax.yaxis.set_major_formatter(majorFormatter) 
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#    #plt.title(r"M="+model[0]+"M$_{\odot}$, V="+v+"km s$^{-1}$")
#    plt.title(r"M="+model[0]+"M$_{\odot}$")
    
       
    
    print "Generating fine clic model..."    
    pmodels_fine = gen_fine_clic(p_angles_fine,xi_r_fine,xi_t_fine,dt_t_fine,rs_n,dr_n,s_model,model,depth,inv_f,par,tcut)

        
    ############ find the perturbed models:
    global omega
    G = 6.67259e-8
    mass = (float((model[0]).replace("p",".")))*1.99e33
    theta = np.deg2rad(s_model[:,0])
    radius = s_model[:,1]*6.96e10
    vrot = s_model[:,4]*1e5
    omega = vrot/(radius*np.sin(theta))
    omega = np.average(omega[1:-1])
    

    
    #g_r = (G*mass/(radius**2))-((omega**2)*(radius*np.sin(theta)))*np.sin(theta)
    #g_theta = ((omega**2)*(radius*np.sin(theta)))*np.cos(theta)

    #geff = (g_r**2 + g_theta**2)**0.5        
        
    
    pert_r = np.empty((len(s_model[:,0]),depth))
    pert_t = np.empty((len(s_model[:,0]),depth))
    pmodels = []
    for i in range(depth):
        pert_r[:,i] = s_model[:,1]*(1.+inv_f*xi_r_rot[-i,:])
        pert_t[:,i] = s_model[:,2]*(1.+inv_f*dt_t_rot[-i,:])
        g_pert = (G*mass/((pert_r[:,i]*6.96e10)**2))-((omega**2)*((pert_r[:,i]*6.96e10)*np.sin(theta)))*np.sin(theta)
        for j in range(len(g_pert)):
            if np.log10(g_pert[j])>4.33: g_pert[j]=10**(4.33)
        tmodel = 1.*s_model
        tmodel[:,1] = pert_r[:,i]
        tmodel[:,2] = pert_t[:,i]
        tmodel[:,3] = np.log10(g_pert)
        pmodels.append(tmodel)

    if par=="ODD":
        pert_r = np.empty((2*len(s_model[:,0]),depth))
        pert_t = np.empty((2*len(s_model[:,0]),depth))
        
        pmodels = []
        odd_angles = np.arange(4.5,180.,9)
        for i in range(depth):
            pert_r[0:10,i] = s_model[:,1]+inv_f*xi_r_rot[-i,:]*rs_rot[-i,:]
            pert_t[0:10,i] = s_model[:,2]*(1.+inv_f*dt_t_rot[-i,:])
            pert_r[10:20,i] = s_model[::-1,1]*(1.-inv_f*xi_r_rot[-i,::-1])
            pert_t[10:20,i] = s_model[::-1,2]*(1.-inv_f*dt_t_rot[-i,::-1])
            g_pert = (G*mass/((pert_r[:,i]*6.96e10)**2))-((omega**2)*((pert_r[:,i]*6.96e10)*np.sin(odd_angles)))*np.sin(odd_angles)
            for j in range(len(g_pert)):
                if np.log10(g_pert[j])>4.33: g_pert[j]=10**(4.33)
            tmodel = np.empty((20,5))
            tmodel[:,0] = odd_angles[:]
            tmodel[:,1] = pert_r[:,i]
            tmodel[:,2] = pert_t[:,i]
            #tmodel[0:10,3] = s_model[:,3]
            #tmodel[10:20,3] = s_model[::-1,3]
            tmodel[:,3] = np.log10(g_pert)
            tmodel[0:10,4] = s_model[:,4]
            tmodel[10:20,4] = s_model[::-1,4]
                
            pmodels.append(tmodel)
    
#    plt.plot(pmodels[-1][:,0],s_model[:,1],"--",color="0.5")
#    
#    global old_modes
#    import seaborn as sns
#    sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
#                              'ytick.labelsize': 12,'xtick.labelsize': 12,
#                              'legend.fontsize': 16,'axes.titlesize':18,'font.size':14})
#    #plt.plot(pmodels[2][0:-1,0],np.diff(pmodels[2][:,1])/(pmodels[2][1,0]-pmodels[2][0,0]),label=r"$\ell=$"+find_name(model,vel,par,mode).strip().split()[0])
#    #plt.plot(pmodels[2][:,0],pmodels[2][:,1],label="old way")
#    #plt.plot(pmodels_fine[2][0:-1,0],np.diff(pmodels_fine[2][:,1])/(pmodels_fine[2][1,0]-pmodels_fine[2][0,0]))
#    plt.plot(pmodels_fine[2][:,0],pmodels_fine[2][:,2],"-",lw=1.5,label=r"$\ell=$"+find_name(model,vel,par,mode).strip().split()[0])
#    #plt.plot(s_model[:,0],s_model[:,1],label="Static")
#    #plt.vlines(90,9100,9300)
#    o_lims = plt.ylim()
#    #plt.vlines(90,min(pmodels_fine[2][:,2])-100,max(pmodels_fine[2][:,2])+100)
#    plt.grid()
#    plt.xlim(0,180)
#    #plt.ylim(o_lims[0],o_lims[1])
#    plt.xlabel("Colatitude [deg]")
#    plt.ylabel("Radius [M$_{\odot}$]")
#    #plt.ylabel("Temperature [K]")
#    plt.legend(loc="best")
#    #plt.title("Mode: " + find_name(model,vel,par,mode))
    #plt.title(r"Perturbed T$_{\mathrm{eff}}$ - M="+model[0]+"M$_{\odot}$, V="+v+"km s$^{-1}$")
    
    
    return modeloc,pmodels,pmodels_fine
    
def gen_fine_clic(p_angles_fine,xi_r,xi_t,dt_t,rs,dr_n,s_model,model,depth,inv_f,par,tcut):
    clic_theta_grid = 200
    clic_theta = np.linspace(0,180,clic_theta_grid)
    global clic_xi_r,clic_xi_t,clic_dt_t
    
    clic_xi_r = np.zeros((depth,clic_theta_grid))
    clic_xi_t = np.zeros((depth,clic_theta_grid))
    clic_dt_t = np.zeros((depth,clic_theta_grid))
    clic_rs = np.zeros((depth,clic_theta_grid))
    if tcut:
        cut=True
        print "--- CUT ---"
    else:
        cut=False
    for d in range(depth):
        if cut:
            xi_r[d,0:12] = 1.*xi_r[d,12]
            xi_t[d,0:12] = 1.*xi_t[d,12]
            dt_t[d,0:12] = 1.*dt_t[d,12]
            rs[d,0:12] = 1.*rs[d,12]
            if par=="EVEN":
                xi_r[d,-9::] = 1.*xi_r[d,-10]
                xi_t[d,-9::] = 1.*xi_t[d,-10]
                dt_t[d,-9::] = 1.*dt_t[d,-10]
                rs[d,-9::] = 1.*rs[d,-10]
        clic_xi_r[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),xi_r[d,:])
        clic_xi_t[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),xi_t[d,:])
        clic_dt_t[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),dt_t[d,:])
        clic_rs[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),rs[d,:])
    
    
    for i in range(depth):
        if par=="EVEN":
            clic_xi_r[i,100::] = clic_xi_r[i,::-1][100::]
            clic_xi_t[i,100::] = clic_xi_t[i,::-1][100::]
            clic_dt_t[i,100::] = clic_dt_t[i,::-1][100::]
        else:
            clic_xi_r[i,100::] = -1.*clic_xi_r[i,::-1][100::]
            clic_xi_t[i,100::] = -1.*clic_xi_t[i,::-1][100::]
            clic_dt_t[i,100::] = -1.*clic_dt_t[i,::-1][100::]
            
        clic_rs[i,100::] = clic_rs[i,::-1][100::]
            

    #plt.plot(clic_theta[0:200],clic_rs[-2,0:200])

    ############ find the perturbed models:
    # Take static ROTORC model to CLIC grid:
    global s_model_c
    s_model_c = np.zeros((clic_theta_grid,len(s_model[0,:])))
    s_model_c[:,0] = clic_theta
    for i in range(1,len(s_model[0,:])):
        s_model_c[0:100,i] = np.interp(clic_theta[0:100],s_model[:,0],s_model[:,i])
        s_model_c[100::,i] = s_model_c[::-1,i][100::]

    G = 6.67259e-8
    mass = (float((model[0]).replace("p",".")))*1.99e33
    theta = np.deg2rad(s_model_c[:,0])
    radius = s_model_c[:,1]*6.96e10
    vrot = s_model_c[:,4]*1e5
    omega = vrot/(radius*np.sin(theta))
    omega = np.average(omega[50:-50])
    

    a_r = scint.trapz(clic_dt_t[-1,:])
    
#    inv_f = 1.
#    if a_r<0:
#        print "fine -T area"
#        inv_f = -1.    
        
    global pmodels_fine
    pert_r = np.empty((len(s_model_c[:,0]),depth))
    pert_t = np.empty((len(s_model_c[:,0]),depth))
    pmodels_fine = []
    for i in range(depth):
        #pert_r[:,i] = s_model_c[:,1]*(1.+inv_f*clic_xi_r[-i,:])
        pert_r[:,i] = s_model_c[:,1]+inv_f*clic_rs[-i,:]*clic_xi_r[-i,:]
        pert_t[:,i] = s_model_c[:,2]*(1.+inv_f*clic_dt_t[-i,:])
        for j in range(len(pert_t[:,i])):
            if pert_t[j,i]<7500.: pert_t[j,i]=7500.
        g_pert = (G*mass/((pert_r[:,i]*6.96e10)**2))-((omega**2)*((pert_r[:,i]*6.96e10)*np.sin(theta)))*np.sin(theta)
        for j in range(len(g_pert)):
            if np.log10(g_pert[j])>4.33: g_pert[j]=10**(4.33)
        tmodel = 1.*s_model_c
        """
        pert_r[0:11,:] = pert_r[11,:]
        pert_r[190:200,:] = pert_r[189,:]
        
        pert_t[0:11,:] = pert_t[11,:]
        pert_t[190:200,:] = pert_t[189,:]
        """
        tmodel[:,1] = pert_r[:,i]
        tmodel[:,2] = pert_t[:,i]
        tmodel[:,3] = np.log10(g_pert)
        pmodels_fine.append(tmodel)
        
    #plt.plot(clic_dt_t[0,:])
    return pmodels_fine
    
def make_contour(r,var,model,vel,par,mode,x):
    global new_r,new_zp
    import seaborn as sns
    sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,'axes.titlesize':18})

    t_var = np.empty(var[x::,:].shape)
    t_r = np.empty(t_var.shape)
    for i in range(len(r[0,:])):
        t_var[:,i] = var[x::,i]
        t_r[:,i] = np.linspace(r[0,i],r[-1,i],len(t_var[:,i]))
        
    var = 1.*t_var
    r = 1.*t_r
    
    new_r = lint.leg_interp(r[0:-2,:],8,"EVEN")
    new_zp = lint.leg_interp(var[0:-2,:],8,"EVEN")
    #levels = np.linspace(np.min(new_zp),np.max(new_zp), 40)
    
    #####REMOVE 0-10deg at the pole:    
    new_r = new_r[:,11::]
    new_zp = new_zp[:,11::]
    theta = np.linspace(0,np.deg2rad(90),89)
    #########################################
    
    #theta = np.linspace(0,np.deg2rad(90),100)
    newt,n_r = np.meshgrid(theta,new_r[:,0])

#    np.savetxt("/home/castaned/x.dat",(new_r*np.sin(newt)))
#    np.savetxt("/home/castaned/y.dat",(new_r*np.cos(newt)))
#    np.savetxt("/home/castaned/z.dat",(new_zp))
    
    plt.contourf((new_r*np.sin(newt)),(new_r*np.cos(newt)),new_zp, 100, cmap=plt.cm.gnuplot,vmax=np.max(new_zp), vmin=np.min(new_zp))
    plt.contourf((new_r*np.sin(newt)),-(new_r*np.cos(newt)),new_zp, 100, cmap=plt.cm.gnuplot,vmax=np.max(new_zp), vmin=np.min(new_zp))
    plt.contourf(-(new_r*np.sin(newt)),-(new_r*np.cos(newt)),new_zp, 100, cmap=plt.cm.gnuplot,vmax=np.max(new_zp), vmin=np.min(new_zp))
    plt.contourf(-(new_r*np.sin(newt)),(new_r*np.cos(newt)),new_zp, 100, cmap=plt.cm.gnuplot,vmax=np.max(new_zp), vmin=np.min(new_zp))
    #CSl = plt.contour((new_r*np.sin(newt)),(new_r*np.cos(newt)),new_zp, 20, colors="k")
    #CS = plt.contourf((new_r*np.cos(newt)), (new_r*np.sin(newt)), new_zp, cmap=plt.cm.Spectral,levels=levels)
    #plt.axes().set_aspect('equal')
    #plt.xlim(np.ceil(r[-1,-1]))
    plt.ylim(-plt.xlim()[1],plt.xlim()[1])
    #plt.xlabel("Radius [R$_{\odot}$]")
    #plt.ylabel("Radius [R$_{\odot}$]")
    plt.axis('off')
    try:
        m = (find_name(model,vel,par,mode)).strip()
    except:
        m = "freq_%.5f" % find_sigma(model,vel,par,mode)
    v = find_vel(model,vel)

    plt.title("M"+model[0]+"_V"+v+" - mode: "+ m)
    return
    
def save_pert(modeloc,model,vel,par,mode,name):
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    global xi_r,xi_t,dt_t,zg,r
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    print os.getcwd()
    f = open(modeloc+"M"+model[0]+"_V"+str(vel)+"_"+name+"_perturbations","w")  
    f.write("M"+model[0]+", V="+str(vel)+", Mode freq = %.5f\n" % find_sigma(model,vel,par,mode))
    f.write("n_angle, r, xi_r, xi_t, dT/T\n")
    for i in range(len(r[:,0])-10,len(r[:,0])):
        for j in range(len(xi_r[-1,:])):
            #print "%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j])
            f.write("%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j])) 
            
    f.close()
    return
         
def list_gen(ells,nmin,nmax,tpe):
    ls = []
    for j in ells:
        for i in range(nmin,nmax+1):
            ls.append(str(j)+" "+tpe+str(i))
            
    return ls
"""

------------------------------------------------------
------------------------------------------------------

"""
def run_pipeline(modes,sub_dir="",minL=False,mm="1p875",vv=0,tph=0,drscale=None,tcut=False):
    global start,model,vel,old_modes,modeloc,pmodels,pmodels_fine,modes_not_found,mds
    #MODE info:
    
    model = [mm]
    vel = vv #index, not velocity!
    #modes = list_gen([0,1,2,3,4,5,6,7,8,9,10],5,5,"p")
    #list_gen([0,1,2,3],1,1,"p")
    #modes = list_gen([2],8,8,"p")
    #["1 p1","1 p2","1 p3","1 p4","1 p5", "1 p6", "1 p7", "1 p8", "1 p9"]
    #["0 1H","0 2H","0 3H","0 4H","0 5H","0 6H","0 7H","0 8H","0 9H","0 10H"]
    #par = emode(mm,vv,modes[0]).parity
    
    
    clic = False
    new_mags = False
    save_perturbation = False
    only_mags = False

    
    incl = list(np.arange(0,95,10))
    #incl = [90]
    
    
    reese = False
    force_f = False
    f_freq = 1.22
       
    phase = tph
    ampl = 1./2.2
    tlmax = 4.14 #radians (sigma*t of phased Lmax)
    
    #### contour plotting option:
    plot_contour = True
    excl = 20 #How many core zones to exclude
    ####
    
    
    
    freqs = [6.61043]
    
    dt_grid = 250
    
    old_modes = list(modes)
    mode_by_freq = False
    mode_by_name = False
    if type(modes[0]) == float:
        mode_by_freq = True
    
    if type(modes[0]) == str:
        mode_by_name = True
        
    if mode_by_freq==True:
            modes = find_index(freqs,model,vel,par)
    
    parity = []
    modesx = []
    if mode_by_name==True:
        for i in range(len(modes)):
            if int(modes[i].split()[0])%2==0:
                par = "EVEN"
            else:
                par = "ODD"
            
            parity.append(par)
            idx_tmp = index_by_name(model,vel,par,modes[i])
            print modes[i],"->",idx_tmp
            modesx.append(idx_tmp)
    
            
    start = timer()
    klm = 0
    rem = 0
    
    modes_not_found = []
    for mode in modesx:
        par = emode(mm,vv,mds[modesx.index(mode)]).parity
        if mode == 999:
            modes_not_found.append(old_modes[modes.index(mode)])
            continue
        
        rem += 1
        
        if mode_by_name == True:
            par = parity[klm]
            klm +=1
            
        if force_f==True:
            sigma=f_freq
        else:
            #sigma = find_sigma(model,vel,par,mode)
            sigma = emode(mm,vv,old_modes[modesx.index(mode)]).frequency
          
            
        #Find the perturbed models!
    
        modeloc,pmodels,pmodels_fine = run_visc_pert(model,vel,mode,par,sigma,reese,force_f,minL,phase,ampl,tlmax,tcut)
        modeloc += sub_dir
        if not os.path.exists(modeloc):
            os.makedirs(modeloc)
        
        print modeloc
        print mds[modesx.index(mode)]
        nrzone = len(dt_t[:,0])    
        
        if plot_contour == True:
            make_contour(r,-zp*cs,model,vel,par,mode,excl)
            
        if save_perturbation==True:
            save_pert(modeloc,model,vel,par,mode,old_modes[modesx.index(mode)])
        
        for i in range(len(pmodels)):
            np.savetxt(modeloc+"model_MODE_"+str(mode)+"_r"+str(nrzone - i),pmodels[i],'%.16e')
            np.savetxt(modeloc+"fine_model_MODE_"+str(mode)+"_r"+str(nrzone - i),pmodels_fine[i],'%.16e')
        
        rzone = 2 # Radial zone from the surface
        rzone = 490
        print rzone
        
        if clic==True:
            #Ugly way of running pyCLIC:
            clic_folder = homedir+"pyCLIC/test/"
            
            if only_mags==False:
                subprocess.call(['mkdir',modeloc+'tmp'])
                os.chdir(modeloc+'tmp')
                tmpstr = "cp -r "+clic_folder+"* ."
                subprocess.call(tmpstr,shell=True)
                
                for i in incl:
                    
                    if new_mags==True:
                        if (os.path.isfile(modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==True:
                            subprocess.call(['rm',modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])
                        else:
                            print "No mag. file to delete!" 
                        if (os.path.isfile(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==True:
                            subprocess.call(['rm',modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])
                        else:
                            print "No walraven_mag. file to delete!"
    
                    
                    
                    subprocess.call(['cp',modeloc+"model_MODE_"+str(mode)+"_r"+str(rzone),modeloc+'tmp'])
                    subprocess.call(['cp',modeloc+"fine_model_MODE_"+str(mode)+"_r"+str(rzone),modeloc+'tmp'])
                        
                    if (os.path.isfile(modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==False:
                        myfile = open(modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode), "w")
                        myfile.write('#%5s %5s %9s %11s %11s %11s %11s %11s\n' % ("incl","n_r","lum","u","b","v","r","i"))
                        myfile.close()
                        subprocess.call(['cp',modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc+'tmp'])
                    else:
                        subprocess.call(['cp',modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc+'tmp'])
                        
                    if (os.path.isfile(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==False:
                        myfile = open(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode), "w")
                        myfile.write('#%5s %5s %9s %11s %11s %11s %11s %11s\n' % ("incl","n_r","lum","w","u","l","b","v"))
                        myfile.close()
                        subprocess.call(['cp',modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc+'tmp'])
                    else:
                        subprocess.call(['cp',modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc+'tmp'])
                    
                    
                    
                    os.chdir(modeloc+'tmp')
                    #sys.path.append(modeloc+'tmp')
                    #print sys.path
                    import pyCLIC_main as pyclic
                    #sys.path.append(vis_path)
                    os.chdir(vis_path)
                    import color_magnitudes as cmag
                    os.chdir(modeloc+'tmp')
                    
                    
                    pyclic.run_CLIC("model_MODE_"+str(mode)+"_r"+str(rzone),[i],False,3000.,7499.,10.0,par,dt_grid,fine_model=True,savexi=True)
                    
                    cmag.calc_mags('outputflux_i'+str(i)+'.final',[i],mode,rzone)
                    cmag.calc_walraven('outputflux_i'+str(i)+'.final',[i],mode,rzone)
    
    
                    subprocess.call(['cp','outputflux_i'+str(i)+'.final',modeloc+"outputflux_i"+str(i)+"_MODE_"+str(mode)+"_r"+str(rzone)])
                    subprocess.call(['cp',"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc])
                    subprocess.call(['cp',"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode),modeloc])
                    
                    subprocess.call(['rm',"model_MODE_"+str(mode)+"_r"+str(rzone)])
                    subprocess.call(['rm','outputflux_i'+str(i)+'.final'])    
                    subprocess.call(['rm',"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])    
                    subprocess.call(['rm',"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])
                    
                    os.chdir(vis_path)
                    
                subprocess.call(['rm',"-r",modeloc+'tmp'])
                    
            else:
                os.chdir(vis_path)
                import color_magnitudes as cmag
                os.chdir(modeloc)
                for i in incl:
                    
                    if new_mags==True:
                        if (os.path.isfile(modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==True:
                            subprocess.call(['rm',modeloc+"magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])
                        else:
                            print "No mag. file to delete!" 
                        if (os.path.isfile(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==True:
                            subprocess.call(['rm',modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)])
                        else:
                            print "No walraven_mag. file to delete!"
                
    
                    if (os.path.isfile("magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==False:
                        myfile = open("magnitudes_MODE_"+str(mode), "w")
                        myfile.write('#%5s %5s %9s %11s %11s %11s %11s %11s\n' % ("incl","n_r","lum","u","b","v","r","i"))
                        myfile.close()
        
                        
                    if (os.path.isfile(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode)))==False:
                        myfile = open(modeloc+"walraven_magnitudes"+"_i"+str(i)+"_MODE_"+str(mode), "w")
                        myfile.write('#%5s %5s %9s %11s %11s %11s %11s %11s\n' % ("incl","n_r","lum","w","u","l","b","v"))
                        myfile.close()
                    
                    
                    #cmag.calc_mags("outputflux_i"+str(i)+"_MODE_"+str(mode)+"_r"+str(rzone),[i],mode,rzone)
                    cmag.calc_walraven("outputflux_i"+str(i)+"_MODE_"+str(mode)+"_r"+str(rzone),[i],mode,rzone)
    
                os.chdir(vis_path)
        print mode
#        plt.plot(0.02*xi_t[:,3]/np.max(m4p4_xi_r[-10,:]),"-",alpha=0.8,label=old_modes[modes.index(mode)])
#        plt.grid()
#        plt.xlabel("Radial Zone")
#        plt.ylabel(r"$\delta$T/T")
#        plt.ylabel(r"$\xi_{r}$")
#        plt.title(r"M=1.875M$\odot$ - V=0 - i=40deg")
#        plt.legend(loc="best")
        
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 12,'axes.titlesize':18,'font.size':14})
start = timer()
#varr = [0,4,6,7,8]
#mselect = "1p875"
#kind = "f"
#r_ord = 0
#mds = list_gen([0,2,3,4,5],r_ord,r_ord,kind)
#mmaxarr = []
#global tempmax
#tempmax = np.array([0,0,0])
#for m in mds:
#    mmax = find_dr_ofmax(emode(mselect,0,m),10,range(10))
#    t = mmax[np.argmax(mmax[:,-1])]
#    if abs(t[2])>tempmax[2]:
#        tempmax = t
#    mmaxarr.append(mmax)
global drtempmax
#drtempmax = 2.28534e+00 #n=5
#drtempmax = 3.3914441477217716 # n=2
drtempmax = 3.111 # n=0
#drtempmax = 4.4445 #8.889 # n=5 M2p5
    
for mnopq in [3]:       
    mselect = "2"
    kind = "p"
    r_ord = mnopq
    varr = [0]
    mds = list_gen([6],r_ord,r_ord,kind)
    
    #with sns.color_palette("rainbow",n_colors=9):
    for i in varr:
         run_pipeline(mds,sub_dir="",mm=mselect,vv=i)
#        run_pipeline(mds,sub_dir="minL/",minL=True,mm=mselect,vv=i)
################## 
#    for j in varr:
#        run_pipeline(mds,sub_dir="phi_max/",mm=mselect,vv=j,tph=123.)
#        run_pipeline(mds,sub_dir="phi_min/",minL=True,mm=mselect,vv=j,tph=123.)
#################        
#    for j in varr:
#        run_pipeline(mds,sub_dir="cut_max/",mm=mselect,vv=j,tph=0.,tcut=True)
#        run_pipeline(mds,sub_dir="cut_min/",minL=True,mm=mselect,vv=j,tph=0.,tcut=True)

#if kind !="g":
#    plt.title(r"M="+mselect+"M$_{\odot}$ - $\ell$ = "+mds[0].split()[0] + ", n = "+mds[0].split()[1].strip("p").strip("h").strip("f").strip("F").strip("g"))
#else:
#    plt.title(r"M="+mselect+"M$_{\odot}$ - $\ell$ = "+mds[0].split()[0] + ", n = -"+mds[0].split()[1].strip("p").strip("h").strip("f").strip("F").strip("g"))

#ody = plt.ylim()
#odx = plt.xlim()
#plt.vlines(10,ody[0],ody[1],linestyles="--",alpha=0.5)
#plt.ylim(ody[0],ody[1])
#plt.hlines(0,odx[0],odx[1],linestyles="--",alpha=0.5)
#plt.xlim(odx[0],odx[1])
#plt.title(r"M="+mselect+"M$_{\odot}$ - V = 0 km/s, n = "+mds[0].split()[1].strip("p").strip("h").strip("f").strip("F").strip("g"))
plt.grid()
end = timer()
#print modes_not_found
print "time:",end-start

