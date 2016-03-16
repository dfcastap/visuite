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
import sys
from timeit import default_timer as timer

from main_modules import find_vel,find_index,find_sigma,find_name,index_by_name

vis_path = os.getcwd()

global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/castaned/Documents/"
else:
    homedir = "/home/castaned/"


def run_visc_pert(model,vel,mode,par,sigma,reese,force_f,minL):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    norm_f = True
    scale = 0.02/(sigma**2)
    #scale = 0.02
    depth = 10 #radial zones down from the surface
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
    
    if os.path.isfile(where+'/MODE_'+par+'_'+str(mode))==False:
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
    global xi_r,xi_t,dt_t,zg,r,zp,cs
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    xi_r,xi_t,dt_t,zg,r,zp,sig,cs = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname)
            
    xi_r_n,xi_t_n,dt_t_n,zg_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,zg,norm_f,scale,depth,reese,sig,par)

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
        
    p_angles = np.empty(np.shape(xi_r_n))
    p_angles_fine = np.empty(np.shape(xi_r_fine))
    rot_ang = np.arange(4.5,90.,9)
    xi_r_rot = np.empty((depth,len(rot_ang)))
    xi_t_rot = np.empty((depth,len(rot_ang)))
    dt_t_rot = np.empty((depth,len(rot_ang)))
    for d in range(depth):
        p_angles[d,:] = (np.linspace(10,80,8))*(1.+xi_t_n[d,:])
        p_angles_fine[d,:] = (np.linspace(0,90,100))*(1.+xi_t_fine[d,:])
        xi_r_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],xi_r_fine[d,:]) 
        xi_t_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],xi_t_fine[d,:])
        dt_t_rot[d,:] = np.interp(rot_ang,p_angles_fine[d,:],dt_t_fine[d,:])
        
    
    
#    
#    #plt.plot(p_angles_fine[2,:],xi_t_fine[2,:])
#    plt.plot(p_angles_fine[2,:],lint.leg_interp(dt_t[-depth::,:],8,"OE")[2,:],label=find_name(model,vel,par,mode))
#    #plt.plot(p_angles[2,:],dt_t[-depth+2,:],"o",mfc="white",mec="k",mew=1)
#    plt.grid()
#    plt.xlim(0,90)
#    plt.yticks()
#    plt.xlabel("Colatitude [deg]")
#    plt.ylabel(r"$\delta$ T/T")
#    #plt.title("Mode:" + find_name(model,vel,par,mode))
#    #ax.yaxis.set_major_formatter(majorFormatter) 
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#    plt.title(r"M="+model[0]+"M$_{\odot}$, V="+v+"km s$^{-1}$")
    
       
    
    print "Generating fine clic model..."    
    pmodels_fine = gen_fine_clic(p_angles_fine,xi_r_fine,xi_t_fine,dt_t_fine,s_model,model,depth,inv_f,par)

        
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
            pert_r[0:10,i] = s_model[:,1]*(1.+inv_f*xi_r_rot[-i,:])
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
    
    
#    global old_modes
#    import seaborn as sns
#    sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
#                              'ytick.labelsize': 12,'xtick.labelsize': 12,
#                              'legend.fontsize': 16,'axes.titlesize':18,'font.size':14})
#    #plt.plot(pmodels[2][0:-1,0],np.diff(pmodels[2][:,1])/(pmodels[2][1,0]-pmodels[2][0,0]))
#    #plt.plot(pmodels[2][:,0],(pmodels[2][:,2]),"o")
#    #plt.plot(pmodels_fine[2][0:-1,0],np.diff(pmodels_fine[2][:,1])/(pmodels_fine[2][1,0]-pmodels_fine[2][0,0]))
#    plt.plot(pmodels_fine[2][:,0],pmodels_fine[2][:,2],lw=1.5,label=old_modes[0])
#    #plt.plot(s_model[:,0],s_model[:,2],label="Static")
#    #plt.vlines(90,9100,9300)
#    o_lims = plt.ylim()
#    plt.vlines(90,min(pmodels_fine[2][:,2])-100,max(pmodels_fine[2][:,2])+100)
#    plt.grid()
#    plt.xlim(0,180)
#    plt.ylim(o_lims[0],o_lims[1])
#    plt.xlabel("Colatitude [deg]")
#    #plt.ylabel("Radius [M$_{\odot}$]")
#    plt.ylabel("Temperature [K]")
#    plt.legend(loc="best")
#    #plt.title("Mode: " + find_name(model,vel,par,mode))
#    plt.title(r"Perturbed T$_{\mathrm{eff}}$ - M="+model[0]+"M$_{\odot}$, V="+v+"km s$^{-1}$")
#    
    
    return modeloc,pmodels,pmodels_fine
    
def gen_fine_clic(p_angles_fine,xi_r,xi_t,dt_t,s_model,model,depth,inv_f,par):
    clic_theta_grid = 200
    clic_theta = np.linspace(0,180,clic_theta_grid)
    global clic_xi_r,clic_xi_t,clic_dt_t
    
    clic_xi_r = np.zeros((depth,clic_theta_grid))
    clic_xi_t = np.zeros((depth,clic_theta_grid))
    clic_dt_t = np.zeros((depth,clic_theta_grid))
    cut = False
    for d in range(depth):
        if cut:
            xi_r[d,0:12] = 1.*xi_r[d,12]
            xi_t[d,0:12] = 1.*xi_t[d,12]
            dt_t[d,0:12] = 1.*dt_t[d,12]
            if par=="EVEN":
                xi_r[d,-9::] = 1.*xi_r[d,-10]
                xi_t[d,-9::] = 1.*xi_t[d,-10]
                dt_t[d,-9::] = 1.*dt_t[d,-10]
        clic_xi_r[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),xi_r[d,:])
        clic_xi_t[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),xi_t[d,:])
        clic_dt_t[d,0:100] = np.interp(clic_theta[0:100],np.linspace(0,90,100),dt_t[d,:])
    
    
    for i in range(depth):
        if par=="EVEN":
            clic_xi_r[i,100::] = clic_xi_r[i,::-1][100::]
            clic_xi_t[i,100::] = clic_xi_t[i,::-1][100::]
            clic_dt_t[i,100::] = clic_dt_t[i,::-1][100::]
        else:
            clic_xi_r[i,100::] = -1.*clic_xi_r[i,::-1][100::]
            clic_xi_t[i,100::] = -1.*clic_xi_t[i,::-1][100::]
            clic_dt_t[i,100::] = -1.*clic_dt_t[i,::-1][100::]
            

    #plt.plot(clic_theta[0:200],clic_dt_t[-2,0:200])

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
        pert_r[:,i] = s_model_c[:,1]*(1.+inv_f*clic_xi_r[-i,:])
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
        
#    plt.plot(clic_theta[0:200],pmodels_fine[2][:,2])
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
def run_pipeline(sub_dir="",minL=False):
    global start,model,vel,modes,old_modes,modeloc,pmodels,pmodels_fine,modes_not_found
    #MODE info:
    par = "ODD"
    model = ["2"]
    vel =  8 #index, not velocity!
    #modes = list_gen([0,1,2,3,4,5,6],1,3,"p")
    modes = list_gen([0,1,2,3,4,5,6,7,8,9,10],1,1,"p")
    #["1 p1","1 p2","1 p3","1 p4","1 p5", "1 p6", "1 p7", "1 p8", "1 p9"]
    #["0 1H","0 2H","0 3H","0 4H","0 5H","0 6H","0 7H","0 8H","0 9H","0 10H"]
    
    
    
    clic = True
    new_mags = True
    save_perturbation = True
    only_mags = False

    
    incl = [0,10,20,30,40,50,60,70,80,90]
    #incl = [90]
    
    
    reese = False
    force_f = False
    f_freq = 3.16
    
    #### contour plotting option:
    plot_contour = False
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
    if mode_by_name==True:
        for i in range(len(modes)):
            if int(modes[i].split()[0])%2==0:
                par = "EVEN"
            else:
                par = "ODD"
            
            parity.append(par)
            idx_tmp = index_by_name(model,vel,par,modes[i])
            print modes[i],"->",idx_tmp
            modes[i] = idx_tmp
    
            
    start = timer()
    klm = 0
    rem = 0
    
    modes_not_found = []
    for mode in modes:
        if mode == 999:
            modes_not_found.append(old_modes[rem])
            continue
        
        rem += 1
        
        if mode_by_name == True:
            par = parity[klm]
            klm +=1
            
        if force_f==True:
            sigma=f_freq
        else:
            sigma = find_sigma(model,vel,par,mode)
            
            
        #Find the perturbed models!
    
        modeloc,pmodels,pmodels_fine = run_visc_pert(model,vel,mode,par,sigma,reese,force_f,minL)
        modeloc += sub_dir
        if not os.path.exists(modeloc):
            os.makedirs(modeloc)
        
        nrzone = len(dt_t[:,0])    
        
        if plot_contour == True:
            make_contour(r,-zp*cs,model,vel,par,mode,excl)
            
        if save_perturbation==True:
            save_pert(modeloc,model,vel,par,mode,old_modes[modes.index(mode)])
        
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
                    
                    
                    pyclic.run_CLIC("model_MODE_"+str(mode)+"_r"+str(rzone),[i],False,3000.,12099.,2.0,par,dt_grid,fine_model=True)
                    
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
#plt.plot(xi_r_rot[-1,:],label=r"$\xi_{r}$")
#plt.plot(xi_t_rot[0,:])
#plt.plot(dt_t_rot[-1,:],label=r"$dT/T$")
#plt.legend(loc="best")
#plt.grid()
#freqs = [1.58717,2.05922,2.49359,2.95717,3.46299,3.99529,4.54267,5.09092,5.64618] # l=0, M2p5 V=0
#modes = [51, 69, 78, 85, 92, 100, 106, 114, 122] # l=0, M2p5 V=0

#freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415,4.47259,5.01979] # l=2, M2p5 V=0
#modes = [8,34,52,66,76,84,90]

#freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V=0
#modes = [12,26,40,49,61,72,80]


#freqs = [1.11169,1.30506,1.55074,1.63928,1.97688,2.25056,2.82636,3.39865] # l=6 M2p5 V=0
#modes = [31, 42, 50, 55, 67, 75, 83, 91]

#freqs = [0.78296,1.63480,2.16027,2.66899,3.18255,3.70941] #l=1 M2p5 V=0
#modes = [5, 61,81,89,96,104] #l=1 M2p5 V0

#freqs = [1.08830,1.42856,1.68216,2.06952,2.54134,3.04737,3.57597] #l=3 M2p5 V=0
#modes = [34, 52, 64, 78, 87, 94, 102]

#freqs = [0.85367,0.99406,1.18376,1.43741,1.59521,1.90170,2.18703,2.73819,3.29288,3.85415] # l=5 M2p5 V=0
#freqs = [1.43741,1.59521,1.90170,2.18703,2.73819]

run_pipeline(sub_dir="")
run_pipeline(sub_dir="minL/",minL=True)


end = timer()
print modes_not_found
print "time:",end-start

