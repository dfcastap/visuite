# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 14:37:11 2016

@author: diego
"""

import numpy as np
import pylab as plt
import glob,subprocess,os
import del_dot_xi as ddxi
import pyNRO_1_mode as pyNRO
import scipy.integrate as scint
import legendre_interp as lint
import os


vis_path = os.getcwd()

global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/castaned/Documents/"
else:
    homedir = "/home/castaned/"
    
class emode:
    """ Mode class. Contains relevant information of a mode."""
    
    def __init__(self,mass,vel,mode,par="EVEN"):
        model = [mass]
        mode_by_name = False
        self.mass = mass
        self.vel_idx = vel
        
        if type(mode) == float:
            self.index = find_index(freqs,model,vel,par)
            mode = self.index

        if type(mode) == int:
            self.name = find_name(model,vel,par,mode).strip()
            mode = self.name
            mode_by_name = True
            
        if type(mode) == str:
            mode_by_name = True
            self.name = mode
            
        
        if mode_by_name==True:
            if int(mode.split()[0])%2==0:
                par = "EVEN"
            else:
                par = "ODD"
            
            self.index = index_by_name(model,vel,par,mode)
            if self.index != 999:
                self.eq_vel = find_vel(model,vel)
                self.parity = par
                self.frequency = find_sigma(model,vel,par,self.index)
                self.smix = find_smix(model,vel,par,self.index)
                self.temp_surf = load_temp_surf([self.mass],self.vel_idx,self.parity,self.index)
                self.area = self.__calc_area()
                self.__check_vis()
            else:
                print "mode not found!"

    
    def plot_dr_r(self,nbfs=8,invert=False):
        import seaborn
        import pylab as plt
        seaborn.set(style="whitegrid",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 10,'axes.titlesize':18,'font.size':14})

        #freq_val = QString(str())
        xarr = np.linspace(10,80,nbfs)
        #data = np.genfromtxt(self.models[0]+"/temp_surf"+str(key))
        if self.parity=="EVEN":
            newdata = 1.*self.temp_surf
            leg = lint.pyL.legendre(newdata,nbfs)
        else:
            newdata = np.cos(np.deg2rad(xarr))*self.temp_surf
            leg = lint.pyL.legendre_odd(newdata,nbfs)
        
        plt.xlabel("Polar Angle (deg)")
        plt.ylabel("$\delta$R/R")
        plt.xlim(0,90)
        plt.hlines(0,0,90,color='k',linestyle='--')
        if invert==True:
            newdata = -1.*newdata
            leg[:,1] = -1.*leg[:,1]
        plt.plot(xarr,newdata,"o",color="b",mfc="white",mec="k",mew=1)
        plt.plot(leg[:,0],leg[:,1],label=r"$\ell$="+self.name+"-("+self.mass.replace("p",".")+"M$_{\odot}$-V="+self.eq_vel+")")     
        plt.legend(loc="best")
        #plt.title(self.mass.replace("p",".")+"M$_{\odot}$-V="+self.eq_vel)
        return
        
    def __calc_area(self,nbfs=8):
        xarr = np.linspace(10,80,nbfs)
        
        if self.parity=="EVEN":
            newdata = 1.*self.temp_surf
            leg = lint.pyL.legendre(newdata,nbfs)
        else:
            newdata = np.cos(np.deg2rad(xarr))*self.temp_surf
            leg = lint.pyL.legendre_odd(newdata,nbfs)
            

        if self.parity=="EVEN":
            leg_plus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]>0 and leg[i,0]>=10 and leg[i,0]<=80])
            leg_minus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]<0 and leg[i,0]>=10 and leg[i,0]<=80])
        else:
            leg_plus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]>0 and leg[i,0]>=10 and leg[i,0]<=80])
            leg_minus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]<0 and leg[i,0]>=10 and leg[i,0]<=80])
        if len(leg_plus)==0: leg_plus=np.zeros([2,2])
        if len(leg_minus)==0: leg_minus=np.zeros([2,2])
        a_plustemp = np.abs(scint.trapz(leg_plus[:,1],leg_plus[:,0]))
        a_minustemp = np.abs(scint.trapz(leg_minus[:,1],leg_minus[:,0]))
        if a_plustemp > a_minustemp:
            return a_minustemp/a_plustemp
        else:
            return a_plustemp/a_minustemp
            
    def __check_vis(self):
        static_m = homedir+"ROTORCmodels/visibilities/"
        where =static_m+self.mass+'Msun/V'+self.eq_vel+"/MODE_"+self.parity+"_"+str(self.index)
        lst = sorted(glob.glob(where+"/magnitudes*"))
        lst_min = sorted(glob.glob(where+"/minL/magnitudes*"))
        lst_wal = sorted(glob.glob(where+"/walraven_magnitudes*"))
        lst_wal_min = sorted(glob.glob(where+"/minL/walraven_magnitudes*"))
        self.j_mags = []
        self.j_mags_min = []
        self.w_mags = []
        self.w_mags_min = []
        if lst > 0:
            self.j_mags = np.array([np.genfromtxt(lst[i]) for i in range(len(lst))])
            self.j_mags_min = np.array([np.genfromtxt(lst_min[i]) for i in range(len(lst_min))])
        else:
            print "No Johnson magnitudes available!"
            
        if lst_wal > 0:
            self.w_mags = np.array([np.genfromtxt(lst_wal[i]) for i in range(len(lst))])
            self.w_mags_min = np.array([np.genfromtxt(lst_wal_min[i]) for i in range(len(lst_wal_min))])
        else:
            print "No Walraven magnitudes available!"
            
        return
        
        
        

def find_vel(model,vel):
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
    
    return vels[vel]

def load_temp_surf(model,vel,par,mode):
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_surf = np.genfromtxt(folder+"/pack_temp_surf")
    
    return temp_surf[mode-1]
   
def find_index(freq,model,vel,par):
    global temp_freqs
    
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
    temp_freqs[:] *= 1e5
    temp_freqs = ((np.array(temp_freqs)).astype(int)).tolist()
    idx = []
    for i in range(len(freq)):
        try:
            idx.append(temp_freqs.index(int(freq[i]*1e5))+1)
        except:
            print "No match for freq: ",freq[i]
            
    return idx
    
def find_sigma(model,vel,par,mode):
    global temp_freqs
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
            
    return temp_freqs[mode-1]
    
def find_smix(model,vel,par,mode): 
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_smix = np.genfromtxt(folder+"/temp_smixes")
    
    return temp_smix[mode-1]
    
def find_name(model,vel,par,mode):
    global temp_modes
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_modes = np.genfromtxt(folder+"temp_MODES")
    try:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype='|S8',delimiter=8)
    except:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype=str)
        temp_modes = [temp_modes[i,0]+" "+temp_modes[i,1] for i in range(len(temp_modes[:,0]))]

            
    return temp_modes[mode-1]
    

def index_by_name(model,vel,par,mode):
    global temp_modes
    v = find_vel(model,vel)
    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_modes = np.genfromtxt(folder+"temp_MODES")
    try:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype='|S8',delimiter=8)
    except:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype=str)
        temp_modes = [temp_modes[i,0]+" "+temp_modes[i,1] for i in range(len(temp_modes[:,0]))]
    
    idx = 999
    for i in range(len(temp_modes)):
        tmpmodes = temp_modes[i].split()
        tmpmode = mode.split()
                
        if len(tmpmodes)>1:
            if tmpmode[0] == "0":
                if tmpmodes[0]==tmpmode[0] and (tmpmodes[1].strip("H")).strip("h")==(tmpmode[1].strip("H")).strip("h").strip("p"):
                    idx = i+1
            else:
                if tmpmodes[0]==tmpmode[0] and tmpmodes[1]==tmpmode[1]:
                    idx = i+1
    
    if idx==999:
        print "MODE",mode,"NOT FOUND BY NAME!"
        
    
    return idx

def run_visc_pert(model,vel,mode,par,sigma,reese,force_f):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    norm_f = True
    scale = 0.125/(sigma**2)
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
    s_model = np.genfromtxt(glob.glob(static_m+model[0]+'Msun/'+model[0]+'Msun_V'+v+"*")[0])
    
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    global xi_r,xi_t,dt_t,zg,r,zp,cs
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    xi_r,xi_t,dt_t,zg,r,zp,sig,cs = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname)
            
    xi_r_n,xi_t_n,dt_t_n,zg_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,zg,norm_f,scale,depth,reese,sig)

    a_r = scint.trapz(dt_t_n[-1,:])
    
    inv_f = 1.
    if a_r<0:
        print "-T area"
        #xi_r_n *= -1.
        #xi_t_n *= -1.
        #dt_t_n *= -1.
        #zg_n *= -1.
        inv_f = -1.

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
        
    

    
#    plt.plot(p_angles_fine[2,:],xi_r_fine[2,:])
#    plt.plot(p_angles[2,:],xi_r_n[2,:],"o",mfc="white",mec="k",mew=1)
#    plt.grid()
#    plt.xlim(0,90)
#    plt.yticks()
#    plt.xlabel("Colatitude [deg]")
#    plt.ylabel(r"$\xi_{r}$")
#    plt.title("Mode:" + find_name(model,vel,par,mode))
#    #ax.yaxis.set_major_formatter(majorFormatter) 
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
       
    
    print "Generating fine clic model..."    
    pmodels_fine = gen_fine_clic(p_angles_fine,xi_r_fine,xi_t_fine,dt_t_fine,s_model,model,depth,inv_f)

        
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
    
    
#    import seaborn as sns
#    sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
#                              'ytick.labelsize': 12,'xtick.labelsize': 12,
#                              'legend.fontsize': 16,'axes.titlesize':18,'font.size':14})
#    #plt.plot(pmodels[2][0:-1,0],np.diff(pmodels[2][:,1])/(pmodels[2][1,0]-pmodels[2][0,0]))
#    plt.plot(pmodels[2][:,0],(pmodels[2][:,2]),lw=1.5,marker="o")
#    #plt.plot(pmodels_fine[2][0:-1,0],np.diff(pmodels_fine[2][:,1])/(pmodels_fine[2][1,0]-pmodels_fine[2][0,0]))
#    plt.plot(pmodels_fine[2][:,0],pmodels_fine[2][:,2],lw=1.5)
#    plt.plot(s_model[:,0],s_model[:,2])
#    #plt.vlines(90,9100,9300)
#    o_lims = plt.ylim()
#    plt.vlines(90,o_lims[0],o_lims[1])
#    plt.grid()
#    plt.xlim(0,180)
#    plt.ylim(o_lims[0],o_lims[1])
#    plt.xlabel("Colatitude [deg]")
#    #plt.ylabel("Radius [M$_{\odot}$]")
#    plt.ylabel("Temperature [K]")
#    plt.title("Mode: " + find_name(model,vel,par,mode))
    
    
    return modeloc,pmodels,pmodels_fine
    
def gen_fine_clic(p_angles_fine,xi_r,xi_t,dt_t,s_model,model,depth,inv_f):
    clic_theta_grid = 200
    clic_theta = np.linspace(0,180,clic_theta_grid)
    global clic_xi_r,clic_xi_t,clic_dt_t
    
    clic_xi_r = np.zeros((depth,clic_theta_grid))
    clic_xi_t = np.zeros((depth,clic_theta_grid))
    clic_dt_t = np.zeros((depth,clic_theta_grid))
    
    for d in range(depth):
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
            
    
    #plt.plot(clic_theta[0:200],clic_xi_r[-2,0:200])

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
    
def save_pert(model,vel,par,mode,l,n):
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    global xi_r,xi_t,dt_t,zg,r
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    print os.getcwd()
    f = open("M"+model[0]+"_V"+str(vel)+"_l"+l+"_n"+n+"_perturbations","w")  
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
    
