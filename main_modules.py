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
    """ Mode class. Contains relevant information of a mode.
    Parameters
    ----------
    mass : array_like
        data times (one dimensional, length N)
    vel : array_like (optional)
        data values
    mode : string name, integer index or exact frequency float
    
    par : mode parity for mode with no name
    
    """
    
    def __init__(self,mass,vel,mode,par="EVEN",lpert=False):
        model = [mass]
        mode_by_name = False
        self.mass = mass
        self.vel_idx = vel
        
        if type(mode) == float:
            self.index = find_index(freqs,model,self.vel_idx,par)
            mode = self.index

        if type(mode) == int:
            self.name = find_name(model,self.vel_idx,par,mode).strip()
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
            
            self.index = index_by_name(model,self.vel_idx,par,mode)
            if self.index != 999:
                self.eq_vel = find_vel(model,self.vel_idx)
                self.parity = par
                self.frequency = find_sigma(model,vel,par,self.index)
                self.smix = find_smix(model,vel,par,self.index)
                self.temp_surf = load_temp_surf([self.mass],self.vel_idx,self.parity,self.index)
                #self.area = self.calc_area()
                
                if lpert:
                    self.__check_vis()
                    self.__check_vis_phase()
                    self.__load_perturbations()
                    
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
        plt.show()
        #plt.title(self.mass.replace("p",".")+"M$_{\odot}$-V="+self.eq_vel)
        return
        
    def __load_perturbations(self,invert=False,f_scale=True,force_f=False):
        self.xi_r,self.xi_t,self.dt_t,self.zg,self.r,self.zp,self.sig,self.cs,self.where = load_MODE_file([self.mass],self.vel_idx,self.index,self.parity,self.frequency,f_scale,force_f)
        return
        
    def find_mode_location(self,subdir=""):
        static_m = homedir+"ROTORCmodels/visibilities/"
        where =static_m+self.mass+'Msun/V'+self.eq_vel+"/MODE_"+self.parity+"_"+str(self.index)+subdir
        return where
        
    def calc_area(self,nbfs=8,incl=0,weight=1.):
        import visible_disk as vdisk
        rdef = 10.
        rzone = 490
        modeln = "/fine_model_MODE_"+str(self.index)+"_r"+str(rzone)
        self.cos_xi = vdisk.find_cosxi(self.find_mode_location()+modeln,[incl],self.parity,fine_model=True)
        ncx,weight = np.zeros(300),np.zeros(300)
        ncx[0:100] = (self.cos_xi[199,0:100])[::-1]      
        ncx[100::] = self.cos_xi[0,:]
        weight = 1.*ncx
        weight[weight<0] = 0        
        
        xarr = np.linspace(10,80,nbfs)
        
        
        if self.parity=="EVEN":
            newdata = 1.*self.temp_surf
            leg = lint.pyL.legendre(newdata,nbfs)
        else:
            newdata = np.cos(np.deg2rad(xarr))*self.temp_surf
            leg = lint.pyL.legendre_odd(newdata,nbfs)
            
        tleg = np.zeros((3*len(leg[:,0]),2))

        tleg[0:len(leg[:,0]),0] = -1.*leg[::-1,0]
        tleg[0:len(leg[:,0]),1] = leg[::-1,1]
        tleg[len(leg[:,0]):2*len(leg[:,0]),:] = 1.*leg[:,:]
        tleg[2*len(leg[:,0])::,0] = leg[:,0]+90.
        if self.parity=="EVEN":
            tleg[2*len(leg[:,0])::,1] = leg[::-1,1]
        else:
            tleg[2*len(leg[:,0])::,1] = -1*leg[::-1,1]
        
        tlegmin = 0.01*tleg
        tlegmin[:,1] *= weight
        tlegmin[:,1] = rdef - tlegmin[:,1]
        
        tlegold = 1.*tleg
        tlegold[:,1] += rdef
        tleg[:,1] *= 0.01*weight
        tleg[:,1] += rdef
        
        amax = scint.trapz(tleg[:,1],tleg[:,0])
        amin = scint.trapz(tlegmin[:,1],tlegmin[:,0])
        
        if self.parity=="EVEN":
#            leg_plus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]>0 and leg[i,0]>=10 and leg[i,0]<=80])
#            leg_minus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]<0 and leg[i,0]>=10 and leg[i,0]<=80])
            leg_plus = 1.*tleg
            leg_plus[:,1] = 0.
            leg_minus = np.copy(leg_plus)
            for i in range(len(leg_plus[:,0])):
                if leg_plus[i,0]>=-90. and leg_plus[i,0]<=180.:
                    if tleg[i,1]>0:
                        leg_plus[i,1]+=tleg[i,1]
                    else:
                        leg_minus[i,1]+=tleg[i,1]
                    
            #leg_plus = np.array([tleg[i] for i in range(len(tleg[:,0])) if tleg[i,1]>0 and tleg[i,0]>=-90.+incl and tleg[i,0]<=90.+incl])
            #leg_minus = np.array([tleg[i] for i in range(len(tleg[:,0])) if tleg[i,1]<0 and tleg[i,0]>=-90.+incl and tleg[i,0]<=90.+incl])
        else:
#            leg_plus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]>0 and leg[i,0]>=10 and leg[i,0]<=80])
#            leg_minus = np.array([leg[i] for i in range(len(leg[:,0])) if leg[i,1]<0 and leg[i,0]>=10 and leg[i,0]<=80])
            leg_plus = 1.*tleg
            leg_plus[:,1] = 0.
            leg_minus = np.copy(leg_plus)
            for i in range(len(leg_plus[:,0])):
                if leg_plus[i,0]>=-90. and leg_plus[i,0]<=180.:
                    if tleg[i,1]<0:
                        leg_plus[i,1]+=tleg[i,1]
                    else:
                        leg_minus[i,1]+=tleg[i,1]
#            leg_plus = np.array([tleg[i] for i in range(len(tleg[:,0])) if tleg[i,1]>0 and tleg[i,0]>=-90.+incl and tleg[i,0]<=90.+incl])
#            leg_minus = np.array([tleg[i] for i in range(len(tleg[:,0])) if tleg[i,1]<0 and tleg[i,0]>=-90.+incl and tleg[i,0]<=90.+incl])
#        
#
#        plt.plot(tleg[:,0],tleg[:,1],label="Weighted by cos(xi)")
#        plt.plot(tlegmin[:,0],tlegmin[:,1],label="Weighted by cos(xi)")
#        plt.plot(tlegold[:,0],tlegold[:,1],"--",alpha=.5,color="k",label="Original perturbation")
#        #plt.fill_between(leg_plus[:,0],0,leg_plus[:,1],alpha=0.5,color="blue")
#        #plt.fill_between(leg_minus[:,0],0,leg_minus[:,1],alpha=0.5,color="red")
#        yold = plt.ylim()
#        plt.vlines(90,yold[0],yold[1],lw=1.)
#        plt.vlines(0,yold[0],yold[1],lw=1.)
#        plt.ylim(yold[0],yold[1])
#        plt.grid()
#        plt.xlim(-90,180)
#        plt.xlabel("Colatitude")
#        plt.ylabel(r"$\xi_{r}$")
#        plt.title(self.name+" - i = "+str(incl)+r"$^{\circ}$")
        return abs(np.log10(amax/amin))
        
        if len(leg_plus)==0: leg_plus=np.zeros([2,2])
        if len(leg_minus)==0: leg_minus=np.zeros([2,2])
        a_plustemp = np.abs(scint.trapz(leg_plus[:,1],leg_plus[:,0]))
        a_minustemp = np.abs(scint.trapz(leg_minus[:,1],leg_minus[:,0]))
        if a_plustemp > a_minustemp:
            return a_minustemp/a_plustemp
        else:
            return a_plustemp/a_minustemp
        return np.abs(scint.trapz(tleg[:,1],tleg[:,0]))
            
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
        
    def __check_vis_phase(self):
        static_m = homedir+"ROTORCmodels/visibilities/"
        where =static_m+self.mass+'Msun/V'+self.eq_vel+"/MODE_"+self.parity+"_"+str(self.index)
        lst = sorted(glob.glob(where+"/phi_max/magnitudes*"))
        lst_min = sorted(glob.glob(where+"/phi_min/magnitudes*"))
        lst_wal = sorted(glob.glob(where+"/phi_max/walraven_magnitudes*"))
        lst_wal_min = sorted(glob.glob(where+"/phi_min/walraven_magnitudes*"))
        self.j_phi_mags = []
        self.j_phi_mags_min = []
        self.w_phi_mags = []
        self.w_phi_mags_min = []
        if lst > 0:
            self.j_phi_mags = np.array([np.genfromtxt(lst[i]) for i in range(len(lst))])
            self.j_phi_mags_min = np.array([np.genfromtxt(lst_min[i]) for i in range(len(lst_min))])
        else:
            print "No Johnson magnitudes available!"
            
        if lst_wal > 0:
            self.w_phi_mags = np.array([np.genfromtxt(lst_wal[i]) for i in range(len(lst))])
            self.w_phi_mags_min = np.array([np.genfromtxt(lst_wal_min[i]) for i in range(len(lst_wal_min))])
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
                    
            if tmpmodes[0]==tmpmode[0] and tmpmode[1].find("f") !=-1:
                if tmpmodes[1].find("f") !=-1 or tmpmodes[1].find("F") !=-1:
                    idx = i+1
                    
            elif tmpmodes[0]==tmpmode[0] and tmpmode[1].find("F") !=-1:
                
                if tmpmodes[1].find("f") !=-1 or tmpmodes[1].find("F") !=-1:
                    idx = i+1

                        
    if idx==999:
        print "MODE",mode,"NOT FOUND BY NAME!"
        
    return idx


    
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
    
def load_MODE_file(model,vel,mode,par,sigma,f_scale,force_f):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    norm_f = True
    if f_scale:
        scale = 0.02/(sigma**2)
    else:
        scale = 0.02
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
    
    #global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    #global xi_r,xi_t,dt_t,zg,r,zp,cs
    #global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    xi_r,xi_t,dt_t,zg,r,zp,sig,cs,xi_dot_g = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname)
            
    #xi_r_n,xi_t_n,dt_t_n,zg_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,zg,norm_f,scale,depth,f_scale,sig,par)

#    a_r = scint.trapz(dt_t_n[-1,:])
#    
#    inv_f = 1.
#    if a_r<0:
#        print "-T area"
#
#        xi_r,xi_t,dt_t = -1.*xi_r,-1.*xi_t,-1.*dt_t
        
    
    
    return xi_r,xi_t,dt_t,zg,r,zp,sig,cs,where
    
def phaser(var,phase,amplitude,tlmax):
    new_var = amplitude*np.abs(np.cos(tlmax + np.deg2rad(phase)))*var
    print "R_cos -- -- ",np.abs(np.cos(tlmax + np.deg2rad(phase)))
    return new_var
    
def phaser2(var,tphase,tampl,rampl):
    import sympy as spy
    from scipy.interpolate import interp1d
    #new_var = amplitude*np.cos(tlmax + np.deg2rad(phase))*var
    theta = np.linspace(0,360,500)
    phases = np.deg2rad(np.linspace(-80,0,500))
    a1 = tampl
    a2 = rampl
    print "Tampl =",a1," -- Rampl =",a2    
    f = []
    for phi in phases:
    #phi = np.pi + sfphi.val   
        
        aa = 4.*a1
        bb = 2.*a2
        ycomp = bb*np.sin(0)+aa*np.sin(np.pi+phi)
        xcomp = bb*np.cos(0)+aa*np.cos(np.pi+phi)
        phi2 = np.rad2deg(np.arctan(ycomp/xcomp))

        psi_T = np.rad2deg(phi+np.pi)
        f.append([phi2+np.rad2deg(phi+np.pi),psi_T])
        #ampl = np.sqrt((bb*np.cos(0)+aa*np.cos(phi))**2+(bb*np.sin(0)+aa*np.sin(phi))**2)
        #ndl = ampl*np.cos(t+phi2+phi)
        
    f = np.array(f)
    fx = interp1d(f[:,0],f[:,1])
    
    psi_T = fx(tphase)
    print "Psi_T = ",psi_T,"deg","Psi_L = ",tphase,"deg"
    x = spy.Symbol('x',real=True)
    f = spy.cos(x+np.deg2rad(tphase))
    fprime = f.diff(x)
    mx = spy.solve(fprime,x)[-1]
    mx = float(mx)
#    plt.plot(theta,np.cos(np.deg2rad(theta)+np.deg2rad(psi_T)))
#    plt.plot(theta,np.cos(np.deg2rad(theta)+np.deg2rad(tphase)))
#    plt.plot(theta,np.cos(np.deg2rad(theta)))
#    plt.vlines(np.rad2deg(float(mx))+180.,-1.2,1.2)
#    plt.vlines(np.rad2deg(float(mx)),-1.2,1.2)
#    plt.ylim(-1.2,1.2)
#    plt.xlim(0,360)
    print "T_cos -- -- ", np.abs(np.cos(mx + np.deg2rad(psi_T)))
      
    return np.abs(np.cos(mx + np.deg2rad(psi_T)))*var,psi_T,phi2+np.rad2deg(phi+np.pi),mx
    
def load_mags(subdir,mass,vel,par,mode,incl):
    v = find_vel([mass],vel)
    
    static_m = homedir+"ROTORCmodels/visibilities/"
    
    where =static_m+mass+'Msun/V'+v+"/MODE_"+par+"_"+str(mode)+subdir+"/"
    
    
    globwal = np.sort((glob.glob(where+"walraven_magnitudes_i*")))
    
    globjh = np.sort(glob.glob(where+"magnitudes_i*"))
    
    print where, len(globwal), len(globjh)
    
    wal = []
    jh = []
    for i in incl:
        try:
            wal.append(np.genfromtxt(where+"walraven_magnitudes_i"+str(i)+"_MODE_"+str(mode)))
        except:
            print "MISSING WALRAV?",mass,str(v),par,str(mode),"i =",i
            
    for i in incl:
        try:
            jh.append(np.genfromtxt(where+"magnitudes_i"+str(i)+"_MODE_"+str(mode)))
        except:
            print "MISSING Johnson?",mass,str(v),par,str(mode),"i =",i     

    
    return wal,jh

def find_dr_ofmax(mode,depth,vels):
    maxvals = []
    for i in range(len(vels)):
        tmode = emode(mode.mass,i,mode.name)
        xi_r,xi_t,dt_t,zg,r,zp,sig,cs,where = load_MODE_file([mode.mass],i,tmode.index,tmode.parity,tmode.frequency,True,False)
    
        if mode.parity == "EVEN":
            xi_r_fine = lint.leg_interp(xi_r[-depth::,:],8,"EVEN")
            xi_t_fine = lint.leg_interp(xi_t[-depth::,:],8,"EVEN")
            dt_t_fine = lint.leg_interp(dt_t[-depth::,:],8,"EVEN")
        else:
            xi_r_fine = lint.leg_interp(xi_r[-depth::,:],8,"OE")
            xi_t_fine = lint.leg_interp(xi_t[-depth::,:],8,"OE")
            dt_t_fine = lint.leg_interp(dt_t[-depth::,:],8,"OE")
            
        rs_n = np.empty(xi_r_fine.shape)
        for d in range(depth):
            rs_n[d,:] = np.interp(np.linspace(0,90,100),np.linspace(10,80,8),r[-d-1,:])    
        
        #plt.plot(np.abs(xi_r_fine[-1,:]*rs_n[-1,:]))
        maxval = np.argmax(np.abs(xi_r_fine[-1,:]*rs_n[-1,:]))
        maxvals.append([maxval,np.linspace(0,90,100)[maxval],(xi_r_fine[-1,:]*rs_n[-1,:])[maxval]])
    maxvals = np.array(maxvals)
    return maxvals