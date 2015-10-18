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
vis_path = os.getcwd()

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
    
def find_index(freq,model,vel,par):
    global temp_freqs
    
    v = find_vel(model,vel)
    folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
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
    folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
            
    return temp_freqs[mode-1]
    
def find_name(model,vel,par,mode):
    global temp_modes
    v = find_vel(model,vel)
    folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    temp_modes = np.genfromtxt(folder+"temp_MODES")
    try:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype='|S8',delimiter=8)
    except:
        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype=str)
        temp_modes = [temp_modes[i,0]+" "+temp_modes[i,1] for i in range(len(temp_modes[:,0]))]

            
    return temp_modes[mode-1]

def run_visc_pert(model,vel,mode,par,sigma):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    norm_f = True
    scale = 0.005
    depth = 10 #radial zones down from the surface
    #---------------------------------------------------
    
    clic_run = True
    global zp
    
    
    v = find_vel(model,vel)
    folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
    
    rotorc_f = "/home/diego/Documents/From_Bob/Delta_Scuti_2010/"+model[0]+"Msun/"+model[0]+"Msun_V"+v+"/"
    
    bob_bin = "/home/diego/Documents/From_Bob/clotho_disc10_bin/"
    
    static_m = "/home/diego/Documents/ROTORCmodels/visibilities/"
    
    
    
    
    #check if visibility file from pulset exits for the selected model:
    
    v_file = glob.glob(rotorc_f+"visibility_file")
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
    
    print "v_file OK!"
    
    ###### FULL MODE FILE GENERATION:
    
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
    tfreq = temp_freqs[mode-1]
    #tfreq = 1.59692
    print pyNRO.run_nro(tfreq,folder,model[0],v,par,mode)
    
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+v):
        os.makedirs(static_m+model[0]+'Msun/'+'V'+v)
        
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+v+"/MODE_"+par+"_"+str(mode)):
        os.makedirs(static_m+model[0]+'Msun/'+'V'+v+"/MODE_"+par+"_"+str(mode))
    
    try:
        where =static_m+model[0]+'Msun/V'+v+"/MODE_"+par+"_"+str(mode)
        subprocess.call(['mv',folder+'MODE_'+par+'_'+str(mode),where+'/MODE_'+par+'_'+str(mode)])
        print "Mode file generation complete!"
    except:
        print "Something broke! Couldn't copy the Full mode file..."
    
    
    ###### del_DOT_xi>
    #modeloc = static_m+model[0]+'Msun/V'+vels[vel]+"/"
    modeloc = where+"/"
    modefname = 'MODE_'+par+'_'+str(mode)
    #modefname = 'MODE13'
    s_model = np.genfromtxt(glob.glob(static_m+model[0]+'Msun/'+model[0]+'Msun_V'+v+"*")[0])
    
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    global xi_r,xi_t,dt_t,zg,r
    global xi_r_n,xi_t_n,dt_t_n,zg_n
    
    xi_r,xi_t,dt_t,zg,r,zp = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname,sigma)
            
    xi_r_n,xi_t_n,dt_t_n,zg_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,zg,norm_f,scale,depth)

    a_r = scint.trapz(xi_r_n[-1,:])
    
    if a_r>0:
        xi_r_n *= -1
        xi_t_n *= -1
        dt_t_n *= -1
        zg_n *= -1
    
    xi_r_rot,xi_t_rot,dt_t_rot,zg_rot = ddxi.to_rotorc(xi_r_n,xi_t_n,dt_t_n,zg_n)
    
        
    ############ find the perturbed models:
    pert_r = np.empty((len(s_model[:,0]),depth))
    pert_t = np.empty((len(s_model[:,0]),depth))
    pmodels = []
    for i in range(depth):
        pert_r[:,i] = s_model[:,1]*(1.+xi_r_rot[-i,:])
        pert_t[:,i] = s_model[:,2]*(1.+dt_t_rot[-i,:])
        tmodel = 1.*s_model
        tmodel[:,1] = pert_r[:,i]
        tmodel[:,2] = pert_t[:,i]
        pmodels.append(tmodel)

    if par=="ODD":
        pert_r = np.empty((2*len(s_model[:,0]),depth))
        pert_t = np.empty((2*len(s_model[:,0]),depth))
        pmodels = []
        odd_angles = np.arange(4.5,180.,9)
        for i in range(depth):
            pert_r[0:10,i] = s_model[:,1]*(1.+xi_r_rot[-i,:])
            pert_t[0:10,i] = s_model[:,2]*(1.+dt_t_rot[-i,:])
            pert_r[10:20,i] = s_model[::-1,1]*(1.-xi_r_rot[-i,::-1])
            pert_t[10:20,i] = s_model[::-1,2]*(1.-dt_t_rot[-i,::-1])
            tmodel = np.empty((20,5))
            tmodel[:,0] = odd_angles[:]
            tmodel[:,1] = pert_r[:,i]
            tmodel[:,2] = pert_t[:,i]
            tmodel[0:10,3] = s_model[:,3]
            tmodel[10:20,3] = s_model[::-1,3]
            tmodel[0:10,4] = s_model[:,4]
            tmodel[10:20,4] =s_model[::-1,4]
            
            pmodels.append(tmodel)
    
    return modeloc,pmodels
    
def make_contour(r,var,model,vel,par,mode):
    import seaborn as sns
    sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,'axes.titlesize':18})
    new_r = lint.leg_interp(r,8,"EVEN")
    new_zp = lint.leg_interp(var,8,"EVEN")
    #levels = np.linspace(np.min(new_zp),np.max(new_zp), 40)
    theta = np.linspace(0,np.deg2rad(90),100)
    newt,n_r = np.meshgrid(theta,new_r[:,0])
    

    CS = plt.contourf((new_r*np.sin(newt)),(new_r*np.cos(newt)),new_zp, 100, cmap=plt.cm.jet,vmax=np.max(new_zp), vmin=np.min(new_zp))
    CSl = plt.contour((new_r*np.sin(newt)),(new_r*np.cos(newt)),new_zp, 20, colors="k")
    #CS = plt.contourf((new_r*np.cos(newt)), (new_r*np.sin(newt)), new_zp, cmap=plt.cm.Spectral,levels=levels)
    #plt.axes().set_aspect('equal')
    #plt.xlim(np.ceil(r[-1,-1]))
    plt.ylim(plt.xlim())
    plt.xlabel("Radius [R$_{\odot}$]")
    plt.ylabel("Radius [R$_{\odot}$]")
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
    
    f = open("M2p5_V0_l"+l+"_n"+n+"_perturbations","w")  
    f.write("M2p5, V=0, Mode freq = %.5f\n" % find_sigma(model,vel,par,mode))
    f.write("n_angle, r, xi_r, xi_t, dT/T\n")
    for i in range(len(r[:,0])-10,len(r[:,0])):
        for j in range(len(xi_r[-1,:])):
            #print "%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j])
            f.write("%i %8.5f %8.5f %8.5f %8.5f\n"%(j+1,r[i,j],xi_r[i,j],xi_t[i,j],dt_t[i,j])) 
            
    f.close()
    return
            
"""

------------------------------------------------------
------------------------------------------------------

"""
#MODE info:
par = "temp"
model = ["2p5"]
vel = 0 #index, not velocity!
modes = [3]
mode_by_freq = False
clic = False
new_mags = False
save_perturbation = False

plot_contour = False
#freqs = [1.58717,2.05922,2.49359,2.95717,3.46299,3.99529,4.54267,5.09092,5.64618] # l=0, M2p5 V=0
#modes = [51, 69, 78, 85, 92, 100, 106, 114, 122] # l=0, M2p5 V=0

#freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415,4.47259,5.01979] # l=2, M2p5 V=0
#modes = [8,34,52,66,76,84]

#freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V=0
#modes = [12,26,40,49,61,72,80]

#freqs = [1.97688,2.25056,2.82636,3.39865,3.97414,4.55198,5.12160] # l=6 M2p5 V=0
#modes = [67,75,83,91]

#freqs = [0.78296,1.63480,2.16027,2.66899,3.18255,3.70941] #l=1 M2p5 V=0
#modes = [5, 61,81,89,96,104] #l=1 M2p5 V0

#freqs = [1.08830,1.42856,1.68216,2.06952,2.54134,3.04737,3.57597] #l=3 M2p5 V=0


#freqs = [1.43741,1.59521,1.90170,2.18703,2.73819,3.29288,3.85415]

dt_grid = 500

if mode_by_freq==True:
    modes = find_index(freqs,model,vel,par)
    

for mode in modes:
    sigma = find_sigma(model,vel,par,mode)
    #Find the perturbed models!
    modeloc,pmodels = run_visc_pert(model,vel,mode,par,sigma)
    if plot_contour == True:
        make_contour(r,xi_r,model,vel,par,mode)
        
    if save_perturbation==True:
        save_pert(model,vel,par,mode,"1","1")
    
    for i in range(len(pmodels)):
        np.savetxt(modeloc+"model_MODE_"+str(mode)+"_r"+str(i+1),pmodels[i],'%.3f')
    
    if new_mags==True:
        try:
            subprocess.call(['rm',modeloc+"magnitudes_MODE_"+str(mode)])
        except:
            print "No mag file to erase!"    
    
    if clic==True:
        #Ugly way of running pyCLIC:
        clic_folder = "/home/diego/Documents/pyCLIC/test/"
        incl = [90.0]
        rzone = 1
        
        for i in incl:
            subprocess.call(['cp',modeloc+"model_MODE_"+str(mode)+"_r"+str(rzone),clic_folder])
                
            if (os.path.isfile(modeloc+"magnitudes_MODE_"+str(mode)))==False:
                myfile = open(modeloc+"magnitudes_MODE_"+str(mode), "w")
                myfile.write('#%5s %5s %9s %11s %11s %11s %11s %11s\n' % ("incl","n_r","lum","u","b","v","r","i"))
                myfile.close()
                subprocess.call(['cp',modeloc+"magnitudes_MODE_"+str(mode),clic_folder])
            else:
                subprocess.call(['cp',modeloc+"magnitudes_MODE_"+str(mode),clic_folder])
            
            os.chdir(clic_folder)
            import main as pyclic
            import color_magnitudes as cmag
            
            pyclic.run_CLIC("model_MODE_"+str(mode)+"_r"+str(rzone),[i],False,3000.,12099.,2.0,par,dt_grid)
            
            cmag.calc_mags('outputflux_i'+str(i)+'.final',[i],mode,rzone)
            subprocess.call(['cp','outputflux_i'+str(i)+'.final',modeloc+"outputflux_i"+str(i)+"_MODE_"+str(mode)+"_r"+str(rzone)])
            subprocess.call(['rm',"model_MODE_"+str(mode)+"_r"+str(rzone)])
            subprocess.call(['rm','outputflux_i'+str(i)+'.final'])
            subprocess.call(['cp',"magnitudes_MODE_"+str(mode),modeloc])
            subprocess.call(['rm',"magnitudes_MODE_"+str(mode)])
            os.chdir(vis_path)

#plt.plot(xi_r_rot[-1,:],label=r"$\xi_{r}$")
#plt.plot(xi_t_rot[0,:])
#plt.plot(dt_t_rot[-1,:],label=r"$dT/T$")
#plt.legend(loc="best")
#plt.grid()
