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
import pandas as pd
vis_path = os.getcwd()

def find_index(freq,model,vel,par):
    global temp_freqs
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
    

def run_visc_pert(model,vel,mode,par):
    # Info for del dot xi calculation:------------------
    # NRO Mode filename:
    #modefname="MODE1"
    norm_f = True
    scale = 0.025
    depth = 10 #radial zones down from the surface
    #---------------------------------------------------
    
    clic_run = True
    
    
    
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
        subprocess.call(["cp","Cmod",folder+"Dmod_"+model[0]+"M_V"+vels[vel]])
        print "Generated new Dmod"
        #print(glob.glob("*_ZAMS"))
        os.chdir(vis_path)
    
    print "v_file OK!"
    
    ###### FULL MODE FILE GENERATION:
    
    temp_freqs = np.genfromtxt(folder+"temp_freqs")
    tfreq = temp_freqs[mode-1]
    #tfreq = 1.59692
    print pyNRO.run_nro(tfreq,folder,model[0],vels[vel],par,mode)
    
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+vels[vel]):
        os.makedirs(static_m+model[0]+'Msun/'+'V'+vels[vel])
        
    if not os.path.exists(static_m+model[0]+'Msun/'+'V'+vels[vel]+"/MODE_"+par+"_"+str(mode)):
        os.makedirs(static_m+model[0]+'Msun/'+'V'+vels[vel]+"/MODE_"+par+"_"+str(mode))
    
    try:
        where =static_m+model[0]+'Msun/V'+vels[vel]+"/MODE_"+par+"_"+str(mode)
        subprocess.call(['mv',folder+'MODE_'+par+'_'+str(mode),where+'/MODE_'+par+'_'+str(mode)])
        print "Mode file generation complete!"
    except:
        print "Something broke! Couldn't copy the Full mode file..."
    
    
    ###### del_DOT_xi>
    #modeloc = static_m+model[0]+'Msun/V'+vels[vel]+"/"
    modeloc = where+"/"
    modefname = 'MODE_'+par+'_'+str(mode)
    #modefname = 'MODE13'
    s_model = np.genfromtxt(glob.glob(static_m+model[0]+'Msun/'+model[0]+'Msun_V'+vels[vel]+"*")[0])
    
    global xi_r_rot,xi_t_rot,dt_t_rot,zg_rot    
    
    xi_r,xi_t,dt_t,zg,r = ddxi.calcdeldotxi(par,model,vel,modeloc,modefname)
            
    xi_r_n,xi_t_n,dt_t_n,zg_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,zg,norm_f,scale,depth)
    
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
        
    
    return modeloc,pmodels
"""

------------------------------------------------------
------------------------------------------------------

"""
#MODE info:
par = "EVEN"
model = ["2p5"]
vel = 0 #index, not velocity!
modes = [61]
mode_by_freq = True
#freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705]
#freqs = [1.58717,2.05922,2.49359,2.95717,3.46299,3.99529,4.54267,5.09092,5.64618] # l=0, M2p5 V=0
# indexes = [51, 69, 78, 85, 92, 100, 106, 114, 122] # l=0, M2p5 V=0

#freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415,4.47259,5.01979] # l=2, M2p5 V=0
# indexes = [9,34,52,66,76,84]

freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V=0
# indexes = 

if mode_by_freq==True:
    modes = find_index(freqs,model,vel,par)
    

for mode in modes:
    #Find the perturbed models!
    modeloc,pmodels = run_visc_pert(model,vel,mode,par)
    
    for i in range(len(pmodels)):
        np.savetxt(modeloc+"model_MODE_"+str(mode)+"_r"+str(i+1),pmodels[i],'%.3f')
        
    #Ugly way of running pyCLIC:
    clic_folder = "/home/diego/Documents/pyCLIC/test/"
    incl = [0.0,50.0,90.0]
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
        
        pyclic.run_CLIC("model_MODE_"+str(mode)+"_r"+str(rzone),[i],False,3000.,12099.,2.0)
        
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
