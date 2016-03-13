# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:39:55 2015

@author: diego
"""

import numpy as np
import pylab as plt
import glob,os,subprocess
import seaborn as sns
import color_magnitudes as cmag
from main_modules import find_vel, find_sigma, index_by_name

dt_grid = 250
vis_path = os.getcwd()

global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/castaned/Documents/"
else:
    homedir = "/home/castaned/"
    
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 16,'axes.titlesize':18,'font.size':14})
wavebands = [3650.,4450.,5510.,6580.,8060.]

wal_wave = [3250,3630,3840,4320,5470]

#def find_vel(model,vel):
#    #--------------------------
#    #M1p875
#    m1vels = ["0","35","62","83","105","125","146","165","187","207"]
#    #M2
#    m2vels = ["0","36","63","84","106","127","148","168","190","211"]
#    #M2p25
#    m3vels = ["0","36","65","87","109","131","152","173","195","217"]
#    #M2p5
#    m4vels = ["0","37p5","67","89","111","134","156","178","200","222"]
#    #M3
#    m5vels = ["0"]
#    if model[0] == "1p875": vels=m1vels
#    if model[0] == "2p25": vels=m3vels
#    if model[0] == "2": vels=m2vels
#    if model[0] == "2p5": vels=m4vels
#    if model[0] == "3": vels=m5vels
#    #---------------------------
#    
#    return vels[vel]
#    
#def find_sigma(model,vel,par,mode):
#    global temp_freqs
#    v = find_vel(model,vel)
#    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
#    temp_freqs = np.genfromtxt(folder+"temp_freqs")
#            
#    return temp_freqs[mode-1]
#    
#def index_by_name(model,vel,par,mode):
#    global temp_modes
#    v = find_vel(model,vel)
#    folder = homedir+"ROTORCmodels/"+par+"/M"+model[0]+"_V"+v+"/"
#    temp_modes = np.genfromtxt(folder+"temp_MODES")
#    try:
#        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype='|S8',delimiter=8)
#    except:
#        temp_modes = np.genfromtxt(folder+"temp_MODES",dtype=str)
#        temp_modes = [temp_modes[i,0]+" "+temp_modes[i,1] for i in range(len(temp_modes[:,0]))]
#    
#    idx = 999
#    for i in range(len(temp_modes)):
#        tmpmodes = temp_modes[i].split()
#        tmpmode = mode.split()
#        if len(tmpmodes)>1:
#            if tmpmodes[0]==tmpmode[0] and tmpmodes[1]==tmpmode[1]:
#                idx = i+1
#    
#    if idx==999:
#        print "MODE",mode,"NOT FOUND BY NAME!"
#        sys.exit(0)
#    
#    return idx

calc_base_output = True                 
fixed_observer = False
observer_i = 0

mass = "1p875"
vel = 6


modes = ["3 p1"]
freqs = [1.63480]

static_i = [0,10,20,30,40,50,60,70,80,90]
#static_i = [90]

colors = ["red","blue","green","magenta","purple"]



incl = [0,10,20,30,40,50,60,70,80,90]
#incl = [90]


mode_by_name = False

if type(modes[0]) == str:
    mode_by_name = True

parity = []
if mode_by_name==True:
    c_modes = list(modes)
    old_modes = list(modes)
    f_modes = []
    for i in range(len(modes)):
        if int(modes[i].split()[0])%2==0:
            par = "EVEN"
        else:
            par = "ODD"
        
        parity.append(par)
        modes[i] = index_by_name([mass],vel,par,modes[i])
        f_modes.append(find_sigma([mass],vel,par,modes[i]))


vel = find_vel([mass],vel)
homedir = "/home/castaned/Documents/"
loc = homedir+"ROTORCmodels/visibilities/"+mass+"Msun/V"+vel
locbase = homedir+"ROTORCmodels/visibilities/"+mass+"Msun/"
#fmag = glob.glob(loc+"/MODE_EVEN_"+str(51)+"/magnitudes*")
clic_folder = homedir+"pyCLIC/test/"

base = []
b_v_base = []
if calc_base_output==True:
    os.chdir(clic_folder)
    import main as pyclic

    for i in range(len(static_i)):
        os.chdir(clic_folder)
        if (os.path.isfile(locbase+"outputflux_"+mass+"Msun_V"+vel+"_CLIC_s_i"+str(static_i[i])))==False:
            pyclic.run_CLIC(mass+"Msun_V"+vel+"_CLIC_s",[static_i[i]],False,3000.,12099.,2.0,"EVEN",dt_grid,fine_model=False)
            subprocess.call(['cp','outputflux_i'+str(static_i[i])+'.final',locbase+"outputflux_"+mass+"Msun_V"+vel+"_CLIC_s_i"+str(static_i[i])])
        
        os.chdir(locbase)
        if (os.path.isfile("walraven_mags_"+"outputflux_"+mass+"Msun_V"+vel+"_CLIC_s_i"+str(static_i[i])))==False:
            cmag.calc_walraven("outputflux_"+mass+"Msun_V"+vel+"_CLIC_s_i"+str(static_i[i]),[static_i[i]],-1,0)
        
        bfile = np.genfromtxt(locbase+"walraven_mags_"+"outputflux_"+mass+"Msun_V"+vel+"_CLIC_s_i"+str(incl[i]))
        base.append(bfile[4])
        b_v_base.append(-2.5*np.log10(bfile[4]/bfile[5]))
        
    base = np.array(base)
    b_v_base = np.array(base)
    
#os.chdir(vis_path)



"""
if ell==0:
    modes = [51, 69, 78, 85, 92, 100, 106, 114, 122, 129] #l=0 M2p5 V0
    freqs = [1.58717,2.05922,2.49359,2.95717,3.46299,3.99529,4.54267,5.09092,5.64618] #l=0 M2p5 V0
    
if ell==1:
    modes = [5, 61, 81, 89, 96, 104, 111, 133] #l=1 M2p5 V0
    freqs = [0.78296,1.63480,2.16027,2.66899,3.18255,3.70941] #l=1 M2p5 V0
    
if ell==2:
    modes = [8, 34, 52, 66, 76, 84, 90, 97, 105, 113, 128] #l=2 M2p5 V0
    freqs = [0.83464,1.17430,1.60199,1.94807,2.36756,2.86427,3.38705,3.92415,4.47259,5.01979]
    
if ell==3:
    modes = [34, 52, 64, 78, 87, 94, 102, 115, 139]
    freqs = [1.08830,1.42856,1.68216,2.06952,2.54134,3.04737,3.57597] #l=3 M2p5 V=0

if ell==4:
    modes = [12,26,40,49,61,72,80,88]
    freqs = [0.85613,1.03514,1.28624,1.54099,1.79625,2.12890,2.64669,3.17963] # l=4 M2p5 V0
    
if ell==6:
    modes = [67,75,83,91]
    freqs = [1.97688,2.25056,2.82636,3.39865]
"""

mag = []
m_mag = []
b_v = []
def load_modes(sub_dir=""):
    global modes, mag, m_mag, parity
    m_mag_return = []
    m_mag = []
    for i in range(len(modes)):
        for j in incl:
            par = parity[i]
            fmag = glob.glob(loc+"/MODE_"+par+"_"+str(modes[i])+sub_dir+"/walraven_magnitudes"+"_i"+str(j)+"_MODE_"+str(modes[i]))
            #print fmag
            mag.append(np.genfromtxt(fmag[0]))
        m_mag.append(mag)
        mag = []
        m_mag_return = list(m_mag)
    return m_mag_return

def plot_amplitude_ratios():
    global ell,modes,colors,bfile,base,par
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        c = colors[ell/2]
        lstyle = "-"
        if par == "ODD":
            c = colors[-ell]
            lstyle = "--"
            
        #b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
        #plot5 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
        #print mag[i][2]
        #plot4 = plt.vlines(freqs[i], 0, (b_v-b_v_base)/b_v_base,color=c,linestyle=lstyle)
        #mag_ratios = mag[i][3:]/mag[i][3]
        mag_ratios = mag[i][3:]/bfile[3:]
        mag_ratios /= mag_ratios[0]
        plt.plot(wavebands,mag_ratios)
        plt.text(wavebands[-1]+140,mag_ratios[-1],r"$\ell$ = "+str(ell) + " - "+ str(i))
        plt.xlabel("Wavelength")
        plt.ylabel("Apmplitude ratio")


#plot_amplitude_ratios()

def lum_plot():
    global ell,modes,colors,bfile,base,par,mag,freqs,incl,static_i,fixed_observer,observer_i,mass,vel,f_modes,parity
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        #c = colors[ell/2]
        #lstyle = "-"
        #if par == parity[i]:
        #    c = colors[-ell]
        #    lstyle = "--"
            
        #b_v = -2.5*np.log10(mag[i][4]/mag[i][5])
        #plot5 = plt.vlines(freqs[i], 0, np.abs(mag[i][2]-base)/base,color=c,linestyle=lstyle)
        #print mag[i][2]
        plt.vlines(f_modes[i], 0, (m_mag[i][0][2]-base[i])/base[i],color="blue")
        
#lum_plot()
#plt.legend([plot1,plot2,plot3],["$\ell = 0$","$\ell = 2$","$\ell = 4$"],loc="best")
#plt.ylim(0,1.5)
#plt.xlim(1,4)
#plt.title("Mode luminosities M=2.5 V=0")
#plt.ylabel("Normalized L")
#plt.xlabel("Frequency")
#plt.show()
        
def lum_incl_plot():
    global ell,modes,colors,bfile,base,par,mag,freqs,incl,static_i,fixed_observer,observer_i,mass,vel
    for i in range(len(modes)):
        #plt.plot([freqs[i],0],[freqs[i],mag[i][0,2]/base])
        #plt.plot([[freqs[i],0],[freqs[i],1]])
        tmp = []
        for j in range(len(incl)):

            if fixed_observer==True:
                tmp.append((m_mag[i][j][4]-base[static_i[observer_i]])/base[static_i[observer_i]])
            else:
                tmp.append((m_mag[i][j][4]-base[j])/base[j])
        
        plt.plot(incl, (np.array(tmp)))
        plt.text(incl[0]+2, np.abs(np.array(tmp))[0],"$\ell$="+str(c_modes[i].split()[0]))
        
            
    plt.hlines(0,0,90)
    plt.xlabel("Inclination (deg)")
    plt.ylabel("$\Delta L/L$")
    if(fixed_observer==True):
        plt.title("M="+mass+", V="+vel+" - $L$ fixed at $i$="+str(static_i[observer_i]))
    else:
        plt.title("M="+mass+", V="+vel+" - $L$ observer")
    plt.grid()
    
def mag_amplitude_plot(minl,maxl):
    global ell,modes,colors,bfile,base,par,mag,freqs,incl,static_i,fixed_observer,observer_i,mass,vel,c_modes
    tmp =[]
    for j in range(len(incl)):
        #tmp.append(1+(np.abs(minl[0][j][4]-maxl[0][j][4])-base[j])/base[j])
        #tmp.append((np.abs(minl[0][j][4]-maxl[0][j][4])))
        tmp.append(np.abs(-2.5*np.log10(minl[0][j][4]/maxl[0][j][4])))
    
    lbl = "$\ell$="+str(c_modes[0].split()[0])
    lsty = "-"
    if int(c_modes[0].split()[0])>5:
        lsty="--"
    plt.plot(incl,np.array(tmp),label=lbl,ls=lsty)
    #plt.text(incl[0]+2, np.array(tmp)[0],lbl)
    plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
    #plt.hlines(0,0,90)
    plt.xlabel("Inclination [deg]")
    plt.ylabel("|m$_{u-\mathrm{max}}$-m$_{u-\mathrm{min}}$| [mag]")
    plt.title("M="+mass+", V="+vel+", n = 1")
    plt.grid()
    

maxl = load_modes()
minl = load_modes("/minL")

mag_amplitude_plot(minl,maxl)
#lum_incl_plot()
#lum_plot()
#plt.grid()