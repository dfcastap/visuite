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
from main_modules import find_vel, find_sigma, index_by_name, list_gen, emode, load_mags

dt_grid = 250
vis_path = os.getcwd()

global homedir
if (os.path.isfile(vis_path+"/lachesis"))==False:
    homedir = "/home/castaned/Documents/"
else:
    homedir = "/home/castaned/"
    
sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 12,'axes.titlesize':16,'font.size':14})
wavebands = [3650.,4450.,5510.,6580.,8060.]

wal_wave = [3250,3630,3840,4320,5470]



def load_modes(sub_dir=""):
    global modes, mag, m_mag, parity
    m_mag_return = []
    m_mag = []
    for i in range(len(modes)):
        for j in incl:
            par = parity[i]
            fmag = glob.glob(loc+"/MODE_"+par+"_"+str(modes[i])+sub_dir+"/walraven_magnitudes"+"_i"+str(j)+"_MODE_"+str(modes[i]))
            #print loc+"/MODE_"+par+"_"+str(modes[i])+sub_dir+"/walraven_magnitudes"+"_i"+str(j)+"_MODE_"+str(modes[i])
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
        plt.text(incl[0]+2, np.abs(np.array(tmp))[0],r"$\ell$="+str(c_modes[i].split()[0]))
        
            
    plt.hlines(0,0,90)
    plt.xlabel(r"Inclination (deg)")
    plt.ylabel(r"$\Delta L/L$")
    if(fixed_observer==True):
        plt.title(r"M="+mass+", V="+vel+" - $L$ fixed at $i$="+str(static_i[observer_i]))
    else:
        plt.title(r"M="+mass+", V="+vel+" - $L$ observer")
    plt.grid()
    
def mag_amplitude_plot(mags,col_idx,lbl,ref=1.,inlog=False):
    global ell,modes,colors,bfile,base,par,mag,freqs,incl,static_i,fixed_observer,observer_i,mass,vel,c_modes

    tmp = np.zeros((len(incl),len(mags)/2))
    for j in range(len(incl)):
        #tmp.append(1+(np.abs(minl[0][j][4]-maxl[0][j][4])-base[j])/base[j])
        #tmp.append((np.abs(minl[0][j][4]-maxl[0][j][4])))
        tmp[j,0]=(np.abs(-2.5*np.log10(mags[0][j][col_idx]/mags[1][j][col_idx])))
        #tmp.append(np.abs(-2.5*np.log10(minl[0][j][4]/bnum))+np.abs(-2.5*np.log10(maxl[0][j][4]/bnum)))
    
    #lbl = r"$\ell$="+str(c_modes[0].split()[0])
    lsty = "-"
    if int(c_modes[0].split()[0])>7:
        lsty="-"
    
    clr = sns.color_palette("Set2", 12)[int(c_modes[0].split()[0])]
    
    if not inlog:
        #plt.plot(incl,np.array(tmp)/ref,label=lbl,color=clr,ls=lsty)
        plt.plot(incl,np.array(tmp)/ref,label=lbl,ls=lsty)
        plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
    else:
        plt.semilogy(incl,np.array(tmp)/ref,label=lbl,color=clr,ls=lsty)
    
    #plt.text(incl[0]+2, np.array(tmp)[0],lbl)
    #plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))

    #plt.xlabel("Inclination [deg]")
    #plt.ylabel("|m$_{u-\mathrm{max}}$-m$_{u-\mathrm{min}}$| [mag]")
    #plt.title("M="+mass+", V="+emode(mass,vel,c_modes[0]).vel+", n = 1")


def mag_amplitude_vs_color_plot(mags,col_idx,yshift,ref=1.,inlog=False):
    global ell,modes,colors,bfile,base,par,mag,freqs,incl,static_i,fixed_observer,observer_i,mass,vel,c_modes
    jh_wave = [3650.,4450.,5510.,6580.,8060.]

    wal_wave = [3250,3630,3840,4320,5470]
    
    m_cols = [3,4,5,6,7]
    
    tmp = np.zeros((len(incl),len(m_cols)))
    for j in range(len(incl)):
        #tmp.append(1+(np.abs(minl[0][j][4]-maxl[0][j][4])-base[j])/base[j])
        #tmp.append((np.abs(minl[0][j][4]-maxl[0][j][4])))
        for i in range(len(m_cols)):
            tmp[j,i]=(np.abs(-2.5*np.log10(mags[0][j][m_cols[i]]/mags[1][j][m_cols[i]])))
        #tmp.append(np.abs(-2.5*np.log10(minl[0][j][4]/bnum))+np.abs(-2.5*np.log10(maxl[0][j][4]/bnum)))
    
    lbl = r"$\ell$="+str(c_modes[0].split()[0])
    lsty = "-"
    if int(c_modes[0].split()[0])>7:
        lsty="--"
    
    mrkrs = [r"$"+str(i)+"$" for i in range(len(incl))]
    #clr = sns.color_palette("coolwarm", 10)
    clr = sns.color_palette("hls", 10)
    
    for i in range(len(incl)):
        if not inlog:
            plt.plot(wal_wave,tmp[i,:]+yshift*i,marker=mrkrs[i],ms=10,color=clr[i],label="i="+str(i*10)+"$^{\circ}$")
            plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
        else:
            plt.semilogy(incl,np.array(tmp)/ref,label=lbl,color=clr,ls=lsty)
    
    
    #plt.text(incl[0]+2, np.array(tmp)[0],lbl)
    #plt.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
    plt.xticks(wal_wave, ["w","u","l","b","v"])
    #plt.xlabel("Inclination [deg]")
    #plt.ylabel("|m$_{u-\mathrm{max}}$-m$_{u-\mathrm{min}}$| [mag]")
    #plt.title("M="+mass+", V="+emode(mass,vel,c_modes[0]).vel+", n = 1")


    
def run_plots(mass, vel, mode):
    global modes, parity, old_modes, f_modes,colors,observer_i,freqs,fixed_observer, loc, locbase,base,b_v_base, incl,c_modes
    calc_base_output = True                 
    fixed_observer = False
    observer_i = 0
    

    
    modes = [mode]
    
    static_i = [0,10,20,30,40,50,60,70,80,90]

    
    colors = ["red","blue","green","magenta","purple"]
      
    
    #incl = [0,10,20,30,40,50,60,70,80,90]
    #incl = list(np.arange(0,95,5))
    
    
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
    clic_folder = homedir+"pyCLIC/test/"
    
    base = []
    b_v_base = []
    if calc_base_output==False:
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
    
    
ii = 8
with sns.color_palette("rainbow",n_colors=9):
    for okc in [6]:
        print str(okc)
        collect = []
        lbls = []
        for ocity in [0,1,2,3,4,5,6,7,8]:
            pfol = "/home/castaned/Documents/ROTORCmodels/visibilities/plots/toBob_apr22/mag_ampl/"
            mag = []
            m_mag = []
            b_v = []
            mass = "1p875"
            vel = ocity
            r_ord = 4
            kind = "g"
            ells = [okc]
            mds = list_gen(ells,r_ord,r_ord,kind)
            p_phase = False
            plog = False
            col_idx = 4
            #f, axarr = plt.subplots(2, sharex=True)
            yshift = 0.0025
            incl = list(np.arange(0,95,10))
            
            for  mm in mds:
                
                em = emode(mass,vel,mm)
                run_plots(em.mass, em.vel_idx, em.name)
                
            #    maxl = load_modes()
            #    minl = load_modes("/minL")
            #    maxl = load_modes("/phi_max")
            #    minl = load_modes("/phi_min")
            
                wal_max,jh_max = load_mags("",em.mass,em.vel_idx,em.parity,em.index,incl)
                wal_min,jh_min = load_mags("/minL",em.mass,em.vel_idx,em.parity,em.index,incl)
                
#                wal_phi_max,jh_phi_max = load_mags("/phi_max",em.mass,em.vel_idx,em.parity,em.index,incl)
#                wal_phi_min,jh_phi_min = load_mags("/phi_min",em.mass,em.vel_idx,em.parity,em.index,incl)
#                
#                if mds.index(mm)==0:
#                    wal_ref = np.abs(-2.5*np.log10(wal_max[0][col_idx]/wal_min[0][col_idx]))
#                    wal_phi_ref = np.abs(-2.5*np.log10(wal_phi_max[0][col_idx]/wal_phi_min[0][col_idx]))
                
                wal_ref = 1.
                wal_phi_ref = 1.
                
                #lbl = r"$\ell$="+str(em.name.split()[0])
                lbl = r"V="+str(em.eq_vel)+" km s$^{-1}$"
                #lbl = r"$\ell$="+str(em.name.split()[0]) #+" - cut"
                
                y = (float(em.eq_vel.replace("p","."))/100.)**2
                
                if not p_phase:
                    mag_amplitude_plot([wal_min,wal_max],col_idx,lbl,ref=wal_ref,inlog=plog)
                    #mag_amplitude_vs_color_plot([wal_min,wal_max],col_idx,yshift)
                    collect.append([np.abs(-2.5*np.log10(wal_max[ii][col_idx]/wal_min[ii][col_idx])),y])
                    
                    lblphase= "0"
                else:
                    mag_amplitude_plot([wal_phi_min,wal_phi_max],col_idx,lbl,ref=wal_phi_ref,inlog=plog)
                    #mag_amplitude_vs_color_plot([wal_phi_min,wal_phi_max],col_idx,yshift)
                    collect.append([np.abs(-2.5*np.log10(wal_phi_max[ii][col_idx]/wal_phi_min[ii][col_idx])),y])
                    lblphase = "123"
     
#                barr = np.array(wal_max)
#                marr = np.array(wal_min)
#                aarr = np.abs(-2.5*np.log10(barr[:,col_idx]/marr[:,col_idx]))
#                xx = np.argmax(aarr)
#                collect.append([aarr[xx],y])
                
                lbls.append(lbl)
                #handles, labels = ax.get_legend_handles_labels()
                #hold = handles[1]
                #lold = handles[1]
                #mnamelst = [r"$\ell=$"+n.split()[0] for n in mmod.list_gen([0,1,2,3,4,5,6,7,8,9,10],3,3,"p") ]
                #plt.legend(mnamelst,loc="best")
                #plt.ylim(0,2.5e-1)
                #lum_incl_plot()
                #lum_plot()
        
#        collect = np.array(collect)
#        est = "-"
#        if int(em.name.split()[0])>5:
#            est = "--"
#        #with sns.color_palette("rainbow",n_colors=5):
#        plt.plot(collect[:,1],collect[:,0],ls=est,label=lbl)
###       plt.ylim(-0.05,0.45)
    
plt.xlabel("Inclination [deg]")
#plt.ylabel("Scaled Mag. Ampl.")
plt.ylabel("Mag. Amplitude [mag]")
#plt.ylabel("Max Mag. Amplitude [mag]")
#plt.xlabel(r"$(V_{\mathrm{eq}} / V_{\mathrm{ref}})^2$")
#plt.ylabel("Mag. Amplitude + i*"+str(yshift)+" [mag]")
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$, V="+emode(mass,vel,mds[0]).eq_vel+"km s$^{-1}$, "+kind+" = "+str(r_ord)+" - L_phase = "+lblphase+"$^{\circ}$")
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$ - $\ell=$"+emode(mass,vel,mds[0]).name+" - L_phase = "+lblphase+"$^{\circ}$")
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$ -  L_phase = "+lblphase+"$^{\circ}$"+" - i = "+str(incl[ii])+"$^{\circ}$")
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$ - i = "+str(incl[ii])+"$^{\circ}$ - n = "+str(r_ord))
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$ - V = "+emode(mass,vel,mds[0]).eq_vel+"km s$^{-1}$, n = "+str(r_ord))
plt.title(r"M="+mass.replace("p",".")+"M$\odot$ - $\ell$ = "+mds[0].split()[0]+", n = -"+str(r_ord))
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$, V="+emode(mass,vel,mds[0]).eq_vel+"km s$^{-1}$, "+mm+ " - L_phase = "+lblphase+"$^{\circ}$")
#plt.title(r"M="+mass.replace("p",".")+"M$\odot$ , n = "+str(r_ord))#+" - i = "+str(incl[ii])+"$^{\circ}$")
plt.grid()

plt.legend(loc="best")
#plt.ylim(0,1.2)
#plt.savefig(pfol+"/M"+mass+"/V"+emode(mass,vel,mds[0]).eq_vel+"_ell-"+str(ells[0])+"-"+str(ells[-1])+"_"+kind+str(r_ord)+"_Lphase_"+lblphase+".png")
#plt.savefig(pfol+"/M"+mass+"_mag_ampl_vs_color/"+str(ells[0])+"_"+kind+str(r_ord)+"/V"+emode(mass,vel,mds[0]).eq_vel+"_ell-"+str(ells[0])+"_"+kind+str(r_ord)+"_Lphase_"+lblphase+".png")

#mnamelst = [r"$\ell=$"+n.split()[0] for n in list_gen([0,1,2,3,4,5,6,7,8,9,10],3,3,"p") ]
#plt.legend(mnamelst,loc="best")