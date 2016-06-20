# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 18:42:32 2016

@author: castaned
"""

import numpy as np
import seaborn as sns
import main_modules as mmod
from matplotlib import cm
import pylab as plt

sns.set(style="white",rc={"figure.figsize": (8, 8),'axes.labelsize': 16,
                              'ytick.labelsize': 12,'xtick.labelsize': 12,
                              'legend.fontsize': 10,'axes.titlesize':18,'font.size':14})
sns.set_style({'legend.frameon': True})
                              

mass = "2p25"
parity = "ODD"
select_modes = []
tag = ["1 p5"]
sym = ["o","s","d","^","*","<",">","h"]
lbl_track = []
xlm = [3.2,4.3]
vels = range(10)
ax = plt.subplot(111)

global tags
tags = np.zeros((len(vels),len(tag)+1))

def plot_point_evol(mass,parity,select_modes,sym,lbl_track,xlm,ax):
    for i in vels:
        print i
        for j in range(90,110):
            tmpmode = mmod.emode(mass,i,j,par=parity)
            c = int(tmpmode.name.split()[0])
            #clr = cm.viridis((15-c)*255/15)
            #clr = sns.color_palette("gnuplot", 18)[c]

#            if tmpmode.frequency>xlm[1]:
#                #print tmpmode.name, tmpmode.index
#                break                   
            
            if tmpmode.frequency>xlm[0]:
                if c!=99:
                    clr = sns.diverging_palette(255, 133, l=60, n=16, center="dark")[c]
                    mkr = sym[c/2]
                if c==99:
                    clr = "k"
                    mkr = "o"
                y = (float(tmpmode.eq_vel.replace("p","."))/100.)**2
                
                for t in tag:
                    if t.split()[0]=="0":
                        if (t.find("h")>0 or t.find("H")>0) and (tmpmode.name.find("h")>0 or tmpmode.name.find("H")>0) and (t.split()[0] == tmpmode.name.split()[0]):
                            print t,tmpmode.name,y
                            tags[i,tag.index(t)+1] = tmpmode.frequency
                            tags[i,0] = y
                        elif (t.find("f")>0 or t.find("F")>0) and (tmpmode.name.find("f")>0 or tmpmode.name.find("F")>0) and (t.split()[0] == tmpmode.name.split()[0]):
                            print t,tmpmode.name,y
                            tags[i,tag.index(t)+1] = tmpmode.frequency
                            tags[i,0] = y
                    elif (t.find("f")>0 or t.find("F")>0) and (tmpmode.name.find("f")>0 or tmpmode.name.find("F")>0) and (t.split()[0] == tmpmode.name.split()[0]):
                        print t,tmpmode.name,y
                        tags[i,tag.index(t)+1] = tmpmode.frequency
                        tags[i,0] = y
                    elif t==tmpmode.name:
                        print t,tmpmode.name,y
                        tags[i,tag.index(t)+1] = tmpmode.frequency
                        tags[i,0] = y
                    
                if len(select_modes)==0:
                    plt.plot(tmpmode.frequency,y,ls="",marker=mkr,color=clr)
                    if i==0:
                        if tmpmode.frequency>xlm[0] and tmpmode.frequency<xlm[1]:
                            #plt.text(tmpmode.frequency,-0.25,str(j),fontsize=10,horizontalalignment='center',color=clr)
                            if c not in lbl_track:
                                lbl_track.append(c)
                                plt.plot([],ls="",marker=mkr,color=clr,label=str(c))
                else:
                    if c in select_modes:
                        clr = sns.diverging_palette(255, 133, l=60, n=len(select_modes), center="dark")[select_modes.index(c)]
                        plt.plot(tmpmode.frequency,y,ls="",marker=mkr,color=clr)
                        if i==0:
                            if tmpmode.frequency>xlm[0] and tmpmode.frequency<xlm[1]:
                                #plt.text(tmpmode.frequency,-0.25,str(j),fontsize=10,horizontalalignment='center',color=clr)
                                if c not in lbl_track:
                                    lbl_track.append(c)
                                    plt.plot([],ls="",marker=mkr,color=clr,label=str(c))
    
    global handles, labels, hl, handles2, labels2
    handles, labels = ax.get_legend_handles_labels()
    labels = [int(i) for i in labels]
    import operator
    hl = sorted(zip(handles, labels),
                key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)
    labels2 = [r"$\ell$="+str(i) for i in labels2]
    
    ax.legend(handles2, labels2,loc='lower center', ncol=8)
    plt.xlabel(r"$\omega/(4\pi G \rho_{\mathrm{ref}})^{1/2}$")
    plt.ylabel(r"$(V_{\mathrm{eq}} / V_{\mathrm{ref}})^2$")
    plt.title(r"M="+mass.replace("p",".")+"M$_{\odot}$")
    curr_ylim = plt.ylim()
    plt.ylim(-0.5,curr_ylim[1]+0.5)
    plt.xlim(xlm[0],xlm[1])
    plt.grid()
    #plt.legend(loc='lower center', ncol=8)
    
    if len(tag)>0:
        for t in tag:
            plt.plot(tags[:,tag.index(t)+1],tags[:,0],"--",color="k",alpha=0.25)
    
def plot_line_evol(mass,parity,select_modes,sym,lbl_track,xlm,ax):

    for j in select_modes:
        mode_fs = []
        mode_vels = []
        
        for i in range(10):
            #print j,i
            tmpmode = mmod.emode(mass,i,j,par=parity)
            mode_fs.append(tmpmode.frequency)
            mode_vels.append(float(tmpmode.eq_vel))
        

        c = int(tmpmode.name.split()[0])
        if c!=99:
            clr = sns.diverging_palette(255, 133, l=60, n=8, center="dark")[c/2]

        mode_fs = np.array(mode_fs)
        mode_vels = np.array(mode_vels)
        plt.plot(mode_fs,(mode_vels/100.)**2,color=clr,marker="o",ms=5)
        
        
        plt.xlabel(r"$\omega/(4\pi G \rho_{\mathrm{ref}})^{1/2}$")
        plt.ylabel(r"$(V_{\mathrm{eq}} / V_{\mathrm{ref}})^2$")
        plt.title(r"M="+mass.replace("p",".")+"M$\odot$")
        curr_ylim = plt.ylim()
        plt.ylim(-0.5,curr_ylim[1])
        plt.xlim(xlm[0],xlm[1])
        plt.grid()
        #plt.legend(loc='lower center', ncol=8)
    
plot_point_evol(mass,parity,select_modes,sym,lbl_track,xlm,ax)
#plot_line_evol(mass,parity,mmod.list_gen([1,5],1,2,"p"),sym,lbl_track,xlm,ax)