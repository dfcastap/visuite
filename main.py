# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 01:04:06 2015

@author: diego
"""

import numpy as np
import glob,subprocess,os
import del_dot_xi as ddxi
import pyNRO_1_mode as pyNRO
old_path = os.getcwd()


#MODE info:
par = "temp"
model = ["2p5"]
vel = 1 #index, not velocity!

#tfreq = 1.59692
tfreq = 0.93025

# Info for del dot xi calculation:------------------
# NRO Mode filename:
modefname="MODE1"
norm_f = True
scale = 0.1
depth = 10 #radial zones down from the surface
#---------------------------------------------------




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
#check if visibility file from pulset exits for the selected model:

v_file = glob.glob(rotorc_f+"visibility_file")
if len(v_file)==0:
    os.chdir(rotorc_f)
    r_mod_name = glob.glob("*_ZAMS")
    print subprocess.call(["cp",r_mod_name[0],"orcmod_pulset"])
    subprocess.call([bob_bin+"pulsetnonadb.exe"])
    #print(glob.glob("*_ZAMS"))
    os.chdir(old_path)
v_file = glob.glob(rotorc_f+"visibility_file")
print "v_file OK!"

###### FULL MODE FILE GENERATION:
print pyNRO.run_nro(tfreq,folder)


###### del_DOT_xi>

#xi_r,xi_t,dt_t,r = ddxi.calcdeldotxi(par,model,vel,modefname)
        
#xi_r_n,xi_t_n,dt_t_n = ddxi.norm_and_scale(xi_r,xi_t,dt_t,norm_f,scale,depth)

#xi_r_rot,xi_t_rot,dt_t_rot = ddxi.to_rotorc(xi_r_n,xi_t_n,dt_t_n)

