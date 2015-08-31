import numpy as np
#import pylab as pl
import subprocess,pexpect,glob,os

oldpath = os.getcwd()


ifreq = 1.59691
tfreq = 1.59692
steps = [1e-5,5e-5,8e-5,1e-4,1.3e-4]
nto_loc = '/home/diego/Documents/NROe/del_dot_xi_mod/nro'

par = "temp"
model = ["2p5"]
vel = 0 #index, not velocity!
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
os.chdir(folder)

print folder
"""
if len(glob.glob("temp_freqs"))==0:
    freqs = []
    rad_nodes = []
    smixes = []
else:
    with open('temp_freqs') as f:
        freqs = f.read().splitlines()
    with open('temp_radnodes') as f:
        rad_nodes = f.read().splitlines()
    with open('temp_smixes') as f:
        smixes = f.read().splitlines()
    ifreq = np.round(float(subprocess.check_output(['tail','-1','NRO_pts']).strip("\n").split("  ")[0]))
    if ifreq <= 1.5:
        step = steps[0]
    if ifreq > 1.5 and ifreq <= 3.:
        step = steps[1]
    if ifreq > 3. and ifreq <= 5.:
        step = steps[2]
    if ifreq > 5. and ifreq <= 8.:
        step = steps[3]
    if ifreq > 8.:
        step = steps[4]
    ifreq = ifreq - 1.5*step
"""
dmod = glob.glob("Dmod*")[0]
print dmod
#dmod = 'Dmod_2M_V0'

nmode = 190
c_flag = 99
times = 0
freqs = []
def run_nro(ifreq):
    global freq,freqs,rad_nodes,smixes,step,file_count,nmode
    global c_flag,index,tempfreq,times,dmod
    global nro_loc
    #if ifreq==0:
    #    ifreq = np.round(float(subprocess.check_output(['tail','-1','NRO_pts']).strip("\n").split("  ")[0]) - 1.5*step,5)
    #    tempfreq = 1.*ifreq
    try:    
        subprocess.call(['mv','cache','cache.bak'])
    except:
        print "No cache!"
    child = pexpect.spawn (nto_loc)
    child.timeout= 100
    fout = file('mylog.txt','w')
    child.logfile_read = fout
    child.expect ('filename')
    child.sendline (dmod)
    child.expect ('level tweak factor')
    child.sendline ('7')
    child.expect ('Option?')
    child.sendline ('3')
    child.expect ('Any change [n value]?')
    child.sendline ('1 '+str(ifreq))
    child.expect ('Any change [n value]?')
    step = steps[0]
    child.sendline ('3 '+str(step))
    child.expect ('Any change [n value]?')
    child.sendline ('4 100')
    child.expect ('Any change [n value]?')
    child.sendline ('15')
    while c_flag!=-1:
        index = child.expect (['CONVERGENCE!', '...end of range!','PROBABLE BAD MODE!',pexpect.EOF, pexpect.TIMEOUT])
        if index==0:
            print child.after
            child.expect("Frequency: ")
            freq = np.float(child.readline())
            #freqs.append(freq)
            if freq == tfreq: c_flag=1
            #child.expect("Radial nodes: ")
            #rad_node = np.str(child.readline()).strip('\r\n')
            #rad_nodes.append(rad_node)
            #child.expect("Smix: ")
            #smix = np.str(child.readline()).strip('\r\n')
            #smixes.append(smix)
            #np.savetxt("temp_freqs",freqs,fmt='%s')
            #np.savetxt("temp_radnodes",rad_nodes,fmt='%s')
            #np.savetxt("temp_smixes",smixes,fmt='%s')
            #file_count = len(freqs)
            #subprocess.call(['cp','temp_equator','temp_equator'+str(file_count)])
            #subprocess.call(['cp','temp_surf','temp_surf'+str(file_count)])
        if index==1:  
            c_flag = 1
        elif index==2:
            print "Bad mode!"
            c_flag = 0
        elif index==4:
            #np.savetxt("temp_freqs",freqs,fmt='%s')
            #np.savetxt("temp_radnodes",rad_nodes,fmt='%s')
            #np.savetxt("temp_smixes",smixes,fmt='%s')
            print "timeout!"
            c_flag = -1
            child.terminate()
            return "nmodes= ", len(freqs)
        
        if c_flag==1:
            child.expect('Option?')
            child.sendline('3')
            
            child.expect('Continue scan (y/n)?')
            child.sendline('n')
            
            child.expect('Option?')
            child.sendline('4')
            c_flag=-1
    
        elif c_flag==0:
            child.expect ('Option?')
            child.sendline ('3')
            child.expect('Eigenfrequency:')
            tempfreq = np.float(ifreq) 
            child.expect ('Any change [n value]?')
            child.sendline ('15')
    
    return "done! freq = "+ freq 
    #child.sendeof()
    #print child.readlines()

print run_nro(ifreq)

