import numpy as np
#import pylab as pl
import subprocess,pexpect,glob,os

oldpath = os.getcwd()    
    
#folder = "/home/diego/Documents/ROTORCmodels/"+par+"/M"+model[0]+"_V"+vels[vel]+"/"

#print folder
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


def run_nro(tfreq,folder,model,vel,par,nmode,homedir):
    global freq,freqs,rad_nodes,smixes,file_count
    global c_flag,index,tempfreq,times,dmod

    nro_loc = homedir+'NROe/del_dot_xi_mod/nro'
    c_flag = 99
    os.chdir(folder)
    dmod = glob.glob("Dmod_"+model+"M_V"+vel)[0]
    print dmod
    steps = [1.1e-4,4.56e-4,8e-4,1e-4,1.3e-4]
    step = steps[0]
    
    #if ifreq==0:
    #    ifreq = np.round(float(subprocess.check_output(['tail','-1','NRO_pts']).strip("\n").split("  ")[0]) - 1.5*step,5)
    #    tempfreq = 1.*ifreq
    try:    
        subprocess.call(['mv','cache','cache.bak'])
    except:
        print "No cache!"
    child = pexpect.spawn (nro_loc)
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
    child.sendline ('1 '+str(tfreq-2.*step))
    child.expect ('Any change [n value]?')
    step = steps[0]
    child.sendline ('3 '+str(step))
    child.expect ('Any change [n value]?')
    child.sendline ('4 201')
    child.expect ('Any change [n value]?')
    child.sendline ('15')
    while c_flag!=-1:
        index = child.expect (['CONVERGENCE!', '...end of range!','PROBABLE BAD MODE!',pexpect.EOF, pexpect.TIMEOUT])
        if index==0:
            print child.after
            child.expect("Frequency: ")
            freq = np.float(child.readline())
            print "f = ",freq
            #freqs.append(freq)
            
            if np.round(freq,4) == np.round(tfreq,4): 
                c_flag=1
            else:
                c_flag=0
            
            c_flag=1
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
            c_flag = 3
            print "end of range..."
            
        if index==2:
            print "Bad mode! Trying again..."
            c_flag = 0
        if index==4:
            #np.savetxt("temp_freqs",freqs,fmt='%s')
            #np.savetxt("temp_radnodes",rad_nodes,fmt='%s')
            #np.savetxt("temp_smixes",smixes,fmt='%s')
            print "timeout!"
            c_flag = -1
            child.terminate()
            return 1
        
        if c_flag==1:
            child.expect('Option?')
            child.sendline('3')
            
            child.expect('Continue scan (y/n)?')
            child.sendline('n')
            
            child.expect('Option?')
            child.sendline('4')
            c_flag=-1
            subprocess.call(['mv','MODE1','MODE_'+par+'_'+str(nmode)])
    
        if c_flag==0:
            child.expect ('Option?')
            child.sendline ('1')
            
            child.expect('Continue scan (y/n)?')
            child.sendline('n')
            
            child.expect('Option?')
            child.sendline('3')
            
            child.expect ('Any change [n value]?')
            child.sendline ('1 '+str(tfreq-step))
            
            child.expect ('Any change [n value]?')
            child.sendline ('15')
            #tempfreq = np.float(tfreq) 
            #child.expect ('Any change [n value]?')
            #child.sendline ('1 '+str(tfreq-step))
            #child.expect ('Any change [n value]?')
            #child.sendline ('15')
            
        if c_flag==3:
            
            child.expect('Option?')
            child.sendline('3')
            
            child.expect ('Any change [n value]?')
            child.sendline ('15')
    
    return 0
    #child.sendeof()
    #print child.readlines()


#print run_nro(tfreq,folder)
os.chdir(oldpath)
