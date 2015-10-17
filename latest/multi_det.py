#!/usr/bin/python

#########################################################
#                                                       #
#              log likelihood calculator                #
#               Synth. Stream approach                  #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
Generates timing residuals and save them as a .res file.
'''
#libraries
import numpy as np
import F_c_F_s
import plotter
import MX_producer
import random
from multiprocessing import Process, Queue
import argparse

#number of grid points
#x=100
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--loc', help='pulsar location file', default='pulsar.loc')
args = parser.parse_args()

#read the pulsar locations
f= open(args.loc,'r')

#pulsar coordinates lists
phi=[]
theta=[]



for line in f:
                        #Read and Split
                        line = line.strip()
                        columns = line.split()
                        #Read pulsar coordinates and put them in theta and phi lists 
                        theta.append(float(columns[0]))
                        phi.append(float(columns[1]))

#print 'Pulsar locations read'

#no of pulsars
n=len(theta)


# Theta phi-grid

t=np.load('theta.grid.npy')
p=np.load('phi.grid.npy')

x=len(t)
#t=np.linspace(0, np.pi, num=x)
#p=np.linspace(0,2*np.pi,num=2*x)

#read the pulsar locations
g= open('source.loc','r')

#pulsar coordinates lists
s_phi=[]
s_theta=[]



for line in g:
                        #Read and Split
                        line = line.strip()
                        columns = line.split()
                        #Read pulsar coordinates and put them in theta and phi lists 
                        s_theta.append(float(columns[0]))
                        s_phi.append(float(columns[1]))

print 'Pulsar locations read'


s_theta=np.array(s_theta)
s_phi=np.array(s_phi)

#include the actual source locations too!

t=np.sort(np.append(t,s_theta))
p=np.sort(np.append(p,s_phi))

ind_t=t.tolist().index(s_theta)
ind_p=p.tolist().index(s_phi)

#make a mesh now
[T,P]=np.meshgrid(t,p)

print 'Calculating Log Likelihood pattern'

#the heart of the program
print "Computing Z"

def z_grid(n_src,ret_val_1):
    global T,P,theta,phi,s_theta,s_phi
    [z,zL,zR,pL,pR]=MX_producer.Z_maker(T,P,theta,phi,n_src)
    i,j = np.unravel_index(z.argmax(), z.shape)
    plotter.llr_contour(z,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD_C")
    plotter.llr_contour(zL,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD_L")
    plotter.llr_contour(zR,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD_R")
    z_L=zL[i,j]
    z_R=zR[i,j]
    p_L=pL[i,j]
    p_R=pR[i,j]
    Mod_F_c=(F_c_F_s.Mod_F_DP_c(theta,phi,T,P)[i,j])
    Mod_F_s=(F_c_F_s.Mod_F_DP_s(theta,phi,T,P)[i,j])
    Y=(np.sqrt(z_L)*Mod_F_s*complex(np.cos(p_L),np.sin(p_L)))/(np.sqrt(z_R)*Mod_F_c*complex(np.cos(p_R),np.sin(p_R)))
    quant1=np.absolute( ((1j*Y)-1)/ (1j*Y+1) )
    print "DP_theta     : ",F_c_F_s.DP_theta(theta,phi,T,P)[i,j]
    print "Mod_F_DP_+   : ",Mod_F_c
    print "Mod_F_DP_x   : ",Mod_F_s
    print "Max LLR      : ",z_L+z_R
    print "rho_L^2      : ",z_L
    print "rho_R^2      : ",z_R
    print "phase_L      : ",p_L
    print "phase_R      : ",p_R
    print "Recov. Loc   : ",T[i,j],P[i,j]
    print "Cos(4 Psi)   : ",np.cos(np.angle(((1j*Y)-1)/ (1j*Y+1)) + F_c_F_s.DP_theta(theta,phi,T,P)[i,j])
    print "Cos(inc)     : ",((1-np.sqrt(quant1))/(1+np.sqrt(quant1)))
    ret_val1.put(z)

def area(area_val,ret_val_2):
    global theta,phi
    theta_area=[]
    phi_area=[]
    for i in range(len(t)):
    	theta_area.append(random.uniform(s_theta-area_val, s_theta+area_val))
    	phi_area.append(random.uniform(s_phi-area_val, s_phi+area_val))
    [T_area,P_area]=np.meshgrid(theta_area,phi_area)
    [z_area,a,b,c,d]=MX_producer.Z_maker(T_area,P_area,theta,phi,1)
    i,j = np.unravel_index(z_area.argmax(), z_area.shape)
    area=float(np.sum((z_area > 0.85*(z_area[i,j]))))/np.prod(z_area.shape)*900/np.sqrt(z_area[i,j])
    ret_val2.put(area)

if __name__=='__main__':
    ret_val1=Queue()
    ret_val2=Queue()
    n_src=1
    p1 = Process(target = z_grid, args=(n_src, ret_val1))
    p1.start()
    area_val=(0.17453292519*3)
    p2 = Process(target = area, args=(area_val,ret_val2))
    p2.start()
    z=ret_val1.get()
    area=ret_val2.get()
    print "Area > 85%   : ",area
#    i,j = np.unravel_index(z.argmax(), z.shape)
    #Post calculations, for recovery of physical parameters
#    z_L=zL[i,j]
#    z_R=zR[i,j]
#    p_L=pL[i,j]
#    p_R=pR[i,j]

#    #print F_c_F_s.F_dot(theta,phi,T,P)
#    Mod_F_c=(F_c_F_s.Mod_F_DP_c(theta,phi,T,P)[i,j])
#    Mod_F_s=(F_c_F_s.Mod_F_DP_s(theta,phi,T,P)[i,j])
#
#    print "DP_theta     : ", F_c_F_s.DP_theta(theta,phi,T,P)[i,j]
#    print "Mod_F_DP_+   : ", (F_c_F_s.Mod_F_DP_c(theta,phi,T,P)[i,j])
#    print "Mod_F_DP_x   : ", (F_c_F_s.Mod_F_DP_s(theta,phi,T,P)[i,j])
#    print "Max LLR      : ",zL[i,j]+zR[i,j]
#    print "rho_L^2      : ",z_L
#    print "rho_R^2      : ",z_R
#    print "phase_L      : ",p_L
#    print "phase_R      : ",p_R
#    print "Recov. Loc   : ",T[i,j],P[i,j]
#
#    #recovering the initial parameters
#    Y=(np.sqrt(z_L)*Mod_F_s*complex(np.cos(p_L),np.sin(p_L)))/(np.sqrt(z_R)*Mod_F_c*complex(np.cos(p_R),np.sin(p_R)))
#    print "Cos(4 Psi)   : ",np.cos(np.angle(((1j*Y)-1)/ (1j*Y+1)) + F_c_F_s.DP_theta(theta,phi,T,P)[i,j])
#
#    quant1=np.absolute( ((1j*Y)-1)/ (1j*Y+1) )
#
#    print "Cos(inc)     : ",((1-np.sqrt(quant1))/(1+np.sqrt(quant1)))
#
#    #plots
#    plotter.llr_contour(zL,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD_L")
#    plotter.llr_contour(zR,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD_R")
#    plotter.llr_contour(z,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"MD")
#    print "Log Likelihood map saved as \'MD.png, MD_L.png, MD_R.png\'"
