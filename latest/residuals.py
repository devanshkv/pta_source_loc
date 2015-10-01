#!/usr/bin/python

#########################################################
#                                                       #
#               Timing Residual Generator               #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
Generates timing residuals and save them as a .res file.
Generates time stamp and save them as .tim file
'''
#libraries
import numpy as np
import F_c_F_s

#read the pulsar locations

f= open('pulsar.loc','r')

#pulsar coordinates lists
phi_p=[]
theta_p=[]



for line in f:
                        #Read and Split
                        line = line.strip()
                        columns = line.split()
                        #Read pulsar coordinates and put them in theta and phi lists.
                        theta_p.append(float(columns[0]))
                        phi_p.append(float(columns[1]))

f.close()
print 'Pulsar locations read'

#read the GW locations

f= open('source.loc','r')
#pulsar coordinates lists
phi_s=[]
theta_s=[]



for line in f:
                        #Read and Split
                        line = line.strip()
                        columns = line.split()
                        #Read pulsar coordinates and put them in theta and phi lists.
                        theta_s.append(float(columns[0]))
                        phi_s.append(float(columns[1]))

f.close()

print 'Source locations read'


theta_p=np.array(theta_p)
phi_p=np.array(phi_p)
theta_s=np.array(theta_s)
phi_s=np.array(phi_s)

#Initial parameters
m1  =float(7.0e08)
m2  =float(2.0e09)
mc  = ((m1*m2)**(3.0/5.0))/((m1+m2)**(1.0/5.0))
f   =float(10e-1)/(2*np.pi)
Dl  =float(1e27)
A   =1#2.0*(mc**(5.0/3.0))*((np.pi*f)**(2.0/3.0))/Dl
inc =float(np.pi/4)
phi0=float(0)
psi =float(0)
m=len(theta_p)

print "Parameters"
print "Chirp Mass   : ",'%e' %mc
print "Amplitude    : ",'%e' %A
print "frequency    : ",'%e' %f
print "Inclination  : ",inc
print "Initial Phase: ",phi0
print "Pol. angle   : ",psi
print "No. of Pulsar: ",m

print 'Computing source parameters'
n_src=1
a=np.zeros(4*n_src)
j=0
for k in range(n_src):
    a[j+0]=-A*(((1+(np.power(np.cos(inc),2)))*np.cos(2*psi)*np.cos(phi0))-(2*np.cos(inc)*np.sin(2*psi)*np.sin(phi0)))
    a[j+1]=-A*(((1+(np.power(np.cos(inc),2)))*np.sin(2*psi)*np.cos(phi0))+(2*np.cos(inc)*np.cos(2*psi)*np.sin(phi0)))
    a[j+2]=-A*(((1+(np.power(np.cos(inc),2)))*np.cos(2*psi)*np.sin(phi0))+(2*np.cos(inc)*np.sin(2*psi)*np.cos(phi0)))
    a[j+3]=-A*(((1+(np.power(np.cos(inc),2)))*np.sin(2*psi)*np.sin(phi0))-(2*np.cos(inc)*np.cos(2*psi)*np.cos(phi0)))
    j+=4

print "a_i's :", a

#no. of points
n=1000

#time over which the residuals were collected
T=np.linspace(1,500,n)

#delta time
dt=(T[1]-T[0])

#print "F_c   :", F_c_F_s.F_c(theta_p[0],phi_p[0],theta_s[0],phi_s[0])
#print "F_s   :", F_c_F_s.F_s(theta_p[0],phi_p[0],theta_s[0],phi_s[0])
#print "F_plus :", F_c_F_s.F_plus(theta_p[0],phi_p[0],theta_s[0],phi_s[0],psi)
#print "F_cros :", F_c_F_s.F_cros(theta_p[0],phi_p[0],theta_s[0],phi_s[0],psi)


print 'Generating residuals via Sesana et. al 2012'
h_p = A*(1 + (np.cos(inc)**2))
h_c = -2*A*np.cos(inc)

print "h plus : ", h_p
print "h_cros : ", h_c
n_r=np.zeros((n,m))

for j in xrange(m):
    F_plus=F_c_F_s.F_plus(theta_p[j],phi_p[j],theta_s[0],phi_s[0],psi)
    F_cros=F_c_F_s.F_cros(theta_p[j],phi_p[j],theta_s[0],phi_s[0],psi)
    for i in xrange(n):
        n_r[i,j]=-((h_p*F_plus*(np.sin((2*np.pi*f*T[i]) + phi0)-np.sin(phi0)))-(h_c*F_cros*(np.cos((2*np.pi*f*T[i]) + phi0)-np.cos(phi0))))


#r_dp=np.zeros((n,m))
#for j in xrange(m):
#    F_plus=F_c_F_s.F_plus(theta_p[j],phi_p[j],theta_s[0],phi_s[0],psi)
#    F_cros=F_c_F_s.F_cros(theta_p[j],phi_p[j],theta_s[0],phi_s[0],psi)
#    F=F_plus + 1j*F_cros
#    for i in xrange(n):
#        H=((h_p*(np.sin((2*np.pi*f*T[i]) + phi0)-np.sin(phi0)))+1j*(h_c*(np.cos((2*np.pi*f*T[i]) + phi0)-np.cos(phi0))))
#        r_dp[i,j]=np.real(H*(np.conj(F)))
#
#print r_dp
#print n_r

h1=np.zeros((n,m))
h2=np.zeros((n,m))
h3=np.zeros((n,m))
h4=np.zeros((n,m))
r=np.zeros((n,m))
res=np.zeros((n,m))

for l in xrange(n_src):
    for j in xrange(m):
        F_c=F_c_F_s.F_c(theta_p[j],phi_p[j],theta_s[0],phi_s[0])
        F_s=F_c_F_s.F_s(theta_p[j],phi_p[j],theta_s[0],phi_s[0])
        for i in xrange(n):
            h1[i,j]=F_c*np.sin(T[i])#(2*np.pi*f*T[i])
            h2[i,j]=F_s*np.sin(T[i])#(2*np.pi*f*T[i])
            h3[i,j]=F_c*np.cos(T[i])#(2*np.pi*f*T[i])
            h4[i,j]=F_s*np.cos(T[i])#(2*np.pi*f*T[i])
            r[i,j]= a[0]*h1[i,j] + a[1]*h2[i,j] + a[2]*h3[i,j] + a[3]*h4[i,j]


r_r=0.5*np.power(r,2).sum(axis=0).sum()/n
print "0.5 (r || r) = ",r_r
#print n_r,"\n"
#print r,"\n"
#print n_r-r

np.save("single.res",r)
np.save("single.tim",T)
