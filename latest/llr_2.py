#!/usr/bin/python

#########################################################
#                                                       #
#              log likelihood calculator                #
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

#number of grid points
x=100

#read the pulsar locations
f= open('pulsar.loc','r')

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

f.close()
print 'Pulsar locations read'

#no of pulsars
n=len(theta)


# Theta phi-grid
t=np.linspace(0, np.pi, num=x)
p=np.linspace(0,2*np.pi,num=2*x)

g=open('source.loc','r')

s_theta=[]
s_phi=[]

for line in g:
                        #Read and Split
                        line = line.strip()
                        columns = line.split()
                        #Read pulsar coordinates and put them in theta and phi lists.
                        s_theta.append(float(columns[0]))
                        s_phi.append(float(columns[1]))

#print 'Pulsar locations read'


s_theta=np.array(s_theta)
s_phi=np.array(s_phi)

t=np.sort(np.append(t,s_theta))
p=np.sort(np.append(p,s_phi))

#to sum F_c and F_s over all pulsars

[T,P]=np.meshgrid(t,p)
print 'Calculating Log Likelihood pattern'

#calculate X_j via Sesana et. al 2012 eqn 22
x=MX_producer.X_maker(T,P,theta,phi,1)

#calculate inverse of M_jk via Sesana et. al 2012 eqn 22
Mat=MX_producer.M_maker(T,P,theta,phi,1)

print Mat.shape,x.shape

#maximized log likelihood and revoered a_is's
result= 0.5*np.einsum('inm,ijnm,jnm->nm', x, Mat, x)
recovery= np.einsum('ijnm,jnm->inm', Mat, x)
#
#
i,j = np.unravel_index(result.argmax(), result.shape)
print "Max LLR : ",result[i,j]
print "Source  : ",T[i,j],P[i,j]
print "a_i's   : ",recovery[:,i,j]
plotter.llr_contour(result,T,P,theta,phi,T[i,j],P[i,j],s_theta,s_phi,"llr_2")
print "Log Likelihood map saved as \'llr_2.png\'"
