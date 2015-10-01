#!/usr/bin/python

#########################################################
#                                                       #
#      		F_c and F_s for anntenna response       #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
The function calculates the value of F_c and F_s for given
coordinates in the map and a pulsar.
Variables : 
theta,phi = for the pulsar
t,p = for the grid
'''
#libraries
import numpy as np

def F_s(theta,phi,t,p):
	return(((np.sin(2*theta)*np.sin(t)*np.sin(p-phi))+(np.power(np.sin(theta),2)*np.cos(t)*np.sin(2*(phi-p))))/(2.0*(1.000001-(np.cos(theta)*np.cos(t)+(np.cos(phi-p)*np.sin(t)*np.sin(theta))))))

def F_c(theta,phi,t,p):
	return((-(np.power((np.sin(t)),2)*(np.power(np.sin(theta),2)-2*np.power(np.cos(theta),2)))-(np.sin(2*t)*np.sin(2*theta)*np.cos(p-phi))+((1+np.power(np.cos(t),2))*np.cos(2*(p-phi))*np.power(np.sin(theta),2)))/(4.0*(1.0000001-(np.cos(theta)*np.cos(t)+(np.cos(phi-p)*np.sin(t)*np.sin(theta))))))

#def DP_theta(theta,phi,t,p):
#	Num=np.zeros((p.shape))
#        Den=np.zeros((p.shape))
#	res=np.zeros((p.shape))
#        for i in xrange(len(theta)):
#		Num+=np.multiply(F_s(theta[i],phi[i],t,p),F_c(theta[i],phi[i],t,p))
#		Den+=np.power(F_c(theta[i],phi[i],t,p),2) - np.power(F_s(theta[i],phi[i],t,p),2)
#	res=0.50 * np.arctan(2.0*Num/Den)
#	return(res)

def DP_theta(theta,phi,t,p):
	for_arg=np.zeros((p.shape))
	for i in xrange(len(theta)):
		d_alpha = 1j*F_s(theta[i],phi[i],t,p) + F_c(theta[i],phi[i],t,p)
		for_arg = for_arg + d_alpha*d_alpha
	return(0.5*np.angle(for_arg))

def F_DP_s(theta,phi,t,p,Theta,Phi):
	return(np.multiply(F_s(theta,phi,t,p), np.cos(DP_theta(Theta,Phi,t,p)))-np.multiply(F_c(theta,phi,t,p), np.sin(DP_theta(Theta,Phi,t,p))))

def F_DP_c(theta,phi,t,p,Theta,Phi):
	return(np.multiply(F_c(theta,phi,t,p), np.cos(DP_theta(Theta,Phi,t,p)))+np.multiply(F_s(theta,phi,t,p), np.sin(DP_theta(Theta,Phi,t,p))))

def Mod_F_DP_s(theta,phi,t,p):
        Mod_F_DP_s=np.zeros((p.shape))
        for i in xrange(len(theta)):
               Mod_F_DP_s += (np.power(F_DP_s(theta[i],phi[i],t,p,theta,phi),2))
        return(np.sqrt(Mod_F_DP_s))

def Mod_F_DP_c(theta,phi,t,p):
        Mod_F_DP_c=np.zeros((p.shape))
        for i in xrange(len(theta)):
               Mod_F_DP_c += (np.power(F_DP_c(theta[i],phi[i],t,p,theta,phi),2))
        return(np.sqrt(Mod_F_DP_c))

def F_DP_dot(theta,phi,t,p):
	dot=np.zeros((p.shape))
	for i in xrange(len(theta)):
		dot += np.multiply(F_DP_c(theta[i],phi[i],t,p,theta,phi),F_DP_s(theta[i],phi[i],t,p,theta,phi))
	return(dot)

def F_dot(theta,phi,t,p):
	dot=np.zeros((p.shape))
        for i in xrange(len(theta)):
                dot += np.multiply(F_c(theta[i],phi[i],t,p),F_s(theta[i],phi[i],t,p))
        return(dot)

'''
Variables :
t,p = for the pulsar
theta,phi = for the grid
'''

def F_p_p(t,p,theta,phi,psi):
	q_g=np.array([(np.sin(theta)*np.cos(psi) + np.sin(psi)*np.cos(phi)*np.cos(theta)),(-np.cos(phi)*np.cos(psi)+np.sin(psi)*np.sin(phi)*np.cos(theta)),-(np.sin(psi)*np.sin(theta))])
	p_g=np.array([(-np.sin(phi)*np.sin(psi) + np.cos(psi)*np.cos(phi)*np.cos(theta)),(np.cos(phi)*np.sin(psi)+np.sin(phi)*np.cos(psi)*np.cos(theta)),-(np.cos(psi)*np.sin(theta))])
	n=np.array([np.sin(t)*np.cos(p),-np.sin(t)*np.sin(p), np.cos(t)])
	k=np.array([-np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), -np.cos(theta)])
	a=np.inner(np.transpose(k),n)
	return(np.transpose((np.power(np.inner(np.transpose(p_g),n),2)-np.power(np.inner(np.transpose(q_g),n),2))/(2.0*(1.0000001*np.ones(a.shape)+a))))

def F_c_p(t,p,theta,phi,psi):
    q_g=np.array([(np.sin(theta)*np.cos(psi) + np.sin(psi)*np.cos(phi)*np.cos(theta)),(-np.cos(phi)*np.cos(psi)+np.sin(psi)*np.sin(phi)*np.cos(theta)),-(np.sin(psi)*np.sin(theta))])
    p_g=np.array([(-np.sin(phi)*np.sin(psi) + np.cos(psi)*np.cos(phi)*np.cos(theta)),(np.cos(phi)*np.sin(psi)+np.sin(phi)*np.cos(psi)*np.cos(theta)),-(np.cos(psi)*np.sin(theta))])
    n=np.array([np.sin(t)*np.cos(p),-np.sin(t)*np.sin(p), np.cos(t)])
    k=np.array([-np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), -np.cos(theta)])
    a=np.inner(np.transpose(k),n)
    return(np.transpose(np.inner(np.transpose(p_g),n)*np.inner(np.transpose(q_g),n)/(1.0000001*np.ones(a.shape)+a)))

def F_c_pc(t,p,theta,phi,psi):
    return(F_p_p(t,p,theta,phi,psi)*np.cos(2*psi) + F_c_p(t,p,theta,phi,psi)*np.sin(2*psi))
def F_s_pc(t,p,theta,phi,psi):
    return(-F_c_p(t,p,theta,phi,psi)*np.cos(2*psi) + F_p_p(t,p,theta,phi,psi)*np.sin(2*psi))

def mod_F(t,p,theta,phi):
	n=np.array([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)])
        k=np.array([-np.sin(theta)*np.cos(phi), -np.sin(theta)*np.sin(phi), -np.cos(theta)])
        a=np.inner(np.transpose(k),n)
	return((1.0-np.transpose(np.inner(np.transpose(n),k))**2)/(1.0000001*np.ones(a.shape)+a))

#def F_plus(t_p,p_p,t_g,p_g,psi):
#    u=np.array([(np.cos(t_g)*np.cos(p_g)),(np.cos(t_g)*np.sin(p_g)),(-np.sin(t_g))])
#    v=np.array([(-np.sin(p_g)),(np.cos(p_g)),(-np.sin(t_g))])
#    p=(u*np.cos(psi)) - v*np.sin(psi)
#    q=(v*np.cos(psi)) + u*np.sin(psi)
#    return(np.transpose((np.power(np.inner(np.transpose(p_g),n),2)-np.power(np.inner(np.transpose(q_g),n),2))/(2.0*(1.0000001*np.ones(a.shape)+a))))
#
#def F_cross(t_p,p_p,t_g,p_g,psi):
#    u=np.array([(np.cos(t_g)*np.cos(p_g)),(np.cos(t_g)*np.sin(p_g)),(-np.sin(t_g))])
#    v=np.array([(-np.sin(p_g)),(np.cos(p_g)),(-np.sin(t_g))])
#    p=(u*np.cos(psi)) - v*np.sin(psi)
#    q=(v*np.cos(psi)) + u*np.sin(psi)
#    return(np.transpose(np.inner(np.transpose(p_g),n)*np.inner(np.transpose(q_g),n)/(1.0000001*np.ones(a.shape)+a))
