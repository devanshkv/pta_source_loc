#########################################################
#                                                       #
#           F_c and F_s for anntenna response           #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
The function calculates the value of F_c and F_s for given
coordinates in the map and a pulsar.
Variables :.
t_p,p_p = for the pulsar
t_s,p_s = for the grid
'''

#libraries
from numpy import *

def F_plus(t_p,p_p,t_s,p_s,psi):
    u=array([cos(t_s)*cos(p_s),cos(t_s)*sin(p_s),-sin(t_s)])
    v=array([sin(p_s),-cos(p_s),zeros(p_s.shape)])
    n=array([sin(t_p)*cos(p_p),sin(t_p)*sin(p_p),cos(t_p)])
    k=array([-sin(t_s)*cos(p_s),-sin(t_s)*sin(p_s),-cos(t_s)])
    p=(u*cos(psi))-(v*sin(psi))
    q=(u*sin(psi))+(v*cos(psi))
    a=inner(transpose(k),n)
    return(transpose((power(inner(transpose(p),n),2)-power(inner(transpose(q),n),2))/(2.0*(1.0000001*ones(a.shape)+a))))

def F_cros(t_p,p_p,t_s,p_s,psi):
    u=array([cos(t_s)*cos(p_s),cos(t_s)*sin(p_s),-sin(t_s)])
    v=array([sin(p_s),-cos(p_s),zeros(p_s.shape)])
    n=array([sin(t_p)*cos(p_p),sin(t_p)*sin(p_p),cos(t_p)])
    k=array([-sin(t_s)*cos(p_s),-sin(t_s)*sin(p_s),-cos(t_s)])
    p=(u*cos(psi))-(v*sin(psi))
    q=(u*sin(psi))+(v*cos(psi))
    a=inner(transpose(k),n)
    return(transpose(inner(transpose(p),n)*inner(transpose(q),n)/(1.0000001*ones(a.shape)+a)))

def F_c(t_p,p_p,t_s,p_s):
    u=array([cos(t_s)*cos(p_s),cos(t_s)*sin(p_s),-sin(t_s)])
    v=array([sin(p_s),-cos(p_s),zeros(p_s.shape)])
    n=array([sin(t_p)*cos(p_p),sin(t_p)*sin(p_p),cos(t_p)])
    k=array([-sin(t_s)*cos(p_s),-sin(t_s)*sin(p_s),-cos(t_s)])
    a=inner(transpose(k),n)
    res=transpose((power(inner(transpose(u),n),2)-power(inner(transpose(v),n),2))/(2.0*(ones(a.shape)+a)))
    return(res)

def F_s(t_p,p_p,t_s,p_s):
    u=array([cos(t_s)*cos(p_s),cos(t_s)*sin(p_s),-sin(t_s)])
    v=array([sin(p_s),-cos(p_s),zeros(p_s.shape)])
    n=array([sin(t_p)*cos(p_p),sin(t_p)*sin(p_p),cos(t_p)])
    k=array([-sin(t_s)*cos(p_s),-sin(t_s)*sin(p_s),-cos(t_s)])
    a=inner(transpose(k),n)
    return(-transpose(inner(transpose(u),n)*inner(transpose(v),n)/(ones(a.shape)+a)))

def DP_theta(t_p,p_p,t_s,p_s):
    for_arg=zeros((p_s.shape))
    for i in xrange(len(t_p)):
        F_cdp=F_c(t_p[i],p_p[i],t_s,p_s)
        F_sdp=F_s(t_p[i],p_p[i],t_s,p_s)
        d_alpha = 1j*F_sdp + F_cdp
        for_arg = for_arg + multiply(d_alpha,d_alpha)
    return(arctan2(for_arg.imag,for_arg.real))

def F_plus_DP(t_p,p_p,t_s,p_s,theta):
    return(multiply(F_c(t_p,p_p,t_s,p_s),cos(theta/2)) + multiply(F_s(t_p,p_p,t_s,p_s),sin(theta/2)))

def F_cros_DP(t_p,p_p,t_s,p_s,theta):
    return(-multiply(F_c(t_p,p_p,t_s,p_s),sin(theta/2)) + multiply(F_s(t_p,p_p,t_s,p_s),cos(theta/2)))

def Mod_F_DP_c(t_p,p_p,t_s,p_s):
    dot=zeros((p_s.shape))
    dptheta=DP_theta(t_p,p_p,t_s,p_s)
    for i in xrange(len(t_p)):
        dot += (power(F_plus_DP(t_p[i],p_p[i],t_s,p_s,dptheta),2))
    return(sqrt(dot))

def Mod_F_DP_s(t_p,p_p,t_s,p_s):
    dot=zeros((p_s.shape))
    dptheta=DP_theta(t_p,p_p,t_s,p_s)
    for i in xrange(len(t_p)):
        dot += (power(F_cros_DP(t_p[i],p_p[i],t_s,p_s,dptheta),2))
    return(sqrt(dot))

def DP_theta_new(theta,phi,t,p):
    Num=zeros((p.shape))
    Den=zeros((p.shape))
    res=zeros((p.shape))
    for i in xrange(len(theta)):
            Num+=multiply(F_s(theta[i],phi[i],t,p),F_c(theta[i],phi[i],t,p))
            Den+=power(F_c(theta[i],phi[i],t,p),2) - power(F_s(theta[i],phi[i],t,p),2)
    res=arctan2(2.0*Num,Den)
    return(res)
def F_dot(theta,phi,t,p):
    dot=zeros((p.shape))
    dptheta=DP_theta(theta,phi,t,p)
    for i in xrange(len(theta)):
        dot += multiply(F_plus_DP(theta[i],phi[i],t,p,dptheta),F_cros_DP(theta[i],phi[i],t,p,dptheta))
    return(dot)
