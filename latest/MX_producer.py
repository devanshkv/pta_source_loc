#!/usr/bin/python

#########################################################
#                                                       #
#              Module for LLR Calculation               #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
M & X maker for Sesana et. al. 2012 eqn 22
Z maker for Haris et. al. 2014 eqn 44
'''


import numpy as np
import F_c_F_s


#normalizing factor
time=np.load("single.tim.npy")

norm_cos=0.0
norm_sin=0.0
dt=time[1]-time[0]
f=float(10e-9)
len_time=len(time)


for i in xrange(len(time)):
    norm_sin+=np.power(np.sin(time[i]),2)/len_time#(2*np.pi*f*time[i]),2)
    norm_cos+=np.power(np.cos(time[i]),2)/len_time#(2*np.pi*f*time[i]),2)

print "normalizing fator : ",norm_cos,norm_sin

#elements of inverse of M matrix
def M_maker(T,P,theta,phi,n):
    #F_c^2
    L=np.zeros((P.shape))
    #F_s^2
    M=np.zeros((P.shape))
    #F_c F_s
    N=np.zeros((P.shape))
    for i in xrange(len(theta)):
        L+=(np.power(F_c_F_s.F_c(theta[i],phi[i],T,P),2)*norm_sin)
        M+=(np.power(F_c_F_s.F_s(theta[i],phi[i],T,P),2)*norm_sin)
        N+=((F_c_F_s.F_c(theta[i],phi[i],T,P)*F_c_F_s.F_s(theta[i],phi[i],T,P))*norm_sin)
    Mat_test=np.zeros((4*n,4*n,len(T),len(T)/2 +1))
    for i in xrange(4*n-1):
            if(Mat_test[i,i].all() == 0.0):
                    Mat_test[i,i]=(-M/(N*N - L*M))
                    Mat_test[i+1,i+1]=(-L/(N*N - L*M))
            if(i%2==0):
                    Mat_test[i,i+1]=( N/(N*N - L*M))
                    Mat_test[i+1,i]=( N/(N*N - L*M))
    return(Mat_test)

#X vector
def X_maker(T,P,theta,phi,n):
    res=np.load("single.res.npy")
    time=np.load("single.tim.npy")
    f=float(10e-9)
    y=len(theta)
    z=len(T)
    x=np.zeros((4*n,z,z/2 + 1))
    for k in xrange(n):
        for j in xrange(y):
            F_c=F_c_F_s.F_c(theta[j],phi[j],T,P)
            F_s=F_c_F_s.F_s(theta[j],phi[j],T,P)
            for i in xrange(len(time)):
                x[k+0]+=res[i,j]*F_c*np.sin(time[i])/len_time#(2*np.pi*f*time[i])
                x[k+1]+=res[i,j]*F_s*np.sin(time[i])/len_time#(2*np.pi*f*time[i])
                x[k+2]+=res[i,j]*F_c*np.cos(time[i])/len_time#(2*np.pi*f*time[i])
                x[k+3]+=res[i,j]*F_s*np.cos(time[i])/len_time#(2*np.pi*f*time[i])
    return(x)


# Calculates everthing needed for Haris et. al 2014 eqn 44
def Z_maker(T,P,theta,phi,n):
    #Mod of F_+,x DP
    Mod_F_DP_s=F_c_F_s.Mod_F_DP_s(theta,phi,T,P)
    Mod_F_DP_c=F_c_F_s.Mod_F_DP_c(theta,phi,T,P)
    res=np.load("single.res.npy")
    time=np.load("single.tim.npy")
    #no. of pulsars
    y=len(theta)
    #no, of samples
    x=len(T)
    #4 quantities as in Haris et. al 2014 eqn 44
    z1=np.zeros((P.shape))
    z2=np.zeros((P.shape))
    z3=np.zeros((P.shape))
    z4=np.zeros((P.shape))
    #(\rho_L,R)^2
    zL=np.zeros((P.shape))
    zR=np.zeros((P.shape))
    #needed for F_+,x DP
    denL=np.zeros((P.shape))
    denR=np.zeros((P.shape))
    #two phases
    pL=np.zeros((P.shape))
    pR=np.zeros((P.shape))
    DPtheta=F_c_F_s.DP_theta(theta,phi,T,P)
    print "Computing Z"
    for k in xrange(n):
        for j in xrange(y):
            F_DP_c=F_c_F_s.F_plus_DP(theta[j],phi[j],T,P,DPtheta)
            F_DP_s=F_c_F_s.F_cros_DP(theta[j],phi[j],T,P,DPtheta)
            for i in xrange(len(time)):
                z1 += res[i,j]*F_DP_c*np.cos(time[i])/len_time#/norm_cos
                z2 += res[i,j]*F_DP_c*np.sin(time[i])/len_time#/norm_sin
                z3 += res[i,j]*F_DP_s*np.cos(time[i])/len_time#/norm_cos
                z4 += res[i,j]*F_DP_s*np.sin(time[i])/len_time#/norm_sin
        denL =np.power(Mod_F_DP_c,2)
        denR =np.power(Mod_F_DP_s,2)
        zL+=np.divide(np.power(z1,2)+np.power(z2,2),(2.0*denL*norm_sin))
        zR+=np.divide(np.power(z3,2)+np.power(z4,2),(2.0*denR*norm_cos))
        pL+=-np.arctan2(z1,z2)
        pR+=-np.arctan2(z3,z4)
        return(zL+zR,zL,zR,pL,pR)
