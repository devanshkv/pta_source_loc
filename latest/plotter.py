#!/usr/bin/python

#########################################################
#                                                       #
#           Plots the |F|^2 map for given pulsars       #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################
'''
Following program creates |F|^2 map for given pulsars 
coordinates in the 'pulsar.loc' file and save the map
as multipage.pdf. Number of points in the grid can be
changed by changing x.
Parameters of the countour plot can be changed in the
variable lvls.
'''

import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from math import *
from sph_cart import *
import numpy as np

def antenna_plot_3D(mod_F,T,P,theta,phi,name):
        [x,y,z]=sph2cart(T,P,mod_F)
	[x_p,y_p,z_p]=sph2cart(theta,phi,(0.2+max(mod_F.flatten()))*np.ones(len(theta)))
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	p_surf=ax.plot_surface(x,y,z,rstride=1,cstride=1,linewidth=0,antialiased=True,facecolors=cm.jet(np.sqrt(x*x + y*y + z*z)))
	ax.scatter(x_p,y_p,z_p,c='g',s=30,marker='*')
	m = cm.ScalarMappable(cmap=cm.jet)
	m.set_array(x*x + y*y + z*z)
	plt.colorbar(m)
	ax.set_zlabel('Z axis',fontsize=15)
	ax.set_xlabel('X axis',fontsize=15)
	ax.set_ylabel('Y axis',fontsize=15)
	#ax.set_ylim([-.6,1.2])
	#ax.set_xlim([-1.2,1.2])
	#ax.set_zlim([-1.2,1.2])
	plt.savefig(name,bbox_inches='tight')
	plt.show()
	return plt.show()

def antenna_plot_contour(NetworkF,Theta,Phi,theta,phi,name):
	fig,ax = plt.subplots()
	levels = np.linspace(0,1,100)
	plt.contour(Phi,Theta,NetworkF,colors='k')
	levels = np.linspace(min(NetworkF.flatten()), max(NetworkF.flatten()), 100)
	p= plt.contourf(Phi,Theta,NetworkF,levels=levels)
	plt.plot(phi,theta,'g*',markersize=10)
	plt.xlabel(r'$\phi$'' (radians)')
	plt.ylabel(r'$\theta$'' (radians)')
	plt.colorbar(p,orientation='vertical')
	plt.savefig(name,bbox_inches='tight')
	return plt.show()

def dot_plot_contour(NetworkF,Theta,Phi,theta,phi,th,ph,name):
        fig,ax = plt.subplots()
        levels = np.linspace(0,1,100)
        plt.contour(Phi,Theta,NetworkF,colors='k')
        levels = np.linspace(min(NetworkF.flatten()), max(NetworkF.flatten()), 100)
        p= plt.contourf(Phi,Theta,NetworkF,levels=levels)
        plt.plot(phi,theta,'g*',markersize=10)
	plt.plot(ph,th,'ko',markersize=10)
        plt.xlabel(r'$\phi$'' (radians)')
        plt.ylabel(r'$\theta$'' (radians)')
        plt.colorbar(p,orientation='vertical')
        plt.savefig(name,bbox_inches='tight')
        return plt.show()

def llr_contour(NetworkF,Theta,Phi,theta,phi,t_sl,p_sl,t_sa,p_sa,name):
    fig,ax = plt.subplots()
    levels = np.linspace(0,1,100)
    plt.contour(Phi,Theta,NetworkF,colors='k')
    levels = np.linspace(min(NetworkF.flatten()), max(NetworkF.flatten()), 100)
    p= plt.contourf(Phi,Theta,NetworkF,levels=levels)
    plt.plot(phi,theta,'g*',markersize=10)
    plt.plot(p_sa,t_sa,'w*',markersize=10)
    plt.plot(p_sl,t_sl,'c*',markersize=10)
    plt.xlabel(r'$\phi$'' (radians)')
    plt.ylabel(r'$\theta$'' (radians)')
    plt.colorbar(p,orientation='vertical')
    plt.savefig(name,bbox_inches='tight')
    return plt.show()

def llr_contour_no_show(NetworkF,Theta,Phi,theta,phi,t_sl,p_sl,t_sa,p_sa,name):
        mpl.use('Agg') 
	fig,ax = plt.subplots()
        levels = np.linspace(0,1,100)
        plt.contour(Phi,Theta,NetworkF,colors='k')
        levels = np.linspace(min(NetworkF.flatten()), max(NetworkF.flatten()), 100)
        p= plt.contourf(Phi,Theta,NetworkF,levels=levels)
        plt.plot(phi,theta,'g*',markersize=10)
        plt.plot(p_sa,t_sa,'w*',markersize=10)
        plt.plot(p_sl,t_sl,'c*',markersize=10)
        plt.xlabel(r'$\phi$'' (radians)')
        plt.ylabel(r'$\theta$'' (radians)')
        plt.colorbar(p,orientation='vertical')
        plt.savefig(name,bbox_inches='tight')
#        return plt.show()


def hammer_plot(NetworkF,Theta,Phi,theta,phi):
        from mpl_toolkits.basemap import Basemap, shiftgrid
        fig,ax = plt.subplots()
        RAD = 180/np.pi
        m = Basemap(projection='hammer',lon_0=0,resolution='c')
	Phi, NetworkF = shiftgrid(180,Phi,NetworkF)
        m.contourf(Phi*RAD, Theta*RAD, NetworkF, 100, cmap=plt.cm.jet,latlon=True)
        plt.plot(phi,theta,'go',markersize=5)
        plt.xlabel(r'$\phi$'' (radians)')
        plt.ylabel(r'$\theta$'' (radians)')
        plt.savefig('hammer.png')
        return plt.show()

