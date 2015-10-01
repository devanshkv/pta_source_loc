#--------------------- Spherical to cartesian coordinates --------------------------#
#  VECT=SPH_TO_CART(THETa, PHI ,RADIUS ) converts the spherical to 
#        cartesian coordinates.
#       RADIUS          :  distance from the origin
#       PHI             :  source direction grid, azimut angle from x-axis
#       THETA           :  source direction grid, polar angle from zenith
#       VECT            :  3D vector containing the cartesian coordinates
import numpy as np

#---------- This function convert the spherical to cartesian coordinates -----------#

def sph2cart(theta,phi,radius):
	x_cart=radius*np.sin(theta)*np.cos(phi)
   	y_cart=radius*np.sin(theta)*np.sin(phi)
   	z_cart=radius*np.cos(theta)
   	return np.array([x_cart,y_cart,z_cart])

#------------------------------- end of the function -------------------------------#
