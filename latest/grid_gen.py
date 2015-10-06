#!/usr/bin/python

import numpy as np

#number of grid points
x=100
# Theta phi-grid
t=np.linspace(0, np.pi, num=x)
p=np.linspace(0,2*np.pi,num=2*x)

np.save("theta.grid",t)
np.save("phi.grid",p)

