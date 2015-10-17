#!/usr/bin/python

#########################################################
#                                                       #
#           Places Pulsars On The Grid Randomly         #
#               Written by : Devansh Agarwal            #
#               devansh@iisertvm.ac.in                  #
#                                                       #
#########################################################

#libraries
import random
import math
import numpy
import pylab
import sys

#input for no of pulsars you nee to place.
n = int(sys.argv[1])

#writes the locations in a file.
f = open('pulsar.loc', 'w')

phi=numpy.zeros(n)
theta=numpy.zeros(n)

#seed
for i in xrange(n):
	phi[i] = random.uniform(0,2*math.pi)
	theta[i] = random.uniform(0, math.pi)
	print  >> f, theta[i], phi[i]

#the plot
#pylab.plot(phi,theta, 'ro')
#pylab.grid()
#pylab.ylim([0,math.pi])
#pylab.xlim([0,2*math.pi])
#pylab.xlabel('$\\theta$')
#pylab.ylabel('$\phi$')
#pylab.show()
