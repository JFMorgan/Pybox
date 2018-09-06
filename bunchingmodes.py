###############################

## Authors: Jenny Morgan 

"""
This calculates the bunching different modes as a function of z. Bunching into 
a higher harmonic can be specified in the comand line"
"""


import sys, glob, os
import numpy as np
from numpy import arange
import matplotlib.pyplot as plt
import tables
import math
from math import pi


import scipy
from scipy import special
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
from numpy.lib.scimath import *
from cmath import sin, cos
from scipy import *

import h5py

from puffdata import fdata
from puffdata import puffData
from retrieve import getPow
from retrieve import getFileSlices
from retrieve import getZData


iTemporal = 0
iPeriodic = 1


##################################################################
#
##

def getFiles(baseName): 
  """ getTimeSlices(baseName) gets a list of files

  That will be used down the line...
  """
  filelist=glob.glob(os.getcwd()+os.sep+baseName+'_electrons_*.h5')
  
  dumpStepNos=[]
  for thisFile in filelist:
    thisDump=int(thisFile.split(os.sep)[-1].split('.')[0].split('_')[-1])
    dumpStepNos.append(thisDump)

  for i in range(len(dumpStepNos)):
    filelist[i]=baseName+'_electrons_'+str(sorted(dumpStepNos)[i])+'.h5'
  return filelist


# Function to calculate the bunching into the different OAM modes #
def modebunching(filename, h):
	hf=h5py.File(filename, 'r')
	n1=hf.get('electrons')
	n1=np.array(n1)

	mdata = fdata(filename)
	rho = mdata.vars.rho  
	

	pos=n1[:,:3] 
	xpos=n1[:,0]
	ypos=n1[:,1]
	z2pos=n1[:,2]

	nume=len(pos)
	b=np.zeros([9], dtype=complex)
	phi=np.zeros([nume], dtype=complex)

	
	for i in range (0,nume):
		x=xpos[i]
		y=ypos[i]
		phi[i]=arctan2(y,x)

	for j in range (0,9):
		l=-4+j
		b[j]=1./nume*(sum(exp(1j*z2pos/(h*2.*rho))*exp(-1j*l*phi)))
	
	B=np.zeros([9], dtype=complex)

	for j in range (0,9):
		l=-4+j
		B[j]=b[j]
	return B

# Get the z position of each file #
def getZData(fname):
    h5f = tables.open_file(fname, mode='r')
    zD = h5f.root.electrons._v_attrs.zTotal
    h5f.close()
    return zD


# Get all files and put them through function to get bunching # 
baseName = sys.argv[1]

# If bunching into a harmonic is specified #
if len(sys.argv) == 3:
   harmonic = float(sys.argv[2])
else:
   harmonic = float (1.0)
filelist = getFiles(baseName)

fcount=0
B = np.zeros([len(filelist),9], dtype=complex)
zData=np.zeros([len(filelist)])
for ij in filelist:
	B[fcount]=modebunching(ij, harmonic)
	zData[fcount] = getZData(ij)
	fcount += 1





B=abs(B)

# Create a file which saves the bunching values of the l=1 and l=-1 modes #
h5f=h5py.File("bunching"+'.h5', 'w')
h5f.create_dataset('bunchingminus1', data=B[:,3])
h5f.create_dataset('bunching1', data=B[:,5])
h5f.create_dataset('zData', data=zData)
h5f.close()



font = matplotlib.font_manager.FontProperties(family='serif') # Define the font the plots will use # 

plt.figure(1)
plt.plot(zData,B[:,0], label= "l=-4")
plt.plot(zData,B[:,1], label= "l=-3")
plt.plot(zData,B[:,2], label= "l=-2")
plt.plot(zData,B[:,3], label= "l=-1")
plt.plot(zData,B[:,4], label= "l=0")
plt.plot(zData,B[:,5], label= "l=1")
plt.plot(zData,B[:,6], label= "l=2")
plt.plot(zData,B[:,7], label= "l=3") 
plt.plot(zData,B[:,8], label= "l=4")
plt.legend(prop=font)
plt.xlabel('z (m)', fontproperties=font)
plt.ylabel('|b|', fontproperties=font)
ax= plt.axes()
for label in ax.get_xticklabels():
    label.set_fontproperties(font)
for label in ax.get_yticklabels():
    label.set_fontproperties(font)
plt.savefig('Bunching at different modes for h = ' + str(harmonic) +'.png')

plt.figure(2)
plt.semilogy(zData,B[:,3], label= "l=-1")
plt.semilogy(zData,B[:,4], label= "l=0")
plt.semilogy(zData,B[:,5], label= "l=1")
plt.legend(prop=font)
plt.xlabel('z (m)', fontproperties=font)
plt.ylabel('|b|', fontproperties=font)
ax= plt.axes()
for label in ax.get_xticklabels():
    label.set_fontproperties(font)
for label in ax.get_yticklabels():
    label.set_fontproperties(font)
plt.savefig('Bunching at different modes on a log scale for h =' + str(harmonic) +'.png')

plt.figure(3)
plt.semilogy(zData,B[:,1], label= "l=-3")
plt.semilogy(zData,B[:,2], label= "l=-2")
plt.semilogy(zData,B[:,3], label= "l=-1")
plt.semilogy(zData,B[:,4], label= "l=0")
plt.semilogy(zData,B[:,5], label= "l=1")
plt.semilogy(zData,B[:,6], label= "l=2")
plt.semilogy(zData,B[:,7], label= "l=3")
plt.legend(prop=font)
plt.xlabel('z (m)', fontproperties=font)
plt.ylabel('|b|', fontproperties=font)
ax= plt.axes()
for label in ax.get_xticklabels():
    label.set_fontproperties(font)
for label in ax.get_yticklabels():
    label.set_fontproperties(font)
plt.savefig('Bunching at different modes on a log scale 7 modes for h =' + str(harmonic) + '.png')





