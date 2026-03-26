#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 10:15:37 2023

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
#import sdf_helper as sh
import h5py


################################################################################
#filename = "/home/ben/Documents/SAMS/SAMS_PI/finalState.h5"
#filename = "/home/ben/Documents/SAMS/bug_fixes/SAMS_PI/finalState.h5"
#filename='/home/ben/Documents/SAMS/MPI/diagnostics_step_93.h5'
filename='diagnostics_00010.h5'
#filename='/home/ben/Documents/SAMS/MPI/collisionless_100.h5'

f=h5py.File(filename, "r")
# Print all root level object names (aka keys) 
# these can be group or dataset names 
print("Keys: %s" % f.keys())
# get first object name/key; may or may NOT be a group
a_group_key = list(f.keys())[0]

# get the object type for a_group_key: usually group or dataset
print(type(f[a_group_key])) 

# If a_group_key is a dataset name, 
# this gets the dataset values and returns as a list
data = list(f[a_group_key])
# preferred methods to get dataset values:
#ds_obj = f[a_group_key]      # returns as a h5py dataset object
#ds_arr = f[a_group_key][()]  # returns as a numpy array
################################################################################
"""
################################################################################
# Get PIP data as comparison
filename='PIP_sod_shock_test.h5'
#filename='/home/ben/Documents/SAMS/MPI/collisionless_100.h5'

f2=h5py.File(filename, "r")
# Print all root level object names (aka keys) 
# these can be group or dataset names 
print("Keys: %s" % f2.keys())
# get first object name/key; may or may NOT be a group
a2_group_key = list(f2.keys())[0]

# get the object type for a_group_key: usually group or dataset
print(type(f2[a2_group_key])) 

# If a_group_key is a dataset name, 
# this gets the dataset values and returns as a list
data_pip = list(f2[a2_group_key])

x_pip=np.asarray(f2['xgrid'])
rho_pip=np.squeeze(np.asarray(f2['ro_p']))
vx_pip=np.squeeze(np.asarray(f2['mx_p']))/np.squeeze(np.asarray(f2['ro_p']))
pr_pip=(5.0/3.0-1.0)*(np.squeeze(np.asarray(f2['en_p'])) \
 -0.5*(np.squeeze(np.asarray(f2['bx']))**2+np.squeeze(np.asarray(f2['by']))**2+np.squeeze(np.asarray(f2['bz']))**2) \
  -0.5*(np.squeeze(np.asarray(f2['mx_p']))**2+np.squeeze(np.asarray(f2['my_p']))**2\
  +np.squeeze(np.asarray(f2['mz_p']))**2)/np.squeeze(np.asarray(f2['ro_p'])))
################################################################################

################################################################################
# Get PIP data as comparison
filename='Lare2D_sod_shock.h5'
#filename='/home/ben/Documents/SAMS/MPI/collisionless_100.h5'

f3=h5py.File(filename, "r")
# Print all root level object names (aka keys) 
# these can be group or dataset names 
print("Keys: %s" % f3.keys())
# get first object name/key; may or may NOT be a group
a3_group_key = list(f3.keys())[0]

# get the object type for a_group_key: usually group or dataset
print(type(f3[a3_group_key])) 

# If a_group_key is a dataset name, 
# this gets the dataset values and returns as a list
data_lare = list(f3[a3_group_key])

x_lare=np.asarray(f3['x'])
x_face_lare=np.asarray(f3['x_face'])
rho_lare=np.squeeze(np.asarray(f3['rho']))
vx_lare=np.squeeze(np.asarray(f3['vx']))
pr_lare=np.squeeze(np.asarray(f3['P']))
################################################################################
"""

fig, axs = plt.subplots(2, 2)

#Getting the grid is different
mesh=f.get('MeshCC')
x=mesh.get('x')[:]

rho=np.asarray(f['rho'])
vx=np.asarray(f['vx'])
en_ion=np.asarray(f['energy_ion'])
en_e=np.asarray(f['energy_electron'])
#en_electron=np.asarray(f['energy_electron'])
#en_electron=0.0
pr=en_ion*(5.0/3.0-1.0)*rho
pr_e=en_e*(5.0/3.0-1.0)*rho
T=5.0/3.0*pr/rho/2.0
axs[0,0].plot(x[0:-2],rho[0:-2,0,0],color='b',label='SAMS')
axs[0, 0].set_title('density')
axs[1,0].plot(x[0:-2],vx[0:-2,0,0],color='b')
axs[1, 0].set_title('vx')
axs[0,1].plot(x[0:-2],pr[0:-2,0,0],color='b')
axs[0, 1].set_title('pressure')
axs[1,1].plot(x[0:-2],T[0:-2,0,0],color='b')
axs[1, 1].set_title('temperature')
fig.tight_layout()
try:
    rho_n=np.asarray(f['rho_n'])
    print('rho_n')
    vx_n=np.asarray(f['vx_n'])
    print('vx_n')
    en_n=np.asarray(f['energy_neutral'])
    print('energy_n')
    pr_n=en_n*(5.0/3.0-1.0)*rho_n
    T_n=5.0/3.0*pr_n/rho_n
    axs[0,0].plot(x[0:-2],rho_n[0:-2,0,0],color='r',linestyle='--')
    axs[1,0].plot(x[0:-2],vx_n[0:-2,0,0],color='r',linestyle='--')
    axs[0,1].plot(x[0:-2],pr_n[0:-2,0,0],color='r',linestyle='--')
    axs[1,1].plot(x[0:-2],T_n[0:-2,0,0],color='r',linestyle='--')
except:
    pass

"""
axs[0,0].plot(x_pip,rho_pip,color='k',linestyle='--',label='PIP')
axs[1,0].plot(x_pip,vx_pip,color='k',linestyle='--')
axs[0,1].plot(x_pip,pr_pip,color='k',linestyle='--')
axs[1,1].plot(x_pip,5.0/3.0*pr_pip/rho_pip,color='k',linestyle='--')
#print(rho[:,0,0])
#print(rho_n[:,0,0])

axs[0,0].plot(x_lare,rho_lare,color='g',linestyle='--',label='LaRe')
axs[1,0].plot(x_face_lare,vx_lare,color='g',linestyle='--')
axs[0,1].plot(x_lare,pr_lare,color='g',linestyle='--')
axs[1,1].plot(x_lare,5.0/3.0*pr_lare/rho_lare,color='g',linestyle='--')
"""

axs[0,0].legend()
#c=plt.pcolormesh(rho[0,:,:])
#cb=plt.colorbar()
#plt.plot(rho[:,0,0],color='b')
#plt.plot(rho[:,1,1],color='b')
#plt.plot(rho[:,0,2],color='g')
#plt.plot(rho[:,0,3],color='m')
#plt.plot(rho[:,0,4],color='k')

ax = plt.gca()
#ax.set_xlim([0, 400])
#ax.set_ylim([0, 1])
plt.show()
