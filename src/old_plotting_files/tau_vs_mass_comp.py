# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:30:07 2015

@author: hunterbrooks

plot tau vs mass difference between mine and Alexs

"""


import numpy as np
#folderOfInterest = '/users/hunterbrooks/astronomy/results/1.4c_new_no_neb_emis/'
folderOfInterest = '/users/hunterbrooks/astronomy/results/1.4c/'

my_mass_idx = 4
his_mass_idx = 2

info = np.loadtxt(folderOfInterest+'superInfo.dat') 
predictions = np.load(folderOfInterest+'predictions.dat.npy')   
errors = np.load(folderOfInterest+'errors.dat.npy')
myField = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/superInfo.dat')[:,3]
myID = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/superInfo.dat')[:,4]

alexs_data = np.genfromtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat')
aField = alexs_data[:,7]
aID = alexs_data[:,8]

alex_sorted = np.zeros(len(predictions[:,0]))
alex_sorted_e = np.zeros(len(predictions[:,0]))

def matchbyfield(infield, inid):
    row_index = -1
    for i in range(0, np.shape(aID)[0]):
        if infield == aField[i]:
            if inid == aID[i]:
                row_index = i
        
    if row_index == -1:
        return False
    
    return row_index
    
for ii in range(0, len(myField)):
    his_idx = matchbyfield(myField[ii], myID[ii])
    if his_idx != False:
        print myID[ii], aID[his_idx]
        alex_sorted[ii] = alexs_data[his_idx, his_mass_idx]
        alex_sorted_e[ii] = alexs_data[his_idx, his_mass_idx+1]
        
import matplotlib.pylab as plt


quad_sum_errs = np.sqrt( (alex_sorted[:]/alex_sorted_e[:])**2 + (predictions[:,my_mass_idx]/errors[:,my_mass_idx,0])**2   )
delta_mass = predictions[:,my_mass_idx] - alex_sorted[:]


"""
# THIS IS FOR MASS DIFFERENCE
plt.errorbar(predictions[:,my_mass_idx],
             alex_sorted[:],
                xerr=[np.abs(errors[:,my_mass_idx,0]-predictions[:,my_mass_idx]), np.abs(errors[:,my_mass_idx,1]-predictions[:,my_mass_idx])],
                yerr=alex_sorted_e[:],
                fmt='bo')
                
onetoone_x = np.linspace(np.min(predictions[:,my_mass_idx]),
                         np.max(predictions[:,my_mass_idx]),10)
                         
plt.plot(onetoone_x,onetoone_x,'g--')

plt.plot(onetoone_x,onetoone_x,'r--')


"""

"""
# error bar test
plt.errorbar(delta_mass,
             predictions[:,0],
             xerr=quad_sum_errs,
             yerr=[np.abs(errors[:,0,0]-predictions[:,0]), np.abs(errors[:,0,1]-predictions[:,0])])
"""

plt.scatter(predictions[:,0],delta_mass)

#plt.title("Mass Difference (w/o neb lines or con't) vs Tau")
plt.title("Mass Difference vs Tau")
plt.ylabel("MCSED (with nebular stuffs) mass - Alex mass [M_sol]")
#plt.ylabel("MCSED  (w/o lines or con't) mass - Alex mass [M_sol]")
plt.xlabel('Tau')
             


plt.show()
