# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:52:12 2015

@author: hunterbrooks

Pull data from the 18 column catalog created in interpret_results.py and plot
the masses in that catalog verse the masses in Alex's catalog.

ALEX: youll need to replace the path for variables "hunters_cata" and 
"alexs_data"

"""

import numpy as np

# load necesary parts of data from hunters catalog
hunters_cata = np.loadtxt('/users/hunterbrooks/astronomy/results/1.4c/catalog.txt')
hmass = hunters_cata[:,16]
hmass_errlow = hunters_cata[:,15]
hmass_errhi = hunters_cata[:,17]
hfield = hunters_cata[:, 0]
hid = hunters_cata[:,1]

# load neccesary parts of alexs catalog
# ALEXS COLUMN NAMES:
# RA Dec stellarmass smerr logsfr_cor logsfr_cor_err LAE field id beta betaerr z 12+log(O/H) L_UV half_light_pix ellipticity ellipticity_err "nn distance pix" half_light_kpc "nn distance kpc" hmagcut lyalum lyalumnegerr lyalumposerr sfi ssfr ssfrerr
alexs_data = np.genfromtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat')
amass = alexs_data[:,2]
amass_errlow = alexs_data[:,3] #alex only has 1 mass error value, use it for 
amass_errhi = alexs_data[:,3] # both high error and low error
afield = alexs_data[:,7]
aid = alexs_data[:,8]

# start blank lists that we'll append matching objects too
plt_hmass = []
plt_hlowerr = []
plt_hhierr = []
h_id_matches = [] #used to confirm we matched right obj ID/field
h_field_matches = [] #used to confirm we matched right obj ID/field
plt_amass = []
plt_alowerr = []
plt_ahierr = []
a_id_matches = [] #used to confirm we matched right obj ID/field
a_field_matches = [] #used to confirm we matched right obj ID/field


# loop through all hunters objects, attempting to find a match in alex's catalog
for i in range(0, len(hmass)):
    #get details on the object from my catalog we're working with right now
    current_hid = int(hid[i])
    current_hfield = int(hfield[i])
    match = False #variable for the index of matching obj in alex's catalog
    # loop through all alex's objects attempting to find a match
    for j in range(0, len(amass)):
        current_aid = int(aid[j])
        current_afield = int(afield[j])
        if current_hid == current_aid:
            if current_hfield == current_afield:
                # if the above 2 conditions are satisfied, there is a match
                plt_hmass.append(hmass[i])
                plt_hlowerr.append(hmass_errlow[i])
                plt_hhierr.append(hmass_errhi[i])
                
                plt_amass.append(amass[j])
                plt_alowerr.append(amass_errlow[j])
                plt_ahierr.append(amass_errhi[j])
                
                h_id_matches.append(hid[i])
                h_field_matches.append(hfield[i])
                
                a_id_matches.append(aid[j])
                a_field_matches.append(afield[j])
            

print plt_amass
print

     
import matplotlib.pylab as plt

plt.title('Alex vs. MCSED Mass Comparison')
plt.xlabel('log(MCSED mass) [M_sol] (WITH NEB EMIS)')
plt.ylabel('log(ALEX mass) [M_sol]')

plt.errorbar(   plt_hmass, plt_amass,
                xerr=[plt_hlowerr,plt_hhierr], yerr=[plt_alowerr, plt_ahierr],
                fmt='bo')
                
# make a 1 to 1 line
one_to_one_range = np.linspace(np.min(plt_hmass),np.max(plt_hmass),2)
plt.plot(one_to_one_range, one_to_one_range, 'r--')
plt.show()
