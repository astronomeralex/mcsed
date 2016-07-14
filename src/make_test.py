# -*- coding: utf-8 -*-
"""
Returns average fluxes [cgs/Hz] for nEntires fake objects.  Each fake object is 
created from parameters as defined in variable_house.py.  There is an 
option to create them from random parameters values. Redshift and error
is constant for each object.

Edit within the for loop if you wish to use the same parameter values for each
nEntires. requires random_truths = False.
"""
# =========================================================================== #
import numpy as np
import factory
import math_functions as mathf
import variable_house as vHouse
# get home
import wheres_home
home = wheres_home.getHomeLocation()

# =========================================================================== #
# get parameters
nDim = vHouse.nDim
paramRanges = vHouse.paramRanges
paramNames = vHouse.nameL
if vHouse.walker_SED_type == 'Interp':
    pSpaceGridSize = vHouse.pSpaceGridSize

# =========================================================================== #
# find nFilters
col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')
nFilters = len(col_names)
    
# =========================================================================== #
# create output arrays, shape is descriped in documentation
nEntires = 1
outputFluxes = np.zeros([nEntires, nFilters])
outputErrors = np.zeros([nEntires, nFilters])
outputInfo = np.zeros([nEntires, 10]) #dummy array that keeps main.py happy
outputRedshift = np.zeros([nEntires])
fake_err = 0.10 # fake error of measurement in psuedo-SED
random_truths = True # would you like to randomly select truths?
z = 2.0 # redshift, constant for ALL nEntries
outputTruths = np.zeros([nEntires,nDim]) 
outputMasses = np.zeros([nEntires,2])

# =========================================================================== #
# start loop, making nEntires fake objects
for run_i in range(0, nEntires):
    
    # ======================================================================= #
    # set fake data properties
    if random_truths == False:
        # edit here if random_truths = False!
        fake_params = np.array([1.3, 0.45, 0.2, -0.75])

    if random_truths == True:
        fake_params = np.zeros([vHouse.nDim])
        for pi in range(0,vHouse.nDim):
            fake_params[pi] = np.random.uniform(vHouse.lowerL[pi],
                                              vHouse.higherL[pi],
                                              1)
        
    # ======================================================================= #
    # create fake stellar population
    sp, s_nu, nu_obs = factory.fakeSED(    paramNames, fake_params,
                                            fake_err,
                                            z)
    lambda_e = sp.wavelengths
                         
    outputMasses[run_i, 0], outputMasses[run_i, 1] = sp.stellar_mass, sp.dust_mass
    # ======================================================================= #
    # create an R_nu array, dont need to do this each time but it's a fast
    # process so its really no big deal. Also it keeps things simple
    R_nu, nFilters = factory.R_nu(home,
                                  nu_obs)
    # ======================================================================== #
    # create D_nu array
    D_nu = factory.D_nu(nFilters,
                        nu_obs)
    
    # ======================================================================= #  
    # create NU_obs array
    NU_obs = factory.NU_obs(nFilters,
                            nu_obs)
                           
    # ======================================================================= # 
    # create S_nu array
    S_nu = factory.S_nu(nFilters,
                        s_nu)
                          
    # ======================================================================= #
    # get avg fluxes for fake data
    dataAvgFlux = mathf.getAvgFlux(   nFilters,
                                   
                                      S_nu,
                                      NU_obs,
                                      D_nu,
                                      R_nu,
                                      True)
                                      
    dataAvgFluxErr = fake_err * dataAvgFlux
        
    # ======================================================================= #
    # save the test data in super.dat format
    fake_params = np.reshape(fake_params, [nDim])
    outputTruths[run_i,:] = fake_params
    
    # set redshift of object as first column
    outputRedshift[run_i] = z
    
    for filt_i in range(0, nFilters):
        outputFluxes[run_i, filt_i  ] = dataAvgFlux[filt_i]
        outputErrors[run_i, filt_i  ] = dataAvgFluxErr[filt_i]
        
# =========================================================================== #
# Save outputSuperDat and outputTruths
print
print ' --- This will over right superFluxes, superErrors, superInfo, truths,'
print ' --- and masses .dat in /data/active/.  Press enter to continue.'
raw_input()
np.savetxt(home+'data/active/superFluxes.dat',outputFluxes)
np.savetxt(home+'data/active/superErrors.dat',outputErrors)
np.savetxt(home+'data/active/superRedshift.dat',outputRedshift)
np.savetxt(home+'data/active/superInfo.dat',outputInfo)
np.savetxt(home+'data/active/truths.dat',outputTruths)
np.savetxt(home+'data/active/masses.dat',outputMasses)