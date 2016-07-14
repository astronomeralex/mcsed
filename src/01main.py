# -*- coding: utf-8 -*-
"""
This script runs MCSED on objects and is organized in the same BLOCK method
that variable_house.py is.  A user should have read through variable_house.py
before interacting with main.py.  Heres an overview:

Block 1 - imports, ignore
Block 2 - for running multiple scripts. *INTERACTABLE*
Block 3 - set up how walkers obtain SEDs during random walk
Block 4 - defination of the stellar population/PRErams, user interactible 
Block 5 - begin the for loop to run MCSED on each object
Block 6 - handle fluxes and errors, including error floors
Block 7 - begin actual random walk. pass lnprob() additional inputs here.
Block 8 - save results of walk in .npz
"""

VERSION = '1.4c' #see docs/change_log.txt

# = BLOCK 1 ========================================================= BLOCK 1 =
# Be polite and print a title.
print "\n ** MCSED ** \n"
print "written by Hunter Brooks, hzb5080@psu.edu"
print "Penn State, Dep't of Astronomy"
print "version:", VERSION, "\n"

# remove all old .pyc files to make sure no funny business is going on
print '... Cleaning old .pyc files'
import os
directory = os.listdir('.')
for filename in directory:
    if filename[-3:] == 'pyc':
        os.remove(filename)

# import everything that is needed
print '... Doing initial imports'
import time
start_time = time.time() # store start time 
import numpy as np
import emcee_functions as mcF
import factory
import math_functions as mathF
# get home
import wheres_home
home = wheres_home.getHomeLocation()
import variable_house as vHouse
import sys


# = BLOCK 2 ========================================================= BLOCK 2 =
nScriptsRunning = 1 #up to 100, numbered with two digits each time
superIndexes = np.loadtxt(home+'data/active/superIndexes.dat')

# get name of this script, and which number it is
thisScriptName = sys.argv[0][-9:] # filename should be NNmain.py, eg 04main.py
thisScriptNumb = int(thisScriptName[0:2])
print thisScriptNumb
print 

# split superIndexes into nScriptsRunning chunks
if (nScriptsRunning == 1):
    chunks = np.array_split(superIndexes, nScriptsRunning)[0]
else:
    chunks = np.array_split(superIndexes, nScriptsRunning)
print chunks    

# = BLOCK 3 ========================================================= BLOCK 3 =
print '... Opening necesary MCSED data'
if vHouse.walker_SED_type == 'Interp':
    pSpaceGridSize = vHouse.pSpaceGridSize
    # being loading interpolation data... make a time stamp
    load_data_time = time.time()
    npzopendata = np.load(home+'data/pSpaces/size'+str(pSpaceGridSize)+'_sed.npz')
    sedGrid = npzopendata["sedGrid"]
    npzopendata.close()
    npzopendata = np.load(home+'data/pSpaces/size'+str(pSpaceGridSize)+'_other.npz')
    #valsGrid=npzopendata['valsGrid']  #BROKEN, see create_pSpace.py
    lambda_e=npzopendata['waves']
    stellarMassGrid=npzopendata['stellarMassGrid']
    npzopendata.close()
    import scipy.interpolate as spterp                                           
    stellarmass_interp_f = spterp.RegularGridInterpolator(vHouse.paramRanges, 
                                                  stellarMassGrid, 
                                                  method='linear')   

    walker_SED_notes = 'pSpace grid size: '+str(pSpaceGridSize)
    print '...      Took ', str(time.time() - load_data_time)
if vHouse.walker_SED_type == 'Direct':
    sp = vHouse.defSP()
    lambda_e = sp.wavelengths
    walker_SED_notes = 'Currently, direct does not work well.. be careful'


# = BLOCK 4 ========================================================= BLOCK 4 =
print '... Opening necesary observational data'
# find nFilters and load input data from mcsed/data/active
col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')
filter_lambdas = np.loadtxt(home+'data/active/filter_lambdas.dat')
nFilters = len(col_names)
# load input data, this is the data for the objects YOU want to match
superFluxes = np.loadtxt(home+'data/active/superFluxes.dat')
superErrors = np.loadtxt(home+'data/active/superErrors.dat')
superInfo = np.loadtxt(home+'data/active/superInfo.dat')
superRedshift = np.loadtxt(home+'data/active/superRedshift.dat')
# handle an irregularity that arises when there is only 1 input
if len(np.shape(superFluxes)) == 1: #there is only 1 row
    nSuperCols = np.shape(superFluxes)[0]
    # we need to reshape superdata because it looses this property when super
    # data contains only a single row.
    superFluxes = np.reshape(superFluxes, [1,nSuperCols])
    superErrors = np.reshape(superErrors, [1,nSuperCols])
    superInfo = np.reshape(superInfo, [1,np.shape(superInfo)[1]])
    
nObjects, nSuperCols = np.shape(superFluxes)
# make sure nSuperCols == nFilters
assert nFilters == nSuperCols, "All observations must have nFilters datum, even if some =0"
# make sure no negative errors of average fluxes
assert len(np.where(superErrors < 0 )[0]) == 0, "There was a negative average flux."


# = BLOCK 5 ========================================================= BLOCK 5 =
# being for loop to run MCSED on each object 
print '... Starting runs, superIndexes: ', chunks[thisScriptNumb]

for run_i in chunks:
# for run_i in chunks[thisScriptNumb]: original    
    # this keeps things simple, since run_i is just an index
    run_i = int(run_i)
    
    # check to see if this object is in /OUTPUTS
    if os.path.isfile(home+'OUTPUTS/output'+str(run_i)+'.npz') == True:
        print
        print '... '+str(run_i)+' has already been run.'
        continue
    
    # log starting time for this obj
    run_start_time = time.time()
    dataAvgFlux = np.zeros([nFilters])
    dataAvgFluxErr = np.zeros([nFilters])
    
    # = BLOCK 6 ===================================================== BLOCK 6 =
    z = superRedshift[run_i] #redshift for this object
    assert z >= 0.0, "Just a check to make sure Redshift is positive..."
    dataAvgFlux = np.zeros([nFilters])
    dataAvgFluxErr = np.zeros([nFilters])
    
    # feed all the data into the active variable (dataAvgFlux/Err)...
    for filt_i in range(0, nFilters):
        # pull data from superFluxes/Errors into working data
        dataAvgFlux[filt_i] = superFluxes[run_i, filt_i  ]
        dataAvgFluxErr[filt_i] = superErrors[run_i, filt_i  ]
        
    # ...then run the data/errors through filter machine
    dataAvgFlux, dataAvgFluxErr = mathF.filterMachine(    dataAvgFlux, 
                                                          dataAvgFluxErr,
                                                          z, 
                                                          filter_lambdas,
                                                          run_i)
                                                          
    # make sure that dataAvgFlux is not entirely zeros!
    nNullObs = len(np.where(dataAvgFlux==0)[0])
    if nNullObs == len(dataAvgFlux):
        print
        print '... '+str(run_i)+" has only null observations (entire dataAvgFlux = 0)"
        continue
    
    # = BLOCK 7 ===================================================== BLOCK 7 =
    print
    print '... Customizing pSpace for super row index: ', run_i
    
    # this is used to move data to observed frame
    nu_obs = mathF.get_nu_obs(lambda_e, z)
    # create an R_nu array, dont need to do this each time but its fast and
    # for confusion's sake we just do it each time...
    R_nu, nFilters = factory.R_nu(home,
                                  nu_obs)
    # create D_nu array. D_nu is the space between nu's that corrispond to 
    # sp.wavelengths
    D_nu = factory.D_nu(nFilters,
                        nu_obs)
    # create NU_obs array
    NU_obs = factory.NU_obs(nFilters,
                            nu_obs)
    if vHouse.walker_SED_type == 'Interp':                        
        # create flux_interp_f from luminosity SED grid
        fluxes, nu_obs = mathF.getFlux(sedGrid, lambda_e, z)
        # find the A value for each flux sed
        A_grid, A_compatible = factory.A_grid(fluxes, 
                                              nFilters, NU_obs, D_nu, R_nu,
                                              dataAvgFlux, dataAvgFluxErr)
        # scale fluxes by each SEDs A value
        fluxes = fluxes*A_compatible
        # create flux interp function
        flux_interp_f = spterp.RegularGridInterpolator(vHouse.paramRanges, 
                                                       fluxes, 
                                                       method='linear')                                    
        # create A interp function
        A_interp_f = spterp.RegularGridInterpolator(vHouse.paramRanges,
                                                    A_grid,
                                                    method='linear')
    if vHouse.walker_SED_type == 'Direct':
        flux_interp_f = sp
        A_interp_f = None
        stellarmass_interp_f = None
    
    print '... Conducting random walk'
    # BEGIN emcee RUN
    MCMCout = mcF.conductMCMC(      flux_interp_f,
                                    A_interp_f,
                                    stellarmass_interp_f,
                                    lambda_e, 
                                   
                                    nFilters,
                                    z,
                                    NU_obs,
                                    D_nu, 
                                    R_nu,
                                   
                                    dataAvgFlux,
                                    dataAvgFluxErr)  
                                
    print '...      Took ', str(time.time() - run_start_time)
    
    # = BLOCK 8 ===================================================== BLOCK 8 =
    # save data
    print '... Saving results of row index: ', run_i
    save_time = time.time()
    
    np.savez_compressed(home+'OUTPUTS/output'+str(run_i)+'.npz', 
                        
                        run_date=(time.strftime("%D :: %H:%M:%S")),
                        time_taken = time.time() - run_start_time,
    
                        dataAvgFlux = dataAvgFlux,
                        dataAvgFluxErr = dataAvgFluxErr,
                        z=z,

                        data_nu_obs = nu_obs,
                        D_nu = D_nu[0,:-1],
                        R_nu = R_nu,
                        
                        nDim = vHouse.nDim,
                        nDram = vHouse.nDram,
                        nWalkers = vHouse.nWalkers,
                        nSteps = vHouse.nSteps,
                        nBurnInSteps = vHouse.nBurnInSteps,
                        
                        paramNames = vHouse.nameL,
                        dramNames = vHouse.dRamNames,
                        paramRanges = vHouse.paramRanges,
                        
                        walker_SED_type = vHouse.walker_SED_type,
                        walker_SED_notes = walker_SED_notes,
                        
                        sampler_type = vHouse.sampler_type,
                        chain = MCMCout[0],
                        lnprobability = MCMCout[1],
                        auto_corr = MCMCout[2],
                        acc_frac = MCMCout[3],
                        flatchain = MCMCout[4],
                        flatlnprobability = MCMCout[5],
                        VERSION = VERSION)
                        
    print '...      Took ', str(time.time() - save_time)
