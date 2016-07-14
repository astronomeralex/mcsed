# -*- coding: utf-8 -*-
import numpy as np
import math_functions as mathF
import plotter_functions as pf
import os
import shutil
import variable_house as vHouse

# =============================================================================
folderOfInterest = '/users/hunterbrooks/astronomy/results/1.4c_new_no_neb_emis/'
maxIDX = 350
predictions = np.zeros([maxIDX,5])
errors = np.zeros([maxIDX,5,2])
chisqtrims = np.zeros([maxIDX])
info = np.loadtxt(folderOfInterest+'superInfo.dat')


# =============================================================================
filterLambdas = np.loadtxt(folderOfInterest+'filter_lambdas.dat')
filterNus = mathF.c_ang / filterLambdas


# =============================================================================        
for dataIDX in range(0,maxIDX):
#for dataIDX in range(0,0):
    
    # =============================================================================
    # try to load this index
    try:
        npzdata = np.load(folderOfInterest+'output'+str(dataIDX)+'.npz')
    except IOError:
        print 'skipping ', dataIDX
        continue
    
    # =============================================================================
    # open up data 
    time_take = npzdata['time_taken']
    
    dataAvgFlux = npzdata['dataAvgFlux']
    dataAvgFluxErr = npzdata['dataAvgFluxErr']
    z = npzdata['z']
    
    data_nu_obs = npzdata['data_nu_obs']
    D_nu = npzdata['D_nu']
    R_nu = npzdata['R_nu']
    
    nDim = npzdata['nDim']
    nDram = npzdata['nDram']
    nWalkers = npzdata['nWalkers']
    nSteps = npzdata['nSteps']
    nBurnInSteps = npzdata['nBurnInSteps']
    
    paramNames = npzdata['paramNames']
    #dramNames = npzdata['dramNames']
    paramRanges = npzdata['paramRanges']
    
    walker_SED_type = npzdata['walker_SED_type']
    walker_SED_notes = npzdata['walker_SED_notes']
    
    #sampler_type = npzdata['sampler_type']
    chain = npzdata['chain']
    lnprobability = npzdata['lnprobability']
    auto_corr = npzdata['auto_corr']
    acc_frac = npzdata['acc_frac']
    flatchain = npzdata['flatchain']
    flatlnprobability = npzdata['flatlnprobability']
    
    npzdata.close()
    field = str(int(info[dataIDX, 3]))
    objid= str(int(info[dataIDX,4]))    
    
    # ======================================================================= #
    # get path to object
    obj_folder = folderOfInterest+str(dataIDX)+'/'
    # make a folder for object
    if os.path.exists(obj_folder) == False:
        os.mkdir(obj_folder)
    
    '''
    # ======================================================================= #
    # TEMPORARY FIX to s/n issue
    okidx = np.where(dataAvgFluxErr != 0)[0]
    dataAvgFlux_temp = np.zeros(np.shape(dataAvgFlux))
    for i in range(0, len(okidx)):
        idx = okidx[i]
        dataAvgFlux_temp[idx] = dataAvgFlux[idx]
    dataAvgFlux = dataAvgFlux_temp
    '''
    """
    # ======================================================================= #
    # see if this object has any observational bands below 1250 rest
    min_index = np.where(np.nonzero(dataAvgFlux) == np.min(np.nonzero(dataAvgFlux)))[0]
    min_wavelength = filterLambdas[min_index]    
    print '    min idx: ', min_index
    print '    min wave:', min_wavelength

    if min_wavelength/(1.0+z) < 1250.0:
        print 'OBJECT HAS OBS BELOW 1250 REST'
        continue
    """
    # ======================================================================= #
    # make M-star into log
    flatchain[:,-1] = np.log10(flatchain[:,-1])
    
    # get preidcitons and errors
    flatchain = pf.chiChop(flatchain,flatlnprobability, 3)
    chisqtrims[dataIDX] = np.shape(flatchain)[0]
    predictions[dataIDX,:] = np.median(flatchain, axis=0)
    
    #np.savez('/Users/hunterbrooks/Desktop/trimmed_chains/'+str(field)+'-'+str(objid)+'.npz',flatchain)

    # get errors, 84th and 16th percentiles for upper and lower err
    errors[dataIDX,:,0] = np.percentile(flatchain, q=16, axis=0)
    errors[dataIDX,:,1] = np.percentile(flatchain, q=84, axis=0)
    
    # ======================================================================= #
    """
    pf.plotParamPath(chain, dataIDX,
                     save_opt=obj_folder)
    
    pf.plotLnprobPath(lnprobability, dataIDX, 
                      save_opt=obj_folder)
    
    
    pf.plotTriangle(flatchain, field, objid, dataIDX, 
                    save_opt=obj_folder)
    
    sp = vHouse.defSP()
    a1 =  np.zeros([1,nDim])
    a1[0,:] = predictions[dataIDX,:-1]
    
    pf.plotSEDandResids(a1,
                        filterLambdas, dataAvgFlux, dataAvgFluxErr, data_nu_obs,
                        z, dataIDX, field, objid,
                        save_opt=folderOfInterest+'graphs/'+str(field)+'-'+str(objid)+'-sed.png')

    # ======================================================================= #
    # save observational data
    
    shutil.copy(obj_folder+'triangle.png',
                folderOfInterest+'graphs/'+field+'-'+objid+'-triangle.png')
    """        
    #np.savetxt(obj_folder+'/dataAvgFlux.dat', dataAvgFlux)
    
    '''
    D_nu_m = np.zeros(np.shape(R_nu))
    NU_obs_m = np.zeros(np.shape(R_nu))
    for i in range(0, len(filterLambdas)):
        D_nu_m[i,:-1] = D_nu
        NU_obs_m[i,:] = data_nu_obs
        

    #pf.testInterp( -0.80, -.5, 500, z,
    pf.testInterp( -0.70, -.65, 3, z,
                   
                    len(filterLambdas),
                    NU_obs_m,
                    D_nu_m, 
                    R_nu,
                   
                    dataAvgFlux,
                    dataAvgFluxErr)
    '''

    print '... done ', dataIDX, ' or  ', field, objid

    
    
    
# =============================================================================
np.save(folderOfInterest+'predictions.dat', predictions)
np.save(folderOfInterest+'errors.dat', errors)
np.save(folderOfInterest+'chisqtrims.dat', chisqtrims)

#predictions = np.load(folderOfInterest+'predictions.dat.npy')   
#errors = np.load(folderOfInterest+'errors.dat.npy')
#myField = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/superInfo.dat')[:,3]
#myID = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/superInfo.dat')[:,4]

"""
flags = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/flags.dat')
pf.compareMassRatioVsRedshift(flags[:,5], flags[:,6],
                    predictions[:,4], errors[:,4], 
                    flags[:,0], flags[:,1], 
                    info[:,3], info[:,4],
                    np.loadtxt('/users/hunterbrooks/astronomy/dust/data/active/superRedshift.dat'))
"""

# PRETY MUCH BROKEN
#pf.plotSelfAny2(predictions, errors, 1, 4, 0, KC13=False)

#pf.histogram(predictions, 2, 'Eb', errors)

# MY LIST: 1,2,0,3
# HIS LIST: 4 (My M_star), 12 (12+log(O/H)), 15 (ellipticity+err), 

#alexs_data = np.loadtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat',usecols=12, skiprows=2)
alexs_data = np.genfromtxt('/users/hunterbrooks/astronomy/dust/data/z2/data_house/alex.dat')
aField = alexs_data[:,7]
aID = alexs_data[:,8]

p1_idx = 2

alex_sorted = np.zeros(len(predictions[:,0]))
alex_sorted_e = np.zeros(len(predictions[:,0]))

def matchbyfield(infield, inid):
    index = -1
    for i in range(0, np.shape(aID)[0]):
        if infield == aField[i]:
            if inid == aID[i]:
                index = i
        
    if index == -1:
        return False
    
    return index
    
for ii in range(0, len(myField)):
    his_idx = matchbyfield(myField[ii], myID[ii])
    if his_idx != False:
        alex_sorted[ii] = alexs_data[his_idx, p1_idx]

pf.plotGlobalAny2(     predictions,#alex_sorted, 
                       errors,#alex_sorted_e, 
                       2,#0, #p1_idx, 
                       
                       predictions, #p2_data,
                       errors, #p2_error,
                       2, #p2_idx,
                       
                       r'Half Light Radius vs. $E_b$',
                       'Half Light Radius from alex.dat [Kpc]', #xlabel,
                       r'$E_b$', #ylabel
                       
                       pColor_idx=None,
                       p1_mine=True,
                       p2_mine=True)

