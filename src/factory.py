# -*- coding: utf-8 -*-
"""
This script contains all the functions that make other things from inputs. Most
are mathematical in nature, but all convert inputs to outputs in some way. I 
call these functions "make functions" hence the name of factory.py.

See the comments below each individual function.
"""

# =========================================================================== #
import numpy as np

# SOME CONSTANTS
c_ang = 2.998e18 #ang/s

# =========================================================================== #
def D_nu(nFilters, nu_obs):
    ''' 
    GENERATES D_nu, OR THE SPACE BETWEEN NU VALUES USED IN THE FINDING OF
    AVERAGE FLUX IN A BAND.  THIS IS NECESARY BECAUSE INSTEAD OF DOING AN
    INTEGRAL, WE DO A SUM.
    
    NOTE: D_nu IS UNTRIMMED and the same length as nu_obs.
    '''

    D_nu = np.zeros([nFilters, len(nu_obs)])
    for nu_i in range(0, len(nu_obs)-1):
        D_nu[:,nu_i] = np.abs(nu_obs[nu_i+1] - nu_obs[nu_i])
        
    return D_nu

# =============================================================================
def NU_obs( nFilters, 
              nu_obs):
    ''' 
    This function creates a 2d array where each row is simply a copy of
    nu_obs.  This is done because when calculating average fluxes, matrix math
    is faster than doing each filter individually.  The output has shape
    [nFilters, len(nu_obs)].
    
    NOTE: IS UNTRIMMED.
    '''
    
    NU_obs = np.zeros([nFilters, len(nu_obs)])
    
    for filt_i in range(0, nFilters):
        NU_obs[filt_i,:] = nu_obs
        
    return NU_obs
    
# =============================================================================
def R_nu( home,
              nu_obs):
    '''CREATES A TRANSMISSION ARRAY OF SHAPE [nFilters, len(lamda_e)] WHICH IS 
    JUST nFilters ROWS AND THE LENGTH OF THE FSPS WAVELENGTH OUTPUT COLS.  THE
    REASON THIS ARRAY HAS nFilters ROWS IS THAT MATRIX MATH IS FASTER WHEN
    COMPUTING THE AVERAGE FLUX IN A BAND, SO WE DO THEM ALL AT ONCE.
    
    INPUTS: home = STRING POINTING TO mcsed/ FOLDER
            nu_obs = OBSERVED FRAME NU'S FROM FSPS, SPECIFIC TO AN OBJECT
    
    RETURNS: R_nu = SHAPE [nFilters, len(lamda_e)] WHERE EACH ROW IS A
                    TRANSMISSION ARRAY FOR THE CORRISPONDING ROW FILTER IN
                    columns.dat. EACH ROW IS NORMALIZED TO 1 AND EACH 
                    TRANSMISSION POINT CORRIPSONDS TO A POINT IN nu_obs.
                    R_nu IS UNTRIMMED.
            nFilters = SIMPLY THE NUMBER OF FILTERS IN columns.dat'''
    
    import scipy.interpolate as spterp
    
    # create an R_nu array
    col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')
    nFilters = len(col_names)
    R_nu = np.zeros([nFilters, len(nu_obs)])
    
    for filt_i in range(0,nFilters):
        data = np.loadtxt(home+'data/filters/filter_house/'+col_names[filt_i]+'.res')
        lambdas = data[:,0] #units= Ang
        nus_o = c_ang / lambdas #units=[per sec]
        res_o = data[:,1]
        
        filter_interp_function = spterp.interp1d( nus_o , res_o , 
                                                          bounds_error=False , 
                                                          fill_value=0 )
        res = filter_interp_function(nu_obs)
        res = res / np.max(res)
        
        R_nu[filt_i,:] = res 
        # R_nu has not been trimmed yet!
        
    return R_nu, nFilters
    
# =============================================================================    
def S_nu(nFilters, s_nu):
    
    import numpy as np
    
    S_nu = np.zeros([nFilters, len(s_nu)])

    for filt_i in range(0,nFilters):
        S_nu[filt_i,:] = s_nu
    
    return S_nu
    
# =============================================================================
def fakeSED(   trueNames, trueVals,
                noise_std_dev,
                z):
    ''' RETURNS A OBSERVERD FLUX PSUEDO-SED BASED ON THE INPUT PARAMETERS. 
    THIS SED HAS A RANDOM GAUSSIAN NOISE ASSOCIATED WITH IT.
    
    INPUTS: fake_params = FAKE PARAMETERS
            noise_std_dev = STANDARD DEVIATION OF GAUSSIAN NOISE
            z = REDSHIFT
            
    RETURNS: s_nu = 1D ARRAY OF FLUX DENSITY 
             nu_obs = OBSERVED FRAME FREQUENCY FOR EACH POINT IN s_nu'''
             
    import variable_house as vHouse
    import math_functions as mathf
    
    sp = vHouse.defSP()
    sp, tage = vHouse.updateFSPS(sp, trueNames, trueVals)
    waves, spec = sp.get_spectrum(tage=tage)
    
    # add random gaussian noise
    noise = np.random.normal(1, noise_std_dev, len(waves))
    lum_e = spec * noise
    
    # convert to fluxes, at our redshift
    fluxOUTPUTS = mathf.getFlux(lum_e, waves, z)
    
    return sp, fluxOUTPUTS[0], fluxOUTPUTS[1]
    
# =============================================================================
def A_grid(fluxes,
                 nFilters, NU_obs, D_nu, R_nu, 
                 dataAvgFlux, dataAvgFluxErr):
    ''' For each flux sed in fluxes, find the standardization value and return
    one '''
    
    # A grid has shape as fluxes except for last dimension
    A_grid = np.zeros(np.shape(fluxes)[:-1])
    # A_compatable has exact same shape as it
    A_compatible = np.zeros(np.shape(fluxes))
    
    #See Cartesian in create_pSpace.py.. turn N dimension fluxes into 2d
    import mcsed_utils as myTools
    import variable_house as vHouse
    import math_functions as mathF
    comboParamValues = myTools.cartesian(vHouse.paramRanges)
    comboRows = np.shape(comboParamValues)[0]
    
    # The meat of the function
    for row_idx in range(0,comboRows):
        
        # Trick to get parameter idx from actual parameter value, cost of using
        # cartesian.
        pVals = comboParamValues[row_idx,:]
        pIndexs = np.zeros(vHouse.nDim)
        
        for i in range(0,vHouse.nDim):
            # get corrisponding index of parameter value
            pIndexs[i] = int(np.where(vHouse.paramRanges[i,:] == pVals[i])[0])
            
            
        # because we use an array to index an array, we need to make it a tuple
        pTuple = tuple(pIndexs)
        
        # turn this flux sed into average fluxes in desired bands
        modelAvgFlux = mathF.getAvgFlux(     nFilters,
                                             fluxes[pTuple], #needed to tuple b/c here
                                             NU_obs,
                                             D_nu,
                                             R_nu,
                                             False)


        # actually get A for this flux sed
        A = mathF.getA(modelAvgFlux, dataAvgFlux, dataAvgFluxErr)
        
        # finally, save A for this flux sed
        A_grid[pTuple] = A
        A_compatible[pTuple] = A
    
    
    return A_grid, A_compatible
    
# =============================================================================
def posToStellarMass(fchain, Afunct, Mfunct):
    ''' given a PT sampler.flatchain, A_interp_f, stellarmass_interp_f, this code
    retreaves every parameter position in the chain and gets a stellar mass
    for it.  This is a work around for PT sampler not being able to handle
    blobs '''
    masses = np.zeros(np.shape(fchain)[:-1])  
    for t in range(0,np.shape(fchain)[0]):
        for ws in range(0,np.shape(fchain)[1]):  
            x = fchain[t,ws,:]
            masses = Afunct(x)*Mfunct(x)
                
    return masses