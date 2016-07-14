# -*- coding: utf-8 -*-
"""
This script contains all the functions assoicated with the emcee random walk,
along with the process to do the walk itself.  The important function here, 
especially if the user is using DRams, is lnprob().
"""
# =========================================================================== #
import numpy as np
import variable_house as vHouse
import factory
import math_functions as mathF
import mcsed_utils as myTools
import scipy.interpolate as spterp

sp_temp = vHouse.defSP()
conroys_ages = sp_temp.log_age - 9.0 #Gyr



# =========================================================================== #
# check to make sure parameters are in bounds
def inBoundsCheck(x):
    for pi in range(0,vHouse.nDim):
        if x[pi] > np.max(vHouse.higherL[pi]):
            return False
        if x[pi] < np.min(vHouse.lowerL[pi]):
            return False       
    # if escaped loop without returning False
    return True

           
# =========================================================================== #
# The probability distribution that you would like to sample.
def lnprob( x,
            flux_interp_f,
            lambda_e,
            nFilters,
            z,
            NU_obs,
            D_nu, 
            R_nu,
            dataAvgFlux,
            dataAvgFluxErr,
            A_interp_f=None,
            stellarmass_interp_f=None):
    # check to make sure parameters are in bounds
    if inBoundsCheck(x) == False:
        if vHouse.sampler_type == 'Ensemble':
            return -1.0*np.inf, np.array([None]*vHouse.nDram)
        if vHouse.sampler_type == 'PT':
            return -1.0*np.inf
    '''
    if vHouse.walker_SED_type == 'Direct':
        # stellar population is passed into lnprob as flux_interp_f 
        sp, tage = vHouse.updateFSPS(flux_interp_f, vHouse.nameL, x)
        model_l_lambda = flux_interp_f.get_spectrum(tage=tage)[1]
        model_s_nu = mathF.getFlux(model_l_lambda, lambda_e, z)[0]
        model_S_nu = factory.S_nu(nFilters, model_s_nu)
        modelAvgFlux = mathF.getAvgFlux(     nFilters,
                                             model_S_nu,
                                             NU_obs,
                                             D_nu,
                                             R_nu,
                                             True)
        A = mathF.getA(modelAvgFlux, dataAvgFlux, dataAvgFluxErr)
        chi_sq = myTools.fracsum((dataAvgFlux-A*modelAvgFlux)**2,(dataAvgFluxErr)**2)
        if vHouse.sampler_type == 'Ensemble':
            return (-1.0/2.0)*chi_sq, A*sp.stellar_mass
        if vHouse.sampler_type == 'PT':
            assert True==False, "PT + Direct is currenltly unsupported b/c blobs"
    '''
    # this is a work around because sfh == 1 doesnt interpolate nicely
    if vHouse.walker_SED_type == 'Direct':
    
        # stellar population is passed into lnprob as flux_interp_f
        logtage = x[3]
        clsti =  np.argmin(np.abs(conroys_ages-logtage))
        # test to make sure not exactly halfway between !!
        
        if conroys_ages[clsti] == logtage:
            sp, tage = vHouse.updateFSPS(flux_interp_f, vHouse.nameL, x)
            model_l_lambda = flux_interp_f.get_spectrum(tage=tage)[1]
            model_s_nu = mathF.getFlux(model_l_lambda, lambda_e, z)[0]
              
        if conroys_ages[clsti] > logtage or conroys_ages[clsti] < logtage:
            if conroys_ages[clsti] > logtage:
                # indexes
                lowi = clsti - 1
                hii = clsti
            if conroys_ages[clsti] < logtage:
                # indexes
                lowi = clsti
                hii = clsti + 1
                
            # param vals for nearest neighboring two points
            lowparams = np.array(x,copy=True)
            lowparams[3] = conroys_ages[lowi]
            highparams = np.array(x,copy=True)
            highparams[3] = conroys_ages[hii]
            
            
            # get high and lower SEDs on conroys points
            sp, tage = vHouse.updateFSPS(flux_interp_f, vHouse.nameL, lowparams)
            low_model_l_lambda = flux_interp_f.get_spectrum(tage=tage)[1]
            sp, tage = vHouse.updateFSPS(flux_interp_f, vHouse.nameL, highparams)
            high_model_l_lambda = flux_interp_f.get_spectrum(tage=tage)[1]
            
            #create x and y inputs of the interp function
            interpy = np.zeros([2,np.shape(low_model_l_lambda)[0]])
            interpy[0,:] = low_model_l_lambda
            interpy[1,:] = high_model_l_lambda
            interpx = np.array([conroys_ages[lowi],conroys_ages[hii]])
            logtage_interp_f = spterp.interp1d(interpx, interpy, kind='linear', axis=0)
            
            # use the interp function
            model_l_lambda = logtage_interp_f(logtage)
            model_s_nu = mathF.getFlux(model_l_lambda, lambda_e, z)[0]
            '''
            print '------------------ ------------------- -----------------'
            print lowparams
            print x
            print highparams
            print
            print logtage_interp_f(conroys_ages[lowi])[:5]
            print model_l_lambda[:5]
            print logtage_interp_f(conroys_ages[hii])[:5]
            raw_input()
            '''
            

        model_S_nu = factory.S_nu(nFilters, model_s_nu)
        modelAvgFlux = mathF.getAvgFlux(     nFilters,
                                             model_S_nu,
                                             NU_obs,
                                             D_nu,
                                             R_nu,
                                             True)
                                             
        A = mathF.getA(modelAvgFlux, dataAvgFlux, dataAvgFluxErr)
        chi_sq = myTools.fracsum((dataAvgFlux-A*modelAvgFlux)**2,(dataAvgFluxErr)**2)
        if vHouse.sampler_type == 'Ensemble':
            # get new SP with exact params
            sp, tage = vHouse.updateFSPS(flux_interp_f, vHouse.nameL, x)
            return (-1.0/2.0)*chi_sq, A*sp.stellar_mass
        if vHouse.sampler_type == 'PT':
            assert True==False, "... just dont use this."    
    
    if vHouse.walker_SED_type == 'Interp':
        model_s_nu = flux_interp_f(x)
        # need to reshape from (1,1963) to (1963,) due to interpolation function
        model_s_nu = np.reshape(model_s_nu, np.shape(model_s_nu)[1])
        model_S_nu = factory.S_nu(nFilters, model_s_nu)
        modelAvgFlux = mathF.getAvgFlux(     nFilters,
                                             model_S_nu,
                                             NU_obs,
                                             D_nu,
                                             R_nu,
                                             True)
        chi_sq = myTools.fracsum((dataAvgFlux-modelAvgFlux)**2,(dataAvgFluxErr)**2)
        if vHouse.sampler_type == 'Ensemble':
            return (-1.0/2.0)*chi_sq, A_interp_f(x)*stellarmass_interp_f(x)
        if vHouse.sampler_type == 'PT':
            return (-1.0/2.0)*chi_sq #no blobs


# =========================================================================== #    
def lnprior(x):
    ''' dummy log prior function to keep PT sampler happy'''
    return 0.0


# =========================================================================== #    
def conductMCMC(    flux_interp_f,
                    A_interp_f,
                    stellarmass_interp_f,
                    lambda_e,
                    nFilters,
                    z,
                    NU_obs,
                    D_nu, 
                    R_nu,
                    dataAvgFlux,
                    dataAvgFluxErr): 
    # Setup Emcee variables
    import emcee
    nDim = vHouse.nDim
    nDram = vHouse.nDram
    nWalkers = vHouse.nWalkers
    nSteps = vHouse.nSteps
    nBurnSteps = vHouse.nBurnInSteps
    global lnprob
    if vHouse.sampler_type == 'Ensemble':
        p0 = np.zeros([nWalkers, nDim])
        for pi in range(0,nDim):
            p0[:, pi] = np.random.uniform(vHouse.lowerL[pi],
                                          vHouse.higherL[pi],
                                          [nWalkers])
        # clear random state to begin
        # create sampler
        sampler = emcee.EnsembleSampler(nWalkers, nDim, lnprob, 
                                        args=[  flux_interp_f,
                                                lambda_e, 
                                               
                                                nFilters,
                                                z,
                                                NU_obs,
                                                D_nu, 
                                                R_nu,
                                               
                                                dataAvgFlux,
                                                dataAvgFluxErr, 
                                                
                                                A_interp_f,
                                                stellarmass_interp_f])
        if vHouse.run_burnin_opt == True:        
            # run a burn in
            pos, prob, state, blobs = sampler.run_mcmc(p0, nBurnSteps)
            # Reset the chain to remove the burn-in samples.
            sampler.reset()
            maxlnprob_idx = np.argmax(prob) # get index of maxlnprob
            maxlnprob_pos = pos[maxlnprob_idx]
            # Starting a ball centered at from pos, run emcee
            p0 = np.zeros([nWalkers, nDim])
            for pi in range(0,nDim):
                p0[:,pi] = maxlnprob_pos[pi] + vHouse.ball_radius[pi]*np.random.randn(nWalkers)
            sampler.run_mcmc(pos, nSteps, rstate0=state)
        if vHouse.run_burnin_opt == False:
            sampler.run_mcmc(p0, nSteps)
            
        # reshape blobs to put into sampler.chain shape
        blobs_chain = np.swapaxes(sampler.blobs, 0,1)
        new_chain = np.zeros([nWalkers, nSteps, nDim+nDram])
        new_flatchain = np.zeros([nWalkers*nSteps, nDim+vHouse.nDram])
        new_chain[:,:,:nDim] = sampler.chain
        new_chain[:,:,-nDram] = blobs_chain
        new_flatchain[:,:nDim] = sampler.flatchain
        
        for nw in range(0, nWalkers):
            for ns in range(0, nSteps):
                flatidx = np.where(sampler.chain[nw,ns,:]==sampler.flatchain)[0]
                new_flatchain[flatidx,:nDim]=sampler.chain[nw,ns,:]
                new_flatchain[flatidx,-nDram]=blobs_chain[nw,ns]
            
        return (new_chain,
                sampler.lnprobability,
                sampler.get_autocorr_time(), 
                sampler.acceptance_fraction,
                new_flatchain,
                sampler.flatlnprobability)
    if vHouse.sampler_type == 'PT':
        nTemps = vHouse.nTemps
        p0 = np.zeros([nWalkers, vHouse.nDim])
        for pi in range(0,vHouse.nDim):
            p0[:, pi] = np.random.uniform(vHouse.lowerL[pi],
                                          vHouse.higherL[pi],
                                          [nWalkers])
        # clear random state to begin                
        # create sampler
        sampler = emcee.PTSampler(nTemps, nWalkers, nDim, 
                                  lnprob, lnprior,
                                  loglargs=[    flux_interp_f,
                                                lambda_e, 
                                               
                                                nFilters,
                                                z,
                                                NU_obs,
                                                D_nu, 
                                                R_nu,
                                               
                                                dataAvgFlux,
                                                dataAvgFluxErr])
        # pick random initial
        p0 = np.zeros([nTemps, nWalkers, nDim])
        for i in range(0, nDim):
            p0[:,:,i] = np.random.uniform(  low=vHouse.lowerL[i], 
                                            high=vHouse.higherL[i],
                                            size=(nTemps, nWalkers))
        if vHouse.run_burnin_opt == True:
            # burn in
            for p, lnprob, lnlike in sampler.sample(p0, iterations=nBurnSteps):
                pass
            # reset sampler after run
            sampler.reset()
            # do actual run
            for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                                       lnlike0=lnlike,
                                                       iterations=nSteps):
                pass
        if vHouse.run_burnin_opt == False:
            for p, lnprob, lnlike in sampler.sample(p0, iterations=nSteps):
                pass
        # recover blobs
        blobsflatchain = factory.posToStellarMass(sampler.flatchain, #[nTemps,nWalkers,nSteps,nDim] 
                                         A_interp_f, 
                                         stellarmass_interp_f)
        # add blobsflatchain to flatchain  
        new_flatchain = np.zeros([nTemps, nWalkers*nSteps, nDim+1])
        new_flatchain[:,:,:-1] = sampler.flatchain
        new_flatchain[:,:,-1] = np.log10(blobsflatchain)
        # make a flat lnprob
        flatlnprobability = 0
        return      (sampler.chain,
                    sampler.lnprobability,
                    sampler.get_autocorr_time(), 
                    sampler.acceptance_fraction,
                    new_flatchain,
                    flatlnprobability)