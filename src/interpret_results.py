# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:14:49 2015

@author: hunterbrooks

This code is designed to interpret results from the random walk chains.  It
makes a catalog from the directory with variable name "pathToResults".  This
catalog has columns: 

Field,ID,Redshift,
dust2(16),dust2(50),dust2(84),
delta(16),delta(50),delta(84),
E_b(16),E_b(50),E_b(84),
log10age(16),log10age(50),log10age(84)		
log10mass(16),log10mass(50),log10mass(84)

for a total of 18 columns.  The numbers in parenthesis indicate the percentile
value of the posterior distribution; these numbers are the widths of the
low/high errors bars.

This catalog has however many rows there are .npz files in "pathToResults".

"""

# =========================================================================== #
def plotTriangle(new_flatchain, field, objid, dataIDX, save_opt=None, copy_opt=None):
    import triangle
    import matplotlib.pylab as plt
    vNames = ['tau (dust2)','delta (dust_index)', 'Eb (uvb)', 'log10(age)', 
              'log10(M_star)']
    triangle_figure = triangle.corner(new_flatchain, labels=vNames,
                                      quantiles=[0.16, 0.5, 0.84],
                                      show_titles=True, title_args={"fontsize": 12},
                                      plot_contours=True,
                                      plot_datapoints=False,
                                      verbose=False
                                      )
                                    
    
    triangle_figure.gca().annotate(field+'-'+objid+', hunters index: '+str(dataIDX)+', after trimming', 
                        xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")
    
                      
    if save_opt != None:
        triangle_figure.savefig(save_opt+'triangle.png')
        plt.close('all')
        if copy_opt != None:
            shutil.copyfile(save_opt+'triangle.png', copy_opt+'triangle-'+str(field)+'-'+str(objid)+'.png')
            
    else:
        triangle_figure.show()
        
# =========================================================================== #
def chiChop(fchain, flnprob, nSigmas):
    ''' Takes flatchain and flatlnprob and chops the lnprob() values
    at the lnprob equivilent of whichSigma where whichSigma is the
    gaussian error.  This is outlier chopping in lnprob space'''
    if nSigmas==1:
        c= 0.5
    if nSigmas==2:
        c= 2.0
    if nSigmas==3:
        c= 4.5
    if nSigmas==4:
        c= 8.0
    chisqmax = np.max(flnprob)
    # get ok indexes
    l = np.where(flnprob >= chisqmax-c)[0]
    # fill new flatchain
    new_fchain = np.zeros( [len(l), np.shape(fchain)[1]] )
    
    for i in range(0, len(l)):
        idx = l[i]
        new_fchain[i,:] = fchain[idx,:]
    
    return new_fchain
    
# =========================================================================== #
def plotLnprobPath(lnprobability, field, objid, OBJi, save_opt=None, copy_opt=None):
    nSteps = np.shape(lnprobability)[1]
    nWalkers = np.shape(lnprobability)[0]
    import matplotlib.pylab as plt
    stepidx = np.arange(1,nSteps+1)
    #plots lnprob for each walker over step idx
    for walkeri in range(0,nWalkers):
        plt.plot(stepidx, lnprobability[walkeri,:], alpha=0.1, color='k')
    plt.xlabel('STEP IDX')
    plt.ylabel('lnprob()')
    plt.title('lnprob of walkers after trimming')
    if save_opt != None:
        plt.savefig(save_opt+'lnprobpath.png')
        plt.close('all')
        if copy_opt != None:
            shutil.copyfile(save_opt+'lnprobpath.png', copy_opt+'lnprobpath-'+str(field)+'-'+str(objid)+'.png')
    else:
        plt.show()
      
    
# =========================================================================== #
def plotParamPath(chain, OBJi,
                  save_opt=None):
    import matplotlib.pylab as plt
    import variable_house as vHouse
    nWalkers = np.shape(chain)[0]
    nSteps = np.shape(chain)[1]
    nDim = np.shape(chain)[2]
    stepidx = np.arange(1,nSteps+1)
    
    for param_i in range(0,nDim-vHouse.nDram):
        for walker_i in range(0,nWalkers):
                plt.plot(stepidx, chain[walker_i,:,param_i],
                         color='k',
                         alpha= 0.1)
        
        plt.xlabel('STEP INDEX')
        plt.ylabel(vHouse.nameL[param_i])
        plt.title('param value of walkers after trimming')
        #show entire range of possible parameter values
        plt.ylim([np.min(vHouse.paramRanges[param_i,:]),
                  np.max(vHouse.paramRanges[param_i,:])]) 
        if save_opt != None:
            plt.savefig(save_opt+'param'+str(param_i)+'path.png')
            plt.close('all')
        else:
            plt.show()   
    
# =========================================================================== #
    
    
    
# necesary during loop
import re
import shutil


# get however many output*.npz files are in folder
pathToResults = '/users/hunterbrooks/astronomy/summer_run2/'
import glob
pathToEach = glob.glob(pathToResults+'output*') #path to each .npz file
nResults = len(pathToEach)
print 'There are '+str(nResults)+' output*.npz files in this folder.'


# make a folder for all triangle plots and all other plots
import os
if os.path.exists(pathToResults+'all_triangles/') == False:
    os.mkdir(pathToResults+'all_triangles/')
if os.path.exists(pathToResults+'all_lnprob/') == False:
    os.mkdir(pathToResults+'all_lnprob/')


# load info used later
import numpy as np
info = np.loadtxt(pathToResults+'superInfo.dat')


# create blank catalog
#Field,ID,Redshift,
#dust2(16),dust2(50),dust2(84),
#delta(16),delta(50),delta(84),
#E_b(16),E_b(50),E_b(84),
#log10age(16),log10age(50),log10age(84)		
#log10mass(16),log10mass(50),log10mass(84)
nCatalogCols = 18
catalog = np.zeros([nResults,nCatalogCols])



# loop through each object
for objIDX in range(0, nResults):
    print
    print objIDX
    
    # record path to object as objName
    objName = pathToEach[objIDX]
    print objName
    
    
    # find out what row in my catalogs this object is
    myIndex = int(re.search(r'\d+', objName[-8:]).group())
    
    # myIndex is the index tthat the info for this object is in for superInfo
    print myIndex
    
    # open up data for object
    npzdata = np.load(objName)
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
    dramNames = npzdata['dramNames']
    paramRanges = npzdata['paramRanges']
    walker_SED_type = npzdata['walker_SED_type']
    walker_SED_notes = npzdata['walker_SED_notes']
    sampler_type = npzdata['sampler_type']
    chain = npzdata['chain']
    lnprobability = npzdata['lnprobability']
    auto_corr = npzdata['auto_corr']
    acc_frac = npzdata['acc_frac']
    flatchain = npzdata['flatchain']
    flatlnprobability = npzdata['flatlnprobability']
    npzdata.close()
    field = str(int(info[  myIndex  , 3]))
    objid= str(int(info[  myIndex  , 4]))
    
    # make a folder for this specific object
    folderForThisObj = pathToResults+field+'-'+objid+'/'
    if os.path.exists(folderForThisObj) == False:
        os.mkdir(folderForThisObj)
        
    # make M-star into log
    flatchain[:,-1] = np.log10(flatchain[:,-1])
    
    # trim flatchain in gaussian space to Nsigmas sigmas
    Nsigmas = 3
    flatchain = chiChop(flatchain,flatlnprobability, Nsigmas)
    
    # take 50th percentile of flatchain to be param estimate
    paramEstimates = np.median(flatchain, axis=0)
 
    # get errors, high and low errors 
    paramLowErrors = np.abs(np.percentile(flatchain, q=16, axis=0) - np.percentile(flatchain, q=50, axis=0))
    paramHighErrors= np.abs(np.percentile(flatchain, q=84, axis=0) - np.percentile(flatchain, q=50, axis=0))
    
    # save information in catalog
    catalog[objIDX,0] = int(field)
    catalog[objIDX,1] = int(objid)
    catalog[objIDX,2] = z
    catalog[objIDX,3] = paramLowErrors[0]
    catalog[objIDX,4] = paramEstimates[0]
    catalog[objIDX,5] = paramHighErrors[0]
    catalog[objIDX,6] = paramLowErrors[1]
    catalog[objIDX,7] = paramEstimates[1]
    catalog[objIDX,8] = paramHighErrors[1]
    catalog[objIDX,9] = paramLowErrors[2]
    catalog[objIDX,10] = paramEstimates[2]
    catalog[objIDX,11] = paramHighErrors[2]    
    catalog[objIDX,12] = paramLowErrors[3]
    catalog[objIDX,13] = paramEstimates[3]
    catalog[objIDX,14] = paramHighErrors[3]
    catalog[objIDX,15] = paramLowErrors[4]
    catalog[objIDX,16] = paramEstimates[4]
    catalog[objIDX,17] = paramHighErrors[4]
    
    # make triangle plot and save it
    plotTriangle(flatchain, field, objid, myIndex, 
                 save_opt=folderForThisObj, 
                 copy_opt=pathToResults+'all_triangles/')
                 
    # plot ln path
    plotLnprobPath(lnprobability, field, objid, myIndex, 
                      save_opt=folderForThisObj, 
                      copy_opt=pathToResults+'all_lnprob/')
                      
    # plot param path
    plotParamPath(chain, myIndex,
                  save_opt=folderForThisObj)
                  
    print 'Done: ', objName
                  
                  

# save catalog
np.savetxt(pathToResults+'catalog.txt',catalog)  

'''bc large paramete space, get multimodal posteriors''' 