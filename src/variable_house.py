# -*- coding: utf-8 -*-
"""
This script houses the bulk of the user inputs.  To keep locations easily 
refrenceable, the code is divided into several "BLOCK"s indicated by commented
lines of equal signs. Here's an overview:

Block 1 - imports, ignore
Block 2 - EMPTY
Block 3 - class defination of parameters, ignore
Block 4 - defination of the stellar population/PRErams, user interactible 
Block 5 - defination of PARrams, user interactable *IMPORTANT*
Block 6 - defination of POSTrams, user interactable
Block 7 - defination of DRams, user interactable
Block 8 - just some definations that could be used later, ignore
Block 9 - updateFSPS() routine, see code for details. Read, but ignore.
Block 10- Random walk options *IMPORTANT*

Just as a reminder: 
PREams are parameters that apply to every SED created. 
PARams are parameters that MCSED will explore.
POSTrams are parameters that need to be set as a result of a certain PARams.
DRams are derived characteristics like M_dust or M_star
"""
# = BLOCK 1 ========================================================= BLOCK 1 =
import numpy as np
import fsps
# home location of mcsed/
import wheres_home
home = wheres_home.getHomeLocation()

nActiveThreads = 1

# = BLOCK 2 ========================================================= BLOCK 2 =
# recomment me!
walker_SED_type = 'Direct' #'Interp' or 'Direct'

if walker_SED_type == 'Interp':
    # In this block, set you desired parameter grid size.  This is effecively a
    # density of SEDs to be created over the parameter ranges set in block 5. Since
    # walkers linearly interpolate during their random walk, a higher density of
    # SEDs in desireable.  The parameter space grid is regularly shaped so that
    # this density is the density for all parameters.
    pSpaceGridSize = 13
    
if walker_SED_type == 'Direct':
    pass


# = BLOCK 3 ========================================================= BLOCK 3 =
# defining the parameter class, don't touch this guy.
lowerL = []
higherL = []
nameL = []
paramRanges = []

class param:
    def __init__(self, name, lower, higher):
        self.name = name 
        self.lower = lower
        self.higher = higher
        
        if walker_SED_type == 'Interp':
            self.range= np.linspace(lower, higher, pSpaceGridSize)
            
        if walker_SED_type == 'Direct':
            # make a "range" array being only the upper and lower value becuz
            # param.range is used later on
            self.range= np.linspace(lower, higher, 2)
        
        lowerL.append(lower)
        higherL.append(higher)
        nameL.append(name)
        paramRanges.append(self.range)

            
# = BLOCK 4 ========================================================= BLOCK 4 =
## Define PRErams here
def defSP():
    # define SP and set PRErams
    sp = fsps.StellarPopulation(
                                # Constant for each run
                                dust_type=2, 
                                add_dust_emission=False, 
                                add_neb_emission=True,
                                
                                imf_type=2, #Kroupa 2001

                                sfh=1, #cons't? think so
                                const=1,
                                fburst=0,
                                tau=100,
                                tburst=11,
                                
                                zmet=13
                                )
    return sp


# = BLOCK 5 ========================================================= BLOCK 5 =
## Define PARAMS here!! 
param('dust2', 0 , 2.50)
param('dust_index', -1.0 , 1.0)
param('uvb', 0.0 , 5.0)
param('log10tage', -3.0, 0.0)


# = BLOCK 6 ========================================================= BLOCK 6 =
## Define POSTrams here
def setPOSTrams(sp):
    sp.params["dust1"] = 2.0 * sp.params["dust2"]
    return sp
    

# = BLOCK 7 ========================================================= BLOCK 7 =
## Define Drams here
nDram = 1 # Stellar Mass
dRamNames = ['log10(M_star)']


# = BLOCK 8 ========================================================= BLOCK 8 =
nDim = len(nameL)
paramRanges = np.asarray(paramRanges)


# = BLOCK 9 ========================================================= BLOCK 8 =
def updateFSPS(sp, trueNames, trueVals):
    ''' This function an Stellar Population, after setting all the variables
    in trueNames.  It accounts for variables in the format logE and log10 and
    in addition, has a feature to return the age of a SP in gigayears.  This
    last part is usefull because this is how get_spectrum() requires tage'''
    
    # set PARrams
    for pi in range(0,len(trueNames)):
        
        # the python FSPS function get_spectrum() doesnt use the  dictionary 
        # param "tage", so we log the index and account for it later
        if 'tage' in trueNames[pi]:
            tage_i = pi
            
        # set each parameter, including "tage".  "tage" is also because 
        # sp.dustmass/sp.stellarmass, unlike get_spectrum(), do use the
        # dictionary defination of "tage".  Also, we account for log10/E params
        if 'log10' in trueNames[pi]: 
            paramName = trueNames[pi][5:]
            sp.params[paramName] = 10**trueVals[pi]
        elif 'logE' in trueNames[pi]:
            paramName = trueNames[pi][4:]
            sp.params[paramName] = np.e**trueVals[pi]
        else:
            sp.params[trueNames[pi]] = trueVals[pi]
            
    # set POST rams 
    sp = setPOSTrams(sp)

    # as mentioned above, we also need tage [Gyr] because for some reason 
    # get_spec() doesn't default to sp.param['tage']. 
    if 'log10' in trueNames[tage_i]:
        tage = 10**trueVals[tage_i]
    elif 'logE' in trueNames[tage_i]:
        tage = np.e**trueVals[tage_i]
    else:
        tage = trueVals[tage_i]

    # return
    return sp, tage
    
# = BLOCK 10 ======================================================= BLOCK 10 =
sampler_type = 'Ensemble' # PT or Ensemble
run_burnin_opt = True
# fill ball_radius with ball_radius_frac of parameter value
ball_radius_frac =  0.0001 #of entire parameter range
ball_radius = np.zeros([nDim])
for i in range(0, nDim):
    ball_radius[i]= ball_radius_frac*(higherL[i]-lowerL[i])
nWalkers = 20
nSteps = 40
nBurnInSteps = 40
nTemps = 0

if sampler_type == 'PT':
    print 'WARNING : PT SAMPLER IS STILL UNDER DEVELOPMENT'
