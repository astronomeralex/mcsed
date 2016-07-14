'''
This script is run as a stand alone process; its result is saved and used
later in preperation for the emcee walk.  Before one runs this script, they
must properly set up variable_house.py.

This code contains a script to create SEDs for each point in the parameter 
space under investigation (called sedGrid).  It can also create a dust mass
and stellar mass value for each point.  

Once the user has a SED for each point in parameter space, they will generally
create an interpolation function so that emcee's walkers can walk between
parameter points.  This is done automatically in main.py

All of the grids created in this scipy have size [pSpaceGridSize*nDims,N] where
N is the size of the FSPS object in question.  For a dust/stellar mass, N=1,
but for a SED N = 1963.  Example: if you want to investigate the SEDs created
by varying four parameters on a regular grid of size 20, sedGrid would have
shape = [20,20,20,20,1963]. Note that each entry is a float, so these
can get very big.

Note that the grid is regularly shaped in the parameter dimesnions for speedy
interpolation.

The units of the SEDs are basic FSPS units: L_sol/Hz.  All SEDs generate points
at the same wavelengths, called waves, which is also saved during this script.
Again, these luminosities are converted into fluxes later.

The  only user specified variables are immediatly following the imports.

Results are saved in /data/pSpaces/
'''
# =========================================================================== #
import time
import wheres_home
import numpy as np
import variable_house as vHouse
import mcsed_utils as myTools
print "... Everything imported"

# USER SPECIFIED VARIABLES
save_option = "npz"  # "npz" or "npy" or "pickle"
print_option = True # Would you like to print updates to the screen

# =========================================================================== #
def create_all():
    
    # this takes a while so we like to keep track of how long
    start_time = time.time()
    
    # we need to know where home is in order to get to variable_home.py
    masterfilepath = wheres_home.getHomeLocation()

    # Set up stellar population based on variable house inital conditions
    sp = vHouse.defSP()
    
    # needed for cartesian(), to get all param combo's
    parameterPointsRange = vHouse.paramRanges
    parameterVarNames = vHouse.nameL

    # load vars from parameters
    regGridSize = vHouse.pSpaceGridSize
    nDim = vHouse.nDim

    # get the waves from FSPS
    waves = sp.wavelengths
    
    # create size templates
    massResultSize = []
    sedResultSize = []
    valsResultSize = []
    
    for i in range(0,nDim):
        massResultSize.append(regGridSize)
        sedResultSize.append(regGridSize)
        valsResultSize.append(regGridSize)
        
    # add extra dimension size for these 2
    sedResultSize.append(len(waves))
    valsResultSize.append(nDim)

    valsGrid = np.zeros(valsResultSize) # ignore this when iterating
    dustMassGrid = np.empty(massResultSize)
    stellarMassGrid = np.empty(massResultSize)
    sedGrid = np.empty(sedResultSize)

    print
    print 'empty grids built...'



    #See Cartesian below to understand comboParamValues
    comboParamValues = myTools.cartesian(parameterPointsRange)
    comboRows = np.shape(comboParamValues)[0]
    
    print 'FSPS initalized, will do ' + str(comboRows) + ' iterations.'
    print 'sedGrid will be in shape: ',np.shape(sedGrid)
    print 'valsGrid will be in shape: ',np.shape(valsGrid)
    print 'dustMassGrid will be in shape: ',np.shape(dustMassGrid)
    print 'stellarMassGrid will be in shape: ',np.shape(stellarMassGrid)
    print '         -=- BEGIN -=- '
                                
                                
    # The meat of the function
    for row_idx in range(0,comboRows):
        
        # Trick to get parameter idx from actual parameter value, cost of using
        # cartesian.
        pVals = comboParamValues[row_idx,:]
        pIndexs = np.zeros(nDim)
        
        for i in range(0,nDim):
            # get corrisponding index of parameter value
            pIndexs[i] = int(np.where(parameterPointsRange[i,:] == pVals[i])[0])

        # call math_function function that sets up a stellar population
        sp, tage = vHouse.updateFSPS(sp, parameterVarNames, pVals)
        waves, spec = sp.get_spectrum(tage=tage)
        
        # because we use an array to index an array, we need to make it a tuple
        pTuple = tuple(pIndexs)

        #valsGrid[pTuple] = comboParamValues[row_idx]
        dustMassGrid[pTuple] = sp.dust_mass
        stellarMassGrid[pTuple] = sp.stellar_mass
        sedGrid[pTuple] = spec
    
        if print_option == True:
            print '...% done: ' + str( 100* ((float(row_idx+1))/comboRows) )  
    
    print 'valsGrid shape: ', np.shape(valsGrid)
    print 'dustMassGrid shape: ', np.shape(dustMassGrid)
    print 'stellarMassGrid shape: ', np.shape(stellarMassGrid)
    print 'sedGrid shape: ', np.shape(sedGrid)
    print 'parameterValues created in: ', str(time.time() - start_time)
    print '         -=- DONE -=- '
    print 'Press enter to save. WARNING: This will overwrite any pSpaces with'
    print 'the same name in data/active/pSpaces/'
    raw_input()
    


    if save_option == 'npy':
        save_time = time.time()
        np.save(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_valsGrid.npy', valsGrid)
        np.save(masterfilepath+'data/pSpaces/waves.npy', waves)
        np.save(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_dustMassGrid.npy', dustMassGrid)
        np.save(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_stellarMassGrid.npy', stellarMassGrid)
        print 'saved using np.save in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.npy'
        print '         -=- SAVED -=- '
        
    if save_option == 'npz':
        save_time = time.time()
        np.savez_compressed(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_other.npz',
                            valsGrid=valsGrid, 
                            waves=waves,
                            dustMassGrid=dustMassGrid,
                            stellarMassGrid=stellarMassGrid)
        print 'saved using np.savez_compressed in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_other.npz'
        print '         -=- SAVED -=- '
        
    if save_option == 'pickle':
        save_time = time.time()
        import cPickle as pickle
        pickle.dump(valsGrid, open(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_valsGrid.pickle', 'wb'))
        pickle.dump(dustMassGrid, open(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_dustMassGrid.pickle', 'wb'))
        pickle.dump(stellarMassGrid, open(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_stellarMassGrid.pickle', 'wb'))
        pickle.dump(waves, open(masterfilepath+'data/pSpaces/waves.pickle', 'wb'))
        print 'saved using cPickle.dump() in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.pickle'
        print '         -=- SAVED -=- '
        
        
    
    # Save SED GRID
    if save_option == 'npy':
        save_time = time.time()
        np.save(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.npy', sedGrid)
        print 'saved using np.save in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.npy'
        print '         -=- SAVED -=- '
        
    if save_option == 'npz':
        save_time = time.time()
        np.savez_compressed(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sed.npz',
                            sedGrid=sedGrid)
        print 'saved using np.savez_compressed in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sed.npz'
        print '         -=- SAVED -=- '
        
    if save_option == 'pickle':
        save_time = time.time()
        import cPickle as pickle
        pickle.dump(sedGrid, open(masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.pickle', 'wb'))
        print 'saved using cPickle.dump() in ,' , str(time.time() - save_time)
        print 'saved @: '+masterfilepath+'data/pSpaces/size'+str(vHouse.pSpaceGridSize)+'_sedGrid.pickle'
        print '         -=- SAVED -=- '
        
        

# =============================================================================
create_all()