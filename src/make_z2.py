# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 17:39:48 2015

@author: hunterbrooks

OUTPUT UNITS: CGS

CURRENT CONFIGURATION: currently set up to run with all of the goods south
bands and all of the photom south bands, except 850LPc and 606Wc since 850LP
and 606W are included already.  This setup is a tad different than 
http://iopscience.iop.org/0067-0049/214/2/24/pdf/apjs_214_2_24.pdf table 3
because I've added f70. 41 bands total. Uses photom IRAC bands over goods 
IRAC bands because if memory serves, they were a tad different.

"""
# ======================================================================= #
import numpy as np
import pyfits
import wheres_home
home = wheres_home.getHomeLocation()

# =========================================================================== #
verbose = True

# =========================================================================== #
# we do some Ra/Dec matching in this script so set the thresholds
goodsThresh = 2 #arcsec
swiftThresh = 2 #arcsec
''' The following two lines are a change to the natural position of swift
ra's and dec's because Greg told me that the person he got this data from
said that the previous scientist did this.  See an email from greg for more detials.
added 1.0 arcsec to the RA and subtracted 0.7 arcsec from the Dec.
 1 arc sec = 0.000277777778 degree'''
swift_extra_ra = 1.0 * 0.000277777778
swift_extra_dec = -0.7 *0.000277777778

# =========================================================================== #
# find nFilters dealing with
col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')
nFilters = len(col_names)

# =========================================================================== #
# load data info, but not actually average fluxes, like redshift/ra/dec
main = np.loadtxt(home+'data/z2/data_house/all_pos_beta.dat')
nObjects = np.shape(main)[0]

# =========================================================================== #
# load Henry's flags data
flags = np.loadtxt(home+'data/z2/data_house/flags.dat')

# =========================================================================== #
# function to get distance between two sets of ra/dec on the night sky
def getApproxDist(ra1, dec1, ra2, dec2):
    '''Inputs: degrees'''
    deltaRA = ra2 - ra1
    deltaDEC = dec2 - dec1
    dist = np.sqrt( (np.cos(np.pi/180.0 * dec1) * deltaRA)**2 + (deltaDEC)**2 )
    return dist

# =========================================================================== #
def checkflags(infield, inid):
    row = -1 #place holder
    flagsfields = flags[:,0]
    flagsids = flags[:,1]
    
    for i in range(0, np.shape(flagsfields)[0]): #number of objects in flags
        if infield == flagsfields[i]:
            if inid == flagsids[i]:
                row = i
   # assert row != -1, "Assure a match for this ra/dec is found in Henry's list
    if row == -1:
        print 'WARNING: Not in henry list: ', infield, inid
        return False
        
    else:
        for i in range(1,6): #allows for -1, -2.. -5
            if flags[row,-i] == 1:
                if verbose == True:
                    print "... This object failed henry's flag test"
                return False
            
    return True

# =========================================================================== #
class fitscatalog:
    def __init__(self, name, location, technames, filternames, othernames):

        self.name = name 
        self.loc = location
        # tech names is the actual fits string you must query
        self.tnames = technames
        # filter names is the MCSED equivilent in filter_names.dat
        self.fnames = filternames
        # ra/dec among other things
        self.onames = othernames
        
        # confirm tech and filter names are same
        assert self.tnames == self.fnames, 'Make sure each tech name has its equivilent filter name'
        
        # make data easily accessible
        self.dat = pyfits.open(self.loc)[1].data
        
    

    def distbycoord(self, inra, indec):
        # find nearest match by calcing distance accross sky
        catara = self.dat['ra']
        catadec = self.dat['dec']
        distances = np.zeros(len(catara))
        
        # calculate distance for each object in catalog
        for i in range(0, len(distances)):
            distances[i] = getApproxDist(inra, indec, catara[i], catadec[i])
        # return index of min of distances list
        # make sure only 1 closest match is found by convert to int
        closestidx = int(np.where(distances == np.min(distances))[0])
        
        # return (distance to cloest object[arcsec], index of closest object)
        return distances[closestidx]*3600, closestidx
            
            
    def getdataforobject(self, idx, requestedfilters):
        # request data for index=idx in bands=requestedfilters, if fits
        # file doesnt have a filter, 0 is returned in that spot
        fluxes = np.zeros(len(requestedfilters))/10**(29) #convert from uJy to CGS
        errors = np.zeros(len(requestedfilters))/10**(29) #_ima and _sim
        for i in range(0, len(requestedfilters)):
            try: 
                filterdataforall = self.dat[ requestedfilters[i] ]
                fluxes[i] = filterdataforall[idx]
                errbasename = requestedfilters[i][1:]
                errors[i] = np.max( self.dat[ errbasename+'_sim'], 
                                    self.dat[ errbasename+'_ima'])
            except KeyError:
                # couldn't find requested filter name in catalog filter options
                continue
        return fluxes, errors
            



           
class txtcatalog:
    def __init__(self, name, location, filternames, filterindexes,
                 othernames, otherindexes, skiprows):
        self.name = name
        self.loc = location
        # names and indexes of the filters, same as in filter_names.dat
        self.fnames = np.asarray(filternames, dtype='str')
        self.fidxs = filterindexes
        # just any other information
        self.onames = np.asanyarray(othernames, dtype='str')
        self.oidxs = otherindexes
        # skip rows at begining of text doc?
        self.skiprows = skiprows
        
        # confirm filter names and indexes are same length
        assert len(self.fnames) == len(self.fidxs)
 
        self.dat = np.loadtxt(self.loc, skiprows=self.skiprows)
        
        # purge -99.0's from data
        ax1_idxs, ax2_idxs = np.where(self.dat == -99.0)
        for i in range(0,len(ax1_idxs)):
            self.dat[ ax1_idxs[i],ax2_idxs[i] ] = 0

            
            
            
    def matchbyfield(self, infield, inid):
        index = None # placeholder
        catafields = self.dat[:,0]
        cataids =  self.dat[:,1]
        for i in range(0, np.shape(self.dat)[0]): #number of objects in self.dat
            if infield == catafields[i]:
                if inid == cataids[i]:
                    index = i
        return index
                    
                    
    
    def getdataforobject(self, row, requestedfilters):
        requestedfilters = np.asarray(requestedfilters, dtype='str')
        fluxes = np.zeros(len(requestedfilters))
        errors = np.zeros(len(requestedfilters))
        
        for i in range(0, len(requestedfilters)):
            try:
                indexinindexarray = int(np.where(self.fnames == requestedfilters[i])[0])
                indexindata = self.fidxs[indexinindexarray]
                if 'photom' in self.name:
                    # include math factor, [inverse of 10^(0.4(25-23.9))] to be in
                    # in uJy, the put into CGS
                    factor = 0.3630 / 10**(29)
                if 'COSMOS' in self.name:
                    factor = 0.3630 / 10**(29) # assuming this is the case
                fluxes[i] = self.dat[row, indexindata]*factor
                errors[i] = self.dat[row, indexindata+1]*factor #next col is err
                
            except TypeError:
                # the requested filter string isn't in self.fnames
                continue
            
        return fluxes, errors


# =========================================================================== #
gn = fitscatalog('goods_north',
                 home+'data/z2/data_house/GH_GOODS-North_tables_DR1.fits',
                 ['f3p6','f4p5','f5p8','f8p0','f24','f70','f100','f160','f250','f350','f500'],
                 ['f3p6','f4p5','f5p8','f8p0','f24','f70','f100','f160','f250','f350','f500'],
                 ['- placeholder -'])
                 
gs = fitscatalog('goods_south',
                 home+'data/z2/data_house/GH_GOODS-South_tables_DR1.fits',
                 ['f3p6','f4p5','f5p8','f8p0','f24','f70','f100','f160'],
                 ['f3p6','f4p5','f5p8','f8p0','f24','f70','f100','f160'],
                 ['- placeholder -'])

pn = txtcatalog('photom_north',
                home+'data/z2/data_house/goods_n_photom_final.dat',
                ['U_KPNO','F435W','BJ_subaru','g_keck','VJ_subaru','F606W','rp_subaru','Rs_keck','i_subaru','F775W','zp_subaru','F850LP','F125W','J_subaru','F140W','F160W','H_subaru','Ks_subaru','IRAC1','IRAC2','IRAC3','IRAC4'],
                range(3,45+1,2), #starting @ idx=3, to 47, every other one
                ['field','id','z'],
                [0,1,2],
                1)

ps = txtcatalog('photom_south',
                home+'data/z2/data_house/goods_s_photom_final.dat',
                ['U38','B_eso','V_eso','Rc_eso','I_eso','U_VLT','R_VLT','F435W','F606Wc','F606W','F775W','F814Wc','F850LPc','F850LP','J_isaac','J_wirc','F125W','F140W','F160W','H_isaac','K_wirc','Ks_isaac','IA427','IA445','IA505','IA527','IA550','IA574','IA598','IA624','IA651','IA679','IA738','IA767','IA797','IA856','IRAC1','IRAC2','IRAC3','IRAC4'],
                range(3,81+1,2), #starting @ idx=3, to 47, every other one
                ['field','id','z'],
                [0,1,2],
                1)

cos = txtcatalog('COSMOS',
                home+'data/z2/data_house/cosmos_photom.dat',
                ['u_cfht','BJ_subaru','g_cfht','VJ_subaru','F606Wc','r_cfht','rp_subaru','i_cfht','i_subaru','F814W','z_cfht','zp_subaru','y_eso','J1_nfirm','J2_nfirm','J3_nfirm','F125W','J_wirc','J_eso','F140W','F160W','H1_nfirm','H2_nfirm','H_wirc','H_eso','K_nfirm','Ks_wirc','Ks_eso','IA427','IA464','IA484','IA505','IA527','IA574','IA624','IA679','IA709','IA738','IA767','IA827','IRAC1','IRAC2','IRAC3','IRAC4'],
                range(3,89+1,2), #starting @ idx=3, to 47, every other one
                ['field','id','z'],
                [0,1,2],
                1)                
                

# =========================================================================== #
# create output array, called super.dat, shape per object is:
outputFluxes = np.zeros([nObjects, nFilters])
outputErrors = np.zeros([nObjects, nFilters])
outputRedshift = np.zeros([nObjects])
outputInfo = np.zeros([nObjects, 7])





# =========================================================================== #
# start loop, making nEntires fake objects
for obji in range(0,nObjects):
    

    
    # ======================================================================= #
    # fill info by re arranging main
    outputInfo[obji,0] = main[obji, 2] #redshift
    outputInfo[obji,1] = main[obji, 3] #ra
    outputInfo[obji,2] = main[obji, 4] #dec
    outputInfo[obji,3] = main[obji, 0] #field number
    outputInfo[obji,4] = main[obji, 1] #ID
    outputInfo[obji,5] = main[obji, 5] #slope
    outputInfo[obji,6] = main[obji, 6] #slope error
    
    # fill output Redshift
    outputRedshift[obji] = main[obji, 2]
    
    # get field and ID to be used to match to photom
    field = int(main[obji,0])
    mainID = int(main[obji, 1])

    # get ra and dec of object, used to find closest match
    ra = main[obji,3]
    dec = main[obji,4]
    
    if verbose == True:    
        print '-----------------', obji
        print ' Working Object: ', int(field), int(mainID)
        print 'ra/dec: ', ra , dec
        
    # ======================================================================= #
    # check that none of Henry's flags are up
    if checkflags(field, mainID) == False:
        continue
    
    
    # ======================================================================= #
    '''
    THIS SECTION IS CURRENTLY SKIPPED BECAUSE ALL THE DATA IM USING IS IN THE
    PHOTOM CATALOG.
    
    # get the closest goods index
    gNcloseDist, gNcloseIdx = gn.distbycoord(ra, dec)
    gScloseDist, gScloseIdx = gs.distbycoord(ra, dec)
    
    
    # fill goodsfluxes with closest
    if gNcloseDist > gScloseDist:
        goodsfluxes, goodserrors = gs.getdataforobject(gScloseIdx, requestedfilters)
        if verbose == True:
            print '... min goods distance is in goods SOUTH'
            print '... min goods distance: ', gScloseDist
            print '... ... goods south index: ', gScloseIdx
        
    if min_distto_gN < min_distto_gS:
        goodsfluxes, goodserrors = gn.getdataforobject(gNcloseIdx, requestedfilters)
        if verbose == True:
            print '... min goods distance is in goods NORTH'
            print '... min goods distance: ', gNcloseDist
            print '... ... goods north index: ', gNcloseIdx
    '''
    
    # ======================================================================= #
    # get closest photom index
    pNindex = pn.matchbyfield(field, mainID)
    pSindex = ps.matchbyfield(field, mainID)
    """
    photomfilters = [   'IRAC1',
                        'IRAC2',
                        'J_wirc',
                        'K_wirc',
                        'J_isaac',
                        'H_isaac',
                        'Ks_isaac',
                        'F125W',
                        'F160W',
                        'F140W',
                        'F606W',
                        'F435W',
                        'F775W',
                        'F850LP',
                        'IA427',
                        'IA445',
                        'IA505',
                        'IA527',
                        'IA550',
                        'IA574',
                        'IA598',
                        'IA624',
                        'IA651',
                        'IA679',
                        'IA738',
                        'IA767',
                        'IA797',
                        'IA856',
                        'U38',
                        'B_eso',
                        'V_eso',
                        'Rc_eso',
                        'I_eso',
                        'U_VLT',
                        'R_VLT',
                        'U_KPNO',
                        'BJ_subaru',
                        'VJ_subaru',
                        'rp_subaru',
                        'Rs_keck',
                        'i_subaru',
                        'zp_subaru',
                        'J_subaru',
                        'H_subaru',
                        'Ks_subaru']
    """

    if pNindex != None and pSindex == None:
        photomfluxes, photomerrors = pn.getdataforobject(pNindex, col_names)
        if verbose == True:
            print '... match found in photom NORTH'
            print '... ... pNindex: ', pNindex
        # actually save to array
        outputFluxes[obji,:] = photomfluxes      
        outputErrors[obji,:] = photomerrors
        continue
        
    if pNindex == None and pSindex != None:
        photomfluxes, photomerrors = ps.getdataforobject(pSindex, col_names)
        if verbose == True:
            print '... match found in photom SOUTH'
            print '... ... pSindex: ', pSindex
         # actually save to array
        outputFluxes[obji,:] = photomfluxes      
        outputErrors[obji,:] = photomerrors
        continue
        
    if pNindex == None and pSindex == None:
        if verbose == True:
            print '... Photom match not found'
    
    
    
    
    # ======================================================================= #
    # get closest COSMOS index
    COSindex = cos.matchbyfield(field, mainID)
        
    if COSindex != None:
        COSfluxes, COSerrors = pn.getdataforobject(COSindex, col_names)
        if verbose == True:
            print '... match found in COSMOS'
            print '... ... COSindex: ', COSindex
        # actually save to array
        outputFluxes[obji, :] = COSfluxes      
        outputErrors[obji, :] = COSerrors
        continue
        
    else:
        if verbose == True:
            print '... COSMOS match not found'
     

# =========================================================================== #
# trim duplicate observations, as given by Greg
dup_fields = [73,66,66,60,63,63,52,60,63,62,58,58,62,58,86]
dup_ids = [2347,3257,5168,247,1756,2814,690,2442,5562,135,2430,1665,797,2592,6066]
for i in range(0, np.shape(outputFluxes)[0]):
    for ii in range(0, len(dup_fields)):
        if outputInfo[i,3] == dup_fields[ii]:
            if outputInfo[i,4] == dup_ids[ii]:
                outputFluxes[i,:] = 0
                outputErrors[i,:] = 0
# trim poor photom, as given by Greg
# DONT DO THIS YET
                
# =========================================================================== #
# make doAbleIndexes.dat
nNonZero = np.zeros(np.shape(outputFluxes)[0])
for i in range(0, np.shape(outputFluxes)[0]):
    nNonZero[i] = len(np.where(outputFluxes[i,:] != 0)[0]) 
doAbleIndexes =  np.where(nNonZero != 0)[0]                
                
# =========================================================================== #
# Save outputSuperDat
np.savetxt(home+'data/z2/superFluxes.dat', outputFluxes)
np.savetxt(home+'data/z2/superErrors.dat', outputErrors)
np.savetxt(home+'data/z2/superRedshift.dat', outputRedshift)
np.savetxt(home+'data/z2/superInfo.dat', outputInfo)
np.savetxt(home+'data/z2/superIndexes.dat', doAbleIndexes)

