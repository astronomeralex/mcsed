'''
Creates a z2.dat file that contains all the relevent information for this bulk
of data. 

units = uJy
'''
import numpy as np
import pyfits

# NOTE IRAC1,2,3,4 in the 3rd block are always put into IRAC1,2,3,4 in the 4th block
colNames = np.array([ '3p6',
                      '4p5',
                      '5p8',
                      '8p0', 
                      'f24',
                      'f70',
                      'f100',
                      'f160',
                      'f250',
                      'f350',
                      'f500',
                      
                      'w2',
                      'm2',
                      'w1',
                      'uu',
                      
                      'U_KPNO',
                      'F435W',
                      'BJ_subaru',
                      'g_keck',
                      'VJ_subaru',
                      'F606W',
                      'rp_subaru',
                      'Rs_keck',
                      'i_subaru',
                      'F775W',
                      'zp_subaru',
                      'F850LP',
                      'F125W',
                      'J_subaru',
                      'F140W',
                      'F160W',
                      'H_subaru',
                      'KS_subaru',
                      'IRAC1',
                      'IRAC2',
                      'IRAC3',
                      'IRAC4',
                      
                      'U38',
                      'B_eso',
                      'V_eso',
                      'Rc_eso',
                      'I_eso',
                      'U_VLT',
                      'R_VLT',
                      'F435W',
                      'F606Wc',
                      'F606W',
                      'F775W',
                      'F814Wc',
                      'F850LPc',
                      'F850LP',
                      'J_isaac',
                      'J_wirc',
                      'F125W',
                      'F140W',
                      'F160W',
                      'H_isaac',
                      'K_wirc',
                      'Ks_isaac',
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
                      'IRAC1',
                      'IRAC2',
                      'IRAC3',
                      'IRAC4'])
                      
                      
# =========================================================================== #
goodsThresh = 2 #arcsec
swiftThresh = 2 #arcsec

''' added 0.9 arcsec to the RA and subtracted 0.5 arcsec from the Dec'''
# arc sec = 0.000277777778 degree
swift_extra_ra = 1.0 * 0.000277777778
swift_extra_dec = -0.7 *0.000277777778

#parameter to print everything
verbose = False


'''
def getDist(ra1, dec1, ra2, dec2):
    return np.sqrt( ((ra1-ra2)*np.cos(dec1))**2 + (dec1 - dec2)**2)
'''
def getDist(ra1, dec1, ra2, dec2):
    
    #convert to radians
    ra1 = np.deg2rad(ra1)
    dec1 = np.deg2rad(dec1)
    ra2 = np.deg2rad(ra2)
    dec2 = np.deg2rad(dec2)
    
    #cos(A) = sin(d1)sin(d2) + cos(d1)cos(d2)cos(ra1-ra2)
    return np.rad2deg( np.arccos( np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2) ) )


# =========================================================================== #
dPath = '/Users/hunterbrooks/Astronomy/modelmatching/data/z2/data_house/'

main = np.loadtxt(dPath+'all_pos_beta.dat')
#shape=(323,7) : field number, ID, redshift, RA, Dec, slope, slope_err


# =========================================================================== #
outputColumns = 7 + 2 + 2 + 2 + 2 
output = np.zeros([np.shape(main)[0], outputColumns])



superColumns = outputColumns + (11)*2  + (4)*2 + (22+40)*2
# 15 raw columns, idx=0-14
# 11 goods filters, idx=15-36, south only has 8
# 4 swift filters, idx=37-44
# 22 photom north, idx = 45 - 88
# 40 photom south, idx = 89 - 168

superinfo = np.zeros([np.shape(main)[0], superColumns])
# GOODS: 3p6,4p5,5p8,8p0, 24, 70 ,100, 160, (250, 350, 500) <- south doesnt have last 3
# SWIFT: w2, m2, w1, uu
# PHOTOM North: U_KPNO, F435W, BJ_subaru, g_keck, VJ_subaru, F606W, rp_subaru, Rs_keck,
#               i_subaru, F7775W, zp_subaru, F850LP, F125W, J_subaru, F140W, F160W, H_subaru,
#              KS_subaru, IRAC 1, IRAC 2, IRAC 3, IRAC 4
# PHOTOM South: u38, b_eso, v_eso, rc_eso, i_eso, U_vlt, R_vlt, F435W F606Wc, F606W
#               F775W, f814Wc, F850LPc, F850LP, J_isaac, J_wirc, F125W, F140W, F160W,
#              H_isaac, K_wirc, Ks_isaac, IA427, IA445, IA505, IA527, IA550, IA574
#             IA598, IA624, IA651, IA679, IA738, IA767, IA797, IA856, IRAC 1, IRAC 2, IRAC3, IRAC4



testing = np.zeros([np.shape(main)[0]])

# column titles of outputColumns below:
'''
field,    1 (idx = 0)
id,
redshift, 
ra, 
dec,    5
slope, 
slope_err, 
goods_n idx, 
goods_n dist,
goods s idx,  10 
goods_s dist, 
swift_idx, 
swift distance,
which_photom,       #1 = north, 2 = south, 0 = none
photom_idx       15


'''
# =========================================================================== #
gN = pyfits.open(dPath + 'GH_GOODS-North_tables_DR1.fits')
gNdata = gN[1].data # 
gS = pyfits.open(dPath + 'GH_GOODS-South_tables_DR1.fits')
gSdata = gS[1].data

'''
  f3p6           uJy    IRAC 3.6 um flux density
  err3p6         uJy    Error on IRAC 3.6 um flux density
  
  f4p5           uJy    IRAC 4.5 um flux density
  err4p5         uJy    IRAC 4.5 um flux density
  
  f5p8           uJy    IRAC 5.8 um flux density
  err5p8         uJy    Error on IRAC 5.8 um flux density
  
  f8p0           uJy    IRAC 8.0 um flux density
  err8p0         uJy    Error on IRAC 8.0 um flux density
  
  f24            uJy    MIPS 24 um flux density
  err24_ima      uJy    MIPS 24 um flux error on residual map
  err24_sim      uJy    MIPS 24 um flux error on Monte-Carlo simulations
  
  f70            uJy    MIPS 70 um flux density
  err70_ima      uJy    MIPS 70 um flux error on residual map
  err70_sim      uJy    MIPS 70 um flux error on Monte-Carlo simulations
  
  f100           uJy    PACS 100 um flux density
  err100_ima     uJy    PACS 100 um flux error on residual map
  err100_sim     uJy    PACS 100 um flux error on Monte-Carlo simulations
  
  f160           uJy    PACS 160 um flux density
  err160_ima     uJy    PACS 160 um flux error on residual map
  err160_sim     uJy    PACS 160 um flux error on Monte-Carlo simulations
  
  f250           uJy    SPIRE 250 um flux density
  err250_ima     uJy    SPIRE 250 um flux error on residual map
  err250_sim     uJy    SPIRE 250 um flux error on Monte-Carlo simulations
  
  f350           uJy    SPIRE 350 um flux density
  err350_ima     uJy    SPIRE 350 um flux error on residual map
  err350_sim     uJy    SPIRE 350 um flux error on Monte-Carlo simulations
  
  f500           uJy    SPIRE 500 um flux density
  err500_ima     uJy    SPIRE 500 um flux error on residual map
  err500_sim     uJy    SPIRE 500 um flux error on Monte-Carlo simulations
  
  SOUTH DOESNT HAVE LAST 3 FILTERS
'''


# =========================================================================== #
# THESE ID's are for the all pos beta . dat
pNdata = np.loadtxt(dPath+'goods_n_photom_final.dat', 
                                    skiprows=1)
# remove -99s
temp = np.zeros(np.shape(pNdata))
for rowIDX in range(0,np.shape(pNdata)[0]):
    for colIDX in range(0,np.shape(pNdata)[1]):
        if pNdata[rowIDX, colIDX] != -99.0:
            temp[rowIDX,colIDX] = pNdata[rowIDX, colIDX]
pNdata = temp
'''
Fld   
ID     
z    
U_KPNO      err    1
F435W       err 
BJ_subaru   err 
g_keck      err 
VJ_subaru   err    5
F606W        err
rp_subaru    err 
Rs_keck      err 
i_subaru     err  
F775W        err    10
zp_subaru    err 
F850LP      err   
F125W       err 
J_subaru     err 
F140W       err     15
F160W       err 
H_subaru    err 
Ks_subaru    err  
IRAC1       err    
IRAC2       err    20
IRAC3       err   
IRAC4       err    22
'''
pSdata = np.loadtxt(dPath+'goods_s_photom_final.dat', 
                                    skiprows=1)
# remove -99s
temp = np.zeros(np.shape(pSdata))
for rowIDX in range(0,np.shape(pSdata)[0]):
    for colIDX in range(0,np.shape(pSdata)[1]):
        if pSdata[rowIDX, colIDX] != -99.0:
            temp[rowIDX,colIDX] = pSdata[rowIDX, colIDX]
pSdata = temp
'''
Fld   
ID    
z   
U38     err    1
B_eso    err 
V_eso    err 
Rc_eso   err   
I_eso    err   5
U_VLT    err    
R_VLT    err  
F435W   err   
F606Wc    err  
F606W   err    10
F775W    err  
F814Wc    err 
F850LPc   err 
F850LP   err  
J_isaac   err   15
J_wirc    err   
F125W    err 
F140W    err 
F160W    err   
H_isaac  err   20
K_wirc    err   
Ks_isaac  err
IA427    err 
IA445    err  
IA505   err     25
IA527    err
IA550    err 
IA574    err 
IA598    err 
IA624    err   30
IA651    err   
IA679    err 
IA738    err  
IA767    err  
IA797    err    35
IA856    err    
IRAC1   err   
IRAC2    err  
IRAC3    err   
IRAC4    err   40
''' # 40 bands
# =========================================================================== #
s = pyfits.open(dPath + 'nphot_uselect.fits')
sAdata = s[1].data #swift all data, no south/north stuff

# =========================================================================== #
for betaIDX in range(0,np.shape(main)[0]):
#for betaIDX in range(239,241):
    
    
    field = main[betaIDX,0]
    mainID = main[betaIDX, 1]

    ra = main[betaIDX,3]
    dec = main[betaIDX,4]
    
    output[betaIDX,0] = field
    output[betaIDX,1] = mainID
    output[betaIDX,2] = main[betaIDX,2]
    output[betaIDX,3] = ra
    output[betaIDX,4] = dec
    output[betaIDX,5] = main[betaIDX,5] #slope
    output[betaIDX,6] = main[betaIDX,6] #slope err
    
    if verbose == True:
        print str(betaIDX) + ' ('+str(int(field))+'-'+str(int(mainID))+ '):'
    
    
    
    # get closest for GOODS
    gN_ra = gNdata['ra'] #2710
    gS_ra = gSdata['ra'] #2531
    
    gN_dec = gNdata['dec']
    gS_dec = gSdata['dec']
    
    
    closest_gN = np.zeros(np.shape(gN_ra))
    closest_gS = np.zeros(np.shape(gS_ra))

    #SANITY CHECK: np.sqrt( (0-90)**2 + np.cos(0-0)**2 ) = 90 degrees
    for nIDX in range(0,np.shape(gN_ra)[0]):
        closest_gN[nIDX] = getDist(ra, dec, gN_ra[nIDX], gN_dec[nIDX]) #in degress
        
    for nIDX in range(0,np.shape(gS_ra)[0]):
        closest_gS[nIDX] = getDist(ra, dec, gS_ra[nIDX], gS_dec[nIDX])#in degress
        

    closest_gN_idx = int(np.where(closest_gN == np.min(closest_gN))[0])
    closest_gS_idx = int(np.where(closest_gS == np.min(closest_gS))[0])
    
    output[betaIDX,7] = closest_gN_idx#n idx
    output[betaIDX,8] = closest_gN[closest_gN_idx] * (3600) #n distance in arcseconds
    if verbose == True:
        print '... min goods north distance ', output[betaIDX,8]
    
    output[betaIDX,9] = closest_gS_idx#s idx
    output[betaIDX,10] = closest_gS[closest_gS_idx] * (3600) #ssouth distance in arcseconds
    if verbose == True:
        print '... min goods south distance ', output[betaIDX,10]
    
    # get closest from swift
    s_ra = sAdata['ra']
    s_dec = sAdata['dec']
    
    closest_s = np.zeros(np.shape(s_ra))
    
    for nIDX in range(0,np.shape(s_ra)[0]):
        closest_s[nIDX] = np.sqrt( (ra-s_ra[nIDX]+swift_extra_ra)**2 + np.cos(dec-s_dec[nIDX]+swift_extra_dec)**2 )
    
    closest_s_idx = int(np.where(closest_s == np.min(closest_s))[0])
    
    output[betaIDX,11] = closest_s_idx#s idx
    output[betaIDX,12] = closest_s[closest_s_idx]  * (3600)#s distance in arcseconds
    if verbose == True:
        print '... min swift distance ', output[betaIDX,12]
    
    
    
    
    # get photom data
    pNfield = pNdata[:,0]
    pNid = pNdata[:,1]
    
    pSfield = pSdata[:,0]
    pSid = pSdata[:,1]
    
    
    for pIDX in range(0,len(pNfield)):
        if pNfield[pIDX] == field:
            if pNid[pIDX] == mainID:
                output[betaIDX,13] = 1
                output[betaIDX,14] = pIDX
                if verbose == True:
                    print '... match found in photom north @ ', pIDX
                
    for pIDX in range(0,len(pSfield)):
        if pSfield[pIDX] == field:
            if pSid[pIDX] == mainID:
                output[betaIDX,13] = 2
                output[betaIDX,14] = pIDX
                if verbose == True:
                    print '... match found in photom south @ ', pIDX
                
                
    # BEGIN THE ACTUAL SAVING OF DATA
    superinfo[betaIDX,0] = field
    superinfo[betaIDX,1] = mainID
    superinfo[betaIDX,2] = main[betaIDX,2]
    superinfo[betaIDX,3] = ra
    superinfo[betaIDX,4] = dec
    superinfo[betaIDX,5] = main[betaIDX,5] #slope
    superinfo[betaIDX,6] = main[betaIDX,6] #slope err
    superinfo[betaIDX,7] = closest_gN_idx#n idx
    superinfo[betaIDX,8] = closest_gN[closest_gN_idx] * (3600) #n distance in arcseconds
    superinfo[betaIDX,9] = closest_gS_idx#s idx
    superinfo[betaIDX,10] = closest_gS[closest_gS_idx] * (3600) #ssouth distance in arcseconds
    superinfo[betaIDX,11] = closest_s_idx#s idx
    superinfo[betaIDX,12] = closest_s[closest_s_idx]  * (3600)#s distance in arcseconds
    superinfo[betaIDX,13] = output[betaIDX,13]
    superinfo[betaIDX,14] = output[betaIDX,14]
    
    
    # if north is closer than thresh
    if (closest_gN[closest_gN_idx] * (3600)) <= goodsThresh:# 3p6,4p5,5p8,8p0, 24, 70 ,100, 160, 250, 350, 500
        if verbose == True:
            print '... pulling data from goods north' 
        superinfo[betaIDX,15] = gNdata['f3p6'][closest_gN_idx]
        superinfo[betaIDX,16] = gNdata['err3p6'][closest_gN_idx]
        
        superinfo[betaIDX,17] = gNdata['f4p5'][closest_gN_idx]
        superinfo[betaIDX,18] = gNdata['err4p5'][closest_gN_idx]
        
        superinfo[betaIDX,19] = gNdata['f5p8'][closest_gN_idx]
        superinfo[betaIDX,20] = gNdata['err5p8'][closest_gN_idx]
        
        superinfo[betaIDX,21] = gNdata['f8p0'][closest_gN_idx]
        superinfo[betaIDX,22] = gNdata['err8p0'][closest_gN_idx]
        
        superinfo[betaIDX,23] = gNdata['f24'][closest_gN_idx]
        superinfo[betaIDX,24] = np.max(  [gNdata['err24_ima'][closest_gN_idx] , gNdata['err24_sim'][closest_gN_idx]] )
        
        superinfo[betaIDX,25] = gNdata['f70'][closest_gN_idx]
        superinfo[betaIDX,26] = np.max([ gNdata['err70_ima'][closest_gN_idx] , gNdata['err70_sim'][closest_gN_idx] ])
        
        superinfo[betaIDX,27] = gNdata['f100'][closest_gN_idx]
        superinfo[betaIDX,28] = np.max([ gNdata['err100_ima'][closest_gN_idx] , gNdata['err100_sim'][closest_gN_idx] ])
        
        superinfo[betaIDX,29] = gNdata['f160'][closest_gN_idx]
        superinfo[betaIDX,30] = np.max([ gNdata['err160_ima'][closest_gN_idx] , gNdata['err160_sim'][closest_gN_idx] ])
        
        superinfo[betaIDX,31] = gNdata['f250'][closest_gN_idx]
        superinfo[betaIDX,32] = np.max([ gNdata['err250_ima'][closest_gN_idx] , gNdata['err250_sim'][closest_gN_idx] ])
        
        superinfo[betaIDX,33] = gNdata['f350'][closest_gN_idx]
        superinfo[betaIDX,34] = np.max([ gNdata['err350_ima'][closest_gN_idx] , gNdata['err350_sim'][closest_gN_idx] ])
        
        superinfo[betaIDX,35] = gNdata['f500'][closest_gN_idx]
        superinfo[betaIDX,36] = np.max([ gNdata['err500_ima'][closest_gN_idx] , gNdata['err500_sim'][closest_gN_idx] ])
        # GOODS UNITS = uJy
 
 
     # if SOUTH is closer than thresh
    if (closest_gN[closest_gN_idx] * (3600)) <= goodsThresh:
        if verbose == True:
            print '... pulling data from goods south'
        superinfo[betaIDX,15] = gNdata['f3p6'][closest_gN_idx]
        superinfo[betaIDX,16] = gNdata['err3p6'][closest_gN_idx]
        
        superinfo[betaIDX,17] = gNdata['f4p5'][closest_gN_idx]
        superinfo[betaIDX,18] = gNdata['err4p5'][closest_gN_idx]
    
        superinfo[betaIDX,19] = gNdata['f5p8'][closest_gN_idx]
        superinfo[betaIDX,20] = gNdata['err5p8'][closest_gN_idx]
        
        superinfo[betaIDX,21] = gNdata['f8p0'][closest_gN_idx]
        superinfo[betaIDX,22] = gNdata['err8p0'][closest_gN_idx]
        
        superinfo[betaIDX,23] = gNdata['f24'][closest_gN_idx]
        superinfo[betaIDX,24] = np.max( [gNdata['err24_ima'][closest_gN_idx] , gNdata['err24_sim'][closest_gN_idx]] )
        
        superinfo[betaIDX,25] = gNdata['f70'][closest_gN_idx]
        superinfo[betaIDX,26] = np.max( [gNdata['err70_ima'][closest_gN_idx] , gNdata['err70_sim'][closest_gN_idx]] )
        
        superinfo[betaIDX,27] = gNdata['f100'][closest_gN_idx]
        superinfo[betaIDX,28] = np.max( [gNdata['err100_ima'][closest_gN_idx] , gNdata['err100_sim'][closest_gN_idx]] )
        
        superinfo[betaIDX,29] = gNdata['f160'][closest_gN_idx]
        superinfo[betaIDX,30] = np.max( [gNdata['err160_ima'][closest_gN_idx] , gNdata['err160_sim'][closest_gN_idx]] )
        # GOODS UNITS = uJy


    # ab mag = -2.5log10(flux) + 23.9 if flux is in uJy
    # flux_error = np.abs(flux * -0.4 * ln(10) * mag_err)
    if (closest_s[closest_s_idx]  * (3600)) <= swiftThresh: #w2, m2, w1, uu
        if verbose == True:
            print '... pulling data from swift'
        superinfo[betaIDX,37] = 10**((1./-2.50)*(sAdata.data['W2_ISOMAG'][closest_s_idx]-23.9))
        superinfo[betaIDX,38] = 10**(-23.9*-0.4) * np.abs( superinfo[betaIDX,37] * -0.4 * np.log(10) * sAdata.data['W2_ISOMAG_ERR'][closest_s_idx]) 
        
        superinfo[betaIDX,39] = 10**((1./-2.50)*(sAdata.data['M2_ISOMAG'][closest_s_idx]-23.9))
        superinfo[betaIDX,40] = 10**(-23.9*-0.4) * np.abs( superinfo[betaIDX,39] * -0.4 * np.log(10) * sAdata.data['M2_ISOMAG_ERR'][closest_s_idx]) 
        
        superinfo[betaIDX,41] = 10**((1./-2.50)*(sAdata.data['W1_ISOMAG'][closest_s_idx]-23.9))
        superinfo[betaIDX,42] = 10**(-23.9*-0.4) * np.abs( superinfo[betaIDX,41] * -0.4 * np.log(10) * sAdata.data['W1_ISOMAG_ERR'][closest_s_idx]) 
        
        superinfo[betaIDX,43] = 10**((1./-2.50)*(sAdata.data['UU_ISOMAG'][closest_s_idx]-23.9))
        superinfo[betaIDX,44] = 10**(-23.9*-0.4) * np.abs( superinfo[betaIDX,43] * -0.4 * np.log(10) * sAdata.data['UU_ISOMAG_ERR'][closest_s_idx]) 
        
        
        
        
        
    # PHOTOM
    #remember to adjust units, inverse of 10^(0.4(25-23.9))
    factor = 0.36307805477010134246737121236246374566858969058777151
    
    nStartIDX = 44
    if superinfo[betaIDX,13] == 1: #its in the north catalog
        if verbose == True:
            print '... pulling data from photom north'
        for super_photom_idx in range(1, (18*2) + 1):
            superinfo[betaIDX,super_photom_idx + nStartIDX ] = factor * pNdata[ superinfo[betaIDX,14] , super_photom_idx ]

        # put the IRAC filters from north into the south channels.
        for super_photom_idx in range(1, (4*2)+1):
            superinfo[betaIDX, 80 + super_photom_idx] = factor * pNdata[ superinfo[betaIDX,14] , super_photom_idx ]
            
            
    sStartIDX = 88
    if superinfo[betaIDX,13] == 2: #its in the south catalog
        if verbose == True:
            print '... pulling data from photom south'
        for super_photom_idx in range(1, (40*2) + 1):
            superinfo[betaIDX,super_photom_idx + sStartIDX ] = factor * pSdata[ superinfo[betaIDX,14] , super_photom_idx ]
        
        
    
    
    # remove goods 3.6 -> 8.0um in favor of irac bands
    if superinfo[betaIDX,15] != 0 and superinfo[betaIDX,81] != 0:
        superinfo[betaIDX,15] = 0
        superinfo[betaIDX,16] = 0
        if verbose == True:
            print '... 1', betaIDX
        
        
    if superinfo[betaIDX,17] != 0 and superinfo[betaIDX,83] != 0:
        superinfo[betaIDX,17] = 0
        superinfo[betaIDX,18] = 0
        if verbose == True:
            print '...2 ', betaIDX
        
    if superinfo[betaIDX,19] != 0 and superinfo[betaIDX,85] != 0:
        superinfo[betaIDX,19] = 0
        superinfo[betaIDX,20] = 0
        if verbose == True:
            print '...3 ', betaIDX
        
    if superinfo[betaIDX,21] != 0 and superinfo[betaIDX,87] != 0:
        superinfo[betaIDX,21] = 0
        superinfo[betaIDX,22] = 0
        if verbose == True:
            print '...4 ', betaIDX
            print '---------- <- Removed double IRAC bands'
        
    
    
    # prints number of filters per row
    testing[betaIDX] = len(np.where(superinfo[betaIDX,15:] != 0)[0])
        
# =========================================================================== #
np.savetxt('/Users/hunterbrooks/Astronomy/modelmatching/data/z2/super_helper.dat', output)
np.savetxt('/Users/hunterbrooks/Astronomy/modelmatching/data/z2/super.dat', superinfo)
np.savetxt('/Users/hunterbrooks/Astronomy/modelmatching/data/z2/columns.dat', colNames, fmt='%s')