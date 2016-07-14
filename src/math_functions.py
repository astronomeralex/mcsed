# -*- coding: utf-8 -*-
"""
Users don't need to touch this script, but if you're interested...

THIS SCRIPT CONTAINS:
    - getFlux(): used to get observed fllux from FSPS luminosity output
    - getd_L(): used to get luminosity distance as a function of redshift
    - getAvgFlux(): simply returns the average flux in each band in 
                    mcsed/data/active/filter_names.dat
    - getA(): returns the scalar value that standarizes the chi squared
              equation
              
Note: when any variable (other than A) has a capital letter in its name 
      (e.g. S_nu as opposed to s_nu) it means that it is a 2d parameter
"""
# =============================================================================
import numpy as np
import mcsed_utils as utils
import wheres_home
home = wheres_home.getHomeLocation()
col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')

# =============================================================================
# SOME CONSTANTS
c_ang = 2.998e18 #ang/s
c_km = 2.998e5 #km/s
lum_sol = 3.83e33 #erg/s
h = 0.70 #hubble uncertainty, is a %. From Hogg's Cosmography Paper
omega_m = 0.3
omega_k = 0
omega_lambda = 0.7 #these need to sum to 1!

# =============================================================================
def getFlux(lum_e, lambda_e, z):
    '''
    Inputs:
        - lum_e = luminosity straight from fsps [L_sol / Hz]
        - lambda_e = fsps wavelengths [ang]
        - z = redshift
        
    Returns:
        - s_nu = flux density in observed frame [cgs/Hz]
        - nu_obs = reddened observed frame frequency of lambda_e
    '''
    
    nu_obs = get_nu_obs(lambda_e, z)
    
    # get luminosity distance
    d_L = getd_L(z)

    #convert to flux using the d_l found above
    lum_cgs = lum_e * lum_sol
    s_nu = (1.0+z)*((lum_cgs)/(4*np.pi*d_L**2)) 
    # the 2d array, S_nu, is created and trimmed later on
    
    return s_nu, nu_obs
    
# =============================================================================
def get_nu_obs(lambda_e, z):
    '''
    Simply reddens a wavelength array based on redshift.
    '''
    nu_e = c_ang / lambda_e
    nu_obs = nu_e / (1.0+z)
    return nu_obs
    

# =============================================================================
def getd_L(z):
    ''' 
    Returns luminosity distance in cm.
    
    Input: z = redshift of object
    
    This algorithum was taken directly from Hogg's 1999 Distance measures in Cosmology
    '''
    
    # def d_H used in finding d_C
    d_H_Mpc = c_km / (100.*h) #units: [Mpc]
    d_H = d_H_Mpc * 3.086e24 #units: [cm]
    
    # def E using in finding d_C
    def E(z_prime):
        return np.sqrt(omega_m*(1+z_prime)**3+omega_k*(1+z_prime)**2+omega_lambda)
   
   # get D_c
    import scipy.integrate as spgate
    def f(z_prime):
        return (E(z_prime))**(-1.0)
    
    integral = spgate.quad(f, 0, z)[0]
    d_C = d_H * integral #units: [cm]
    
    #define d_M
    d_M = d_C #in the current omegas, units: [cm]
    
    #define d_L
    d_L = (1+z)*d_M #units: [cm]
    
    return d_L
    
# =============================================================================
def getAvgFlux( nFilters,
                S_nu,
                NU_obs,
                D_nu,
                R_nu,
                trim_opt):
    
    '''
    Returns average flux in each of the nFilters.  The capital letters in the
    input variable names signifiy that each of them has shape 
    [nFilters, len(lambda_e)] where lambda_e is the wavelength array outputted
    by FSPS.  Because the space between each nu (D_nu) is required in the sum,
    the input arrays can (if trim_option=True) have their last entires trimmed 
    by 1.  This option is handled within the function.  
    
    Inputs:
        - S_nu: ARRAY OF THE OBSERVED SED, REPEATED EACH ROW
        - NU_obs: SIMPLY THE OBSERVED FRAME FREQUENCIES
        - D_nu: DIFFERENCE BETWEEN NU_obs[idx] AND NU_obs[idx+1]
        - R_nu: TRANSMISSION ARRAY.  EACH ROW HAS A TRANSMISSION CURVE THAT
                CORRIPSONDS TO A FILTER. NORMALIZED WITH MAX OF 1. EACH POINT
                CORRIPSPONDS TO THE SAME INDEXED POINT IN NU_obs. 
        - trim_opt: TRIM THE LAST ENTRY OF THE SECOND DIMENSION OF INPUTS
                    USED TO MAKE AN INTEGRAL INTO A SUM USING THE WIDTH OF
                    MANY TINY RECTANCES OF WIDTH D_nu
                
                
    RETURNS:
        - AVERAGE FLUX IN EACH BAND
    '''

    if trim_opt == True:
        #trim all
        S_nu = S_nu[:,:-1]
        D_nu = D_nu[:,:-1]
        NU_obs = NU_obs[:,:-1]
        R_nu = R_nu[:,:-1]

    top = S_nu * R_nu * D_nu / NU_obs
    bot = 1.0  * R_nu * D_nu / NU_obs
    
    sum_top = np.sum(top, axis=1) # denominator (NU_obs) is never zero so we 
    sum_bot = np.sum(bot, axis=1) # don't need sumtool.py
    
    return sum_top/sum_bot
    
# =============================================================================
def getA(modelAvgFlux, dataAvgFlux, dataAvgFluxErr):
    '''
    MATHEMATICAL OPERATION TO RETURN STANDARIZATION FACTOR USED IN
    lnprob().
    
    Solve for this analytically by setting the derivative of Chi Squared = 0.
    '''
    
    """
    # for debugging
    print
    print modelAvgFlux
    print
    print dataAvgFlux
    print
    print dataAvgFluxErr
    print "=========================================="
    """
    

    # find standarization factor
    A_top_sum = utils.fracsum( modelAvgFlux * dataAvgFlux , dataAvgFluxErr**2 )
    A_bot_sum = utils.fracsum( modelAvgFlux**2 , dataAvgFluxErr**2)
    
    return A_top_sum/A_bot_sum
    
    
# =============================================================================
def filterMachine(f, er, z, flambs, obj_idx):
    ''' 
    INPUTS:
        - f: filters
        - er: errors as given from catalogs
        - z: redshift
        - flambds: central wavelengths for the filters
        - obj_idx: index, in superFluxes, of current object
        
    for each band, the error for that band is the maximum of err_ratio
    (a percent of the flux of that band), the actual error ratio, and in some
    cases half the contribution of an emission line if an emission line is in
    that band.this function handles OII, H_beta, H_alpha, and OIII respectively
    
    in addition, this function makes sure an object has atleast 3 observations
    less than 10000, 2 between 10000 and 25000, and one greater than 30000 
    (all observed frame, angstroms) 
    
    another point: we throw out any observation with a S/N < 3
    '''
    # changeable thresholds
    err_ratio = 5.0 #percent
    sig_noise_ratio = 3.0
    
    # define some tools
    nFilt = len(f)
    outEr = np.zeros(nFilt)
    outF = np.zeros(nFilt)
    
    # define minimum required counts for acceptance of an object
    uv_max = 3
    opt_max = 2
    ir_max = 1
    
    # define current counts for uv, optical, and ir
    uv_count = 0
    opt_count = 0
    ir_count = 0
    
    # emission wavelenghts
    emis_lambdas = np.array([3725, 4861, 6563, 5007]) #o2, hb, ha, o3
    emis_coeff = np.array([3.0, 4.0, 1.5, 1.0]) #o2, hb, ha, o3, 1 for o3
    #get from henry's list for this object in cgs
    o3flux = np.loadtxt(home+'data/z2/data_house/flags.dat')[obj_idx,23]*1e-17
    
    # loop through each filter
    for fi in range(0, nFilt):
        
        # make sure there is an observation
        if f[fi] == 0:
            continue
        
        # make sure pass S/N test
        if f[fi]/er[fi] < sig_noise_ratio:
            continue
        
        current_lambda = flambs[fi]
        
        # make sure this filter is not centered below 1250 rest frame
        if current_lambda/(1.0 + z) < 1250:
            continue
    
        # add a count
        if current_lambda <= 10000: #uv, ang
            uv_count = uv_count + 1
        if current_lambda > 10000 and current_lambda <= 25000: #optical, ang
            opt_count = opt_count + 1
        if current_lambda >= 30000: #ir, ang
            ir_count = ir_count + 1
            
        # see if OII, H_beta, H_alpha, OIII are in band, in that order
        # if so, calculate emissionLineContribution
        emissionLineContribution = 0 # default
        # get filter curve
        filter_data = np.loadtxt(home+'data/filters/filter_house/'+col_names[fi]+'.res')
        lambdas = filter_data[:,0] #units= Ang
        res_o = filter_data[:,1] #o standing for origional
        res = res_o/np.max(res_o)
        
        # get all indexes where res is greater than 10%
        applicableIdxs = np.where(res >= 0.1)[0]
        llo = lambdas[applicableIdxs[0]] # since res is ordered can do this
        lhi = lambdas[applicableIdxs[-1]]# since res is ordered can do this
        lrange = lhi - llo
        avgF_l_o3 = o3flux/lrange #in F/Ang
        
        # check if each emission line is wihin lrange
        for ei in range(0, len(emis_lambdas)):
            emisWaveObs = emis_lambdas[ei]*(1.0+z)
            if llo < emisWaveObs and lhi > emisWaveObs: #within lrange defination
                avgF_nu_o3 = (avgF_l_o3*current_lambda**2)/3e18 #F_nu = (F_l * L^2) / c
                avgF_nu_currentLine = avgF_nu_o3 / emis_coeff[ei]
                emissionLineContribution = avgF_nu_currentLine/f[fi] #ratio = line component/entire band
                emissionLineContribution = emissionLineContribution/1e29 # not sure why
                
        # decide on an error
        outEr[fi] = np.max( [f[fi]*(err_ratio/100.0),    # error floor in %
                            er[fi],                          # actual error
                            emissionLineContribution]) # emission line part
        
        # put the flux into proper spot if gotten this far in for loop
        outF[fi] = f[fi]


    # make sure pass count test before returning
    if uv_count < uv_max:
        return np.zeros(nFilt), np.zeros(nFilt)
    if opt_count < opt_max:
        return np.zeros(nFilt), np.zeros(nFilt)
    if ir_count < ir_max:
        return np.zeros(nFilt), np.zeros(nFilt)
        
    # if we've made it this far, we're good
    return outF, outEr
        
    