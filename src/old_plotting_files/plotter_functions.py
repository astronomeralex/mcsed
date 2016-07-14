# -*- coding: utf-8 -*-
# =========================================================================== #
''' ======================================================================= '''
# =========================================================================== #
import numpy as np

# =========================================================================== #
''' ======================================================================= '''
# =========================================================================== #
def plotTriangle(new_flatchain, field, objid, dataIDX, save_opt=None):
    import triangle
    import matplotlib.pylab as plt
    vNames = ['tau (dust2)','delta (dust_index)', 'Eb (uvb)', 'log10(age)', 'log10(M_star)']
    triangle_figure = triangle.corner(new_flatchain, labels=vNames,
                                      quantiles=[0.16, 0.5, 0.84],
                                      show_titles=True, title_args={"fontsize": 12},
                                      plot_contours=True,
                                      plot_datapoints=False,
                                      verbose=False
                                      )
                                    
    
    triangle_figure.gca().annotate(field+'-'+objid+', hunters index: '+str(dataIDX), 
                        xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")
    
                      
    if save_opt != None:
        triangle_figure.savefig(save_opt+'triangle.png')
        plt.close('all')
    else:
        triangle_figure.show()
        
    
        
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================
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
        #show entire range of possible parameter values
        plt.ylim([np.min(vHouse.paramRanges[param_i,:]),
                  np.max(vHouse.paramRanges[param_i,:])]) 
        if save_opt != None:
            plt.savefig(save_opt+'param'+str(param_i)+'path.png')
            plt.close('all')
        else:
            plt.show()
            
 # =========================================================================== #
''' ======================================================================= '''
# =============================================================================
def plotLnprobPath(lnprobability, OBJi, save_opt=None):
    nSteps = np.shape(lnprobability)[1]
    nWalkers = np.shape(lnprobability)[0]
    import matplotlib.pylab as plt
    stepidx = np.arange(1,nSteps+1)
    #plots lnprob for each walker over step idx
    for walkeri in range(0,nWalkers):
        plt.plot(stepidx, lnprobability[walkeri,:], alpha=0.1, color='k')
    plt.xlabel('STEP IDX')
    plt.ylabel('lnprob()')
    if save_opt != None:
        plt.savefig(save_opt+'lnprobpath.png')
        plt.close('all')
    else:
        plt.show()
        
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def plotSEDandResids(aInput,
                     filterLambdas, dataAvgFlux, dataAvgFluxErr, nu_obs,
                     z, OBJi, field, objid,
                     save_opt=None):
                         
    ''' Observed Frame. Lambda vs. Flux_nu. '''
    ax1_legendNames=[]
    ax2_legendNames=[]
    import variable_house as vHouse
    paramNames = vHouse.nameL
    nDim = vHouse.nDim
    import math_functions as mathF
    import matplotlib.pylab as plt
    from matplotlib.ticker import FormatStrFormatter
    plt.figure(figsize=(8*3, 6*3), dpi=80)
    ax1 = plt.subplot(211)
    
    
    # cut indicies around only where we want
    lambda_obs = mathF.c_ang/nu_obs
    a = np.abs(lambda_obs-np.min(filterLambdas))
    b = np.abs(lambda_obs-np.max(filterLambdas))
    loi = int(np.where(a==np.min(a))[0])
    hii = int(np.where(b==np.min(b))[0])
    
    
    # make these tick labels invisible
    plt.setp( ax1.get_xticklabels(), visible=False)
             
    ## share x only
    ax2 = plt.subplot(212, sharex=ax1)
    
    
     # real values
    dataAvgFlux_uJy = dataAvgFlux * 10**(29) # put int uJy
    dataAvgFluxErr_uJy = dataAvgFluxErr * 10**(29) #put into uJy
    ax1.errorbar(filterLambdas, dataAvgFlux_uJy, 
             yerr=dataAvgFluxErr_uJy, 
             fmt='o', color='k', ecolor='k', capsize=0)
    ax1_legendNames.append('Observations')
             
    
    
    aAvgFlux_uJy_m = np.zeros([np.shape(aInput)[0],len(dataAvgFlux)])#matrix meaning holds all
    for ai in range(0,np.shape(aInput)[0]):
        # first SED
        nu_obs, a_s_nu, aAvgFlux = arbitraryTruths(aInput[ai,:], z, dataAvgFlux, 
                                                   dataAvgFluxErr, paramNames)    
        a_s_nu_uJy = 10**(29) * a_s_nu #into uJy
        aAvgFlux_uJy_m[ai,:] = 10**(29) * aAvgFlux #into uJy
        #lambda_obs = mathF.c_ang/nu_obs
        ax1.plot(lambda_obs[loi:hii], a_s_nu_uJy[loi:hii])
        
        # first residuals
        ok_IDXs = np.where(dataAvgFluxErr != 0)[0]
        resids = np.zeros(len(ok_IDXs))
        ok_filterLambdas = np.zeros(len(ok_IDXs))
        ok_aAvgFlux_uJy = np.zeros(len(ok_IDXs))
        
        for ia in range(0, len(ok_IDXs)):
            okIDX = ok_IDXs[ia]
            resids[ia] = (dataAvgFlux_uJy[okIDX] - aAvgFlux_uJy_m[ai,:][okIDX]) / dataAvgFluxErr_uJy[okIDX]
            ok_filterLambdas[ia] = filterLambdas[okIDX]
            ok_aAvgFlux_uJy[ia] = aAvgFlux_uJy_m[ai,:][okIDX]
            
        X2 = getChiSquared(aAvgFlux_uJy_m[ai,:], dataAvgFlux_uJy, dataAvgFluxErr_uJy)
        rX2 = X2/(len(ok_IDXs)-nDim+1)
        rX2s = str(rX2)[:4]
        ax2.scatter(ok_filterLambdas, resids)
        ax2_legendNames.append('ChiSq/Nu_DoF= '+str(rX2s))
        
        ax1.scatter(ok_filterLambdas, ok_aAvgFlux_uJy, marker='s', s=75)
        
    
    
    # plot resid = 0, and plus/minus 1 error
    ax2.plot(ok_filterLambdas, np.zeros(len(ok_filterLambdas)), 'k')
    ax2.plot(ok_filterLambdas, np.ones(len(ok_filterLambdas)), 'k', alpha=0.25)
    ax2.plot(ok_filterLambdas, -(np.ones(len(ok_filterLambdas))), 'k', alpha=0.25)
    
    # make graph nice
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_yscale('log')
    plt.xlabel('$\lambda$, Observed Frame   [ $\AA$ ]')
    ax1.set_ylabel(r'log( $F_{\nu}$ )   [ uJy ]')
    ax2.set_ylabel('(Observation-Model)/Error')
    ax1.set_title(field+'-'+objid+', Index: '+str(OBJi))
    ax1.legend(ax1_legendNames, loc=0)
    ax2.legend(ax2_legendNames, loc=0)
    plt.subplots_adjust(hspace=0.0)
    # auto size plot
    minx = np.min(ok_filterLambdas)
    maxx = np.max(ok_filterLambdas)
    ax1.set_xlim([minx-0.10*minx, maxx+0.10*maxx])
    
    ax1.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    ax1.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    #ax1.set_xticks(np.arange(minx,maxx,1000))
    
    if save_opt != None:
        plt.savefig(save_opt)
    else:
        plt.show()
        
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def arbitraryTruths(aInput, z, dataAvgFlux, dataAvgFluxErr, paramNames):
    import variable_house as vHouse
    import math_functions as mathF
    import factory
    import wheres_home
    home = wheres_home.getHomeLocation()
    # plot arbitrary data
    sp = vHouse.defSP()
    sp, tage = vHouse.updateFSPS(sp , vHouse.nameL, aInput)
    waves, lum_e = sp.get_spectrum(tage=tage)
    a_s_nu, nu_obs = mathF.getFlux(lum_e, waves, z)
    R_nu, nFilters = factory.R_nu(home, nu_obs) #requires that sample filters are in current home                     
    D_nu = factory.D_nu(nFilters, nu_obs)
    NU_obs = factory.NU_obs(nFilters,nu_obs)
    S_nu = factory.S_nu(nFilters,a_s_nu)           
    aAvgFlux = mathF.getAvgFlux(   nFilters,
                                      S_nu,
                                      NU_obs,
                                      D_nu,
                                      R_nu,
                                      True)
    aA = mathF.getA(aAvgFlux, dataAvgFlux, dataAvgFluxErr)
    aAvgFlux = aA*aAvgFlux
    a_s_nu = aA*a_s_nu
    return nu_obs, a_s_nu, aAvgFlux
    
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def getChiSquared(model, data, dataerr):
    import math_functions as mathF
    import mcsed_utils
    A =  mathF.getA(model, data, dataerr)
    top = (A*model-data)**2
    bottom = dataerr**2
    print 'THIS MIGHT BE BROKEN'
    return mcsed_utils.fracsum(top , bottom)
    
 # =========================================================================== #
''' ======================================================================= '''   
# =============================================================================    
def plotSelfAny2(predictions, errors, p1_idx, p2_idx, pColor_idx=None, KC13=False):
    
    import matplotlib.pylab as plt
    import variable_house as vHouse
    import matplotlib.lines as mlines
    
    cm = plt.cm.get_cmap('Reds')
    p1 = []
    p1_l_err = []
    p1_u_err = []
    p2 = []
    p2_l_err = []
    p2_u_err = []
    coloring = []
    
    nonzero = np.where(predictions[:,0] != 0)[0]

    for i in range(0, len(nonzero)):
        idx = nonzero[i]
        p1.append(float(predictions[idx,p1_idx]))
        p1_l_err.append(float(errors[idx,p1_idx,0]))
        p1_u_err.append(float(errors[idx,p1_idx,1]))
        p2.append(float(predictions[idx,p2_idx]))
        p2_l_err.append(float(errors[idx,p2_idx,0]))
        p2_u_err.append(float(errors[idx,p2_idx,1]))
        coloring.append(float(predictions[idx, pColor_idx]))
        
# because errors are a percentile, subtract to median from them    
    p1_l_err = np.asarray(p1)-np.asarray(p1_l_err)
    p2_l_err = np.asarray(p2)-np.asarray(p2_l_err)
    p1_u_err = np.asarray(p1_u_err)-np.asarray(p1)
    p2_u_err = np.asarray(p2_u_err)-np.asarray(p2)
    
    if pColor_idx != None:
    
        sc = plt.scatter(p1, p2, 
                         c=coloring, vmin=np.min(coloring), vmax=np.max(coloring), 
                         s=65, cmap=cm, label='50% of sample', alpha = 0.5)
        plt.errorbar(p1, p2, xerr=(p1_l_err, p1_u_err), yerr=(p2_l_err, p2_u_err), 
                     linestyle="None", color='k', capsize=0, barsabove=False, alpha=0.15)   
                     
                     
    if pColor_idx == None:
        sc = plt.scatter(p1, p2, 
                         s=65, cmap=cm, label='50% of sample', alpha = 0.5)
        plt.errorbar(p1, p2, xerr=(p1_l_err, p1_u_err), yerr=(p2_l_err, p2_u_err), 
                     linestyle="None", color='k', capsize=0, barsabove=False, alpha=0.15)   
                     
                 
    if KC13 == True:
        #Eb =(0.85±0.09)−(1.9±0.4)δ
        delta = np.linspace(-0.6,0.3,2)
        eb = 0.85-1.9*delta
        #kc13 = plt.plot(delta, eb, 'k', label='KC13')
        plt.plot(delta, eb, 'k', label='KC13')
        kc13 = mlines.Line2D([np.min(delta),np.min(eb)], [np.max(delta),np.max(eb)], 
                              color='k', marker='*',
                          markersize=0, label='Blue stars')
        
        
        delta = np.linspace(-1.0,1.0,2)
        eb = 0.85-1.9*delta
        kc13_e = plt.plot(delta, eb, 'k--',label='KC13 Extrapolated')

        
        a = np.array([[0,0]])
        cal = plt.scatter(a[:,0], a[:,1], color='b', marker='^', s=100, label='Calzetti')
        
        a = np.array([[-0.1,3]])
        lmc = plt.scatter(a[:,0], a[:,1], color='b', marker='+', s=100, label='Calzetti')
        
        a = np.array([[-0.45,0]])
        smc = plt.scatter(a[:,0], a[:,1], color='b', marker='x', s=100, label='Calzetti')
        
        
        
        a = np.array([[0,4.1]])
        mw = plt.scatter(a[:,0], a[:,1], color='b', marker='d', s=100, label='Calzetti')
        
        a = np.array([[-0.2,0]])
        reddy = plt.scatter(a[:,0], a[:,1], color='b', marker='s', s=100, label='Calzetti')
        
        a = np.array([[0.1,1]])
        sco = plt.scatter(a[:,0], a[:,1], color='b', marker='*', s=100, label='Calzetti')
        
        a = np.array([[0.46,-0.21]])
        gregL = plt.scatter(a[:,0], a[:,1], color='y', marker='o', s=100, label='Calzetti')
        plt.errorbar(a[:,0], a[:,1], xerr=[0.04], yerr=[0.17],color='y', 
                 linestyle="None",  capsize=0, barsabove=False)   
        
        a = np.array([[0.31,-0.15]])
        gregM = plt.scatter(a[:,0], a[:,1], color='cyan', marker='o', s=100, label='Calzetti')
        plt.errorbar(a[:,0], a[:,1], xerr=[0.039], yerr=[0.19],color='cyan', 
                 linestyle="None",  capsize=0, barsabove=False)   
        
        a = np.array([[0.14,-0.09]])
        gregH = plt.scatter(a[:,0], a[:,1],color='purple', marker='o', s=100, label='Calzetti')
        plt.errorbar(a[:,0], a[:,1], xerr=[0.049], yerr=[0.23],color='purple', 
                 linestyle="None",  capsize=0, barsabove=False)
        
    plt.title(r'$\delta$ vs $E_b$, colored by log10(Stellar Mass)')
    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$E_b$')
    plt.colorbar(sc)
    plt.legend([sc,kc13,cal, lmc, smc, mw, reddy, sco, gregL, gregM, gregH],
               ['50% of sample', 'KC13', 'Calzetti', 'LMC', 'SMC','MW', 'Reddy', 'Scoville', '7.2<logM*<8.1', '8.1<logM*<8.8', '8.8<logM*<10.2'],
                #handles = [kc13],
                loc=0)
    plt.show()
    
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def compareMasses(his, hiserr, mine, myerr, hisfield, hisid, myfield, myid, z):   
    import matplotlib.pylab as plt
         
    a = []
    b = []
    c = []
    d = []
    e = []
    
    for i in range(0, np.shape(mine)[0]):
        if mine[i] != 0:
            #his_id_idx = np.where(hisid == myid[i])[0]
            pos_id_idx = np.where(hisid == myid[i])[0]
            if len(pos_id_idx) == 1:
                his_id_idx = int(pos_id_idx)
            else:
                pos_field_idx = np.where(hisfield == myfield[i])[0]
                for ii in range(0, len(pos_field_idx)):
                    for iii in range(0, len(pos_id_idx)):
                        if pos_field_idx[ii] == pos_id_idx[iii]:
                            his_id_idx = int(pos_id_idx[iii])
    
            #print str(myid[i]), str(myfield[i]), str(hisid[his_id_idx]), str(hisfield[his_id_idx])
            assert myfield[i] == hisfield[his_id_idx]
            assert myid[i] == hisid[his_id_idx]
            if his[his_id_idx] != -999:
                a.append(his[his_id_idx])
                b.append(hiserr[his_id_idx])
                c.append(mine[i]-np.log10(((1.0+z[i]))))#-np.log10(2.75))#(1.0+z[i]))
                d.append(myerr[i,0]-np.log10(((1.0+z[i]))))#-np.log10(2.75))#1.0+z[i]))
                e.append(myerr[i,1]-np.log10(1.0+z[i]))#-np.log10(2.75))#1.0+z[i]))
                
    # since errors are just percentile, need difference from median
    d = np.asarray(c)-np.asarray(d)
    e = np.asarray(e)-np.asarray(c)
    
    plt.errorbar(a, c, xerr=b, yerr=[d,e], fmt='o', 
                 color='b', capsize=0, alpha=0.50)
    
    # plot the y=x line
    x = np.linspace(np.min((a,c)),
                    np.max((a,c)),
                    10)
    plt.plot(x, x, 'k--')
    
    '''
    # plot the best fit line
    m, b = np.polyfit(a, c, 1)
    plt.plot(x, m*x + b, 'r--')
    '''
    
    plt.xlabel("Alex's log10(mass)")
    plt.ylabel('MCSED log10(mass) / (1+z)')
    plt.title("Alex's Masses vs MCSED Masses/(1+z)")#, Red=LineOfBestFit(m: "+str(m)+", b:"+str(b)+"), Black=OneToOneLine")
    #plt.legend(['With Neb. Emis.'],loc=0)#, 'W/o Neb. Emis'], loc=0)
        
    plt.show()
    
    
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def compareMassRatioVsRedshift(his, hiserr, mine, myerr, hisfield, hisid, myfield, myid, z):   
    import matplotlib.pylab as plt
         
    a = []
    b = []
    c = []
    d = []
    e = []
    z_list = []
    
    for i in range(0, np.shape(mine)[0]):
        if mine[i] != 0:
            #his_id_idx = np.where(hisid == myid[i])[0]
            pos_id_idx = np.where(hisid == myid[i])[0]
            if len(pos_id_idx) == 1:
                his_id_idx = int(pos_id_idx)
            else:
                pos_field_idx = np.where(hisfield == myfield[i])[0]
                for ii in range(0, len(pos_field_idx)):
                    for iii in range(0, len(pos_id_idx)):
                        if pos_field_idx[ii] == pos_id_idx[iii]:
                            his_id_idx = int(pos_id_idx[iii])
    
            #print str(myid[i]), str(myfield[i]), str(hisid[his_id_idx]), str(hisfield[his_id_idx])
            assert myfield[i] == hisfield[his_id_idx]
            assert myid[i] == hisid[his_id_idx]
            if his[his_id_idx] != -999:
                a.append(his[his_id_idx])
                b.append(hiserr[his_id_idx])
                c.append(mine[i])#-np.log10(((1.0+z[i]))))
                d.append(myerr[i,0])#-np.log10(((1.0+z[i]))))
                e.append(myerr[i,1])#-np.log10(1.0+z[i]))
                z_list.append(z[i])
                
    # since errors are just percentile, need difference from median
    d = np.asarray(c)-np.asarray(d)
    e = np.asarray(e)-np.asarray(c)
    
    c_a = 10**np.array(c)
    a_a = 10**np.array(a)
    
    ratio = c_a/a_a
    
    plt.errorbar(z_list, ratio, fmt='o', 
                 color='b', capsize=0, alpha=0.50)
                 
    # plot the y=x line
    x = np.linspace(np.min(z_list),
                    np.max(z_list),
                    10)
    plt.plot(x, np.ones(len(x)), 'k--')
                 
    plt.yscale('log')
    
    plt.xlabel("Redshift")
    plt.ylabel("MCSED  / Alex's  [note: ratio of actual masses]")
    plt.title("Redshift vs Mass Ratio (MCSED/Alex)")
    #plt.legend(['With Neb. Emis.'],loc=0)#, 'W/o Neb. Emis'], loc=0)
        
    plt.show()
    
# =========================================================================== #
''' ======================================================================= '''   
# =============================================================================        
def gridSmoothness(nRandSamples, p1_name, p1_idx, p2_name, p2_idx):
    import wheres_home
    home = wheres_home.getHomeLocation()
    import math_functions as mathF
    import factory
    import scipy.interpolate as spterp
    import variable_house as vHouse
    import matplotlib.pylab as plt
    
    # ======================================================================= #
    assert vHouse.walker_SED_type == 'Interp', "This script looks at smoothness of the pSpace!"
    
    pSpaceGridSize = vHouse.pSpaceGridSize
    npzopendata = np.load(home+'data/pSpaces/size'+str(pSpaceGridSize)+'_sed.npz')
    sedGrid = npzopendata["sedGrid"]
    npzopendata.close()
    npzopendata = np.load(home+'data/pSpaces/size'+str(pSpaceGridSize)+'_other.npz')
    #valsGrid=npzopendata['valsGrid']  #BROKEN, see create_pSpace.py
    lambda_e=npzopendata['waves']
    #dustMassGrid=npzopendata['dustMassGrid']
    #stellarMassGrid=npzopendata['stellarMassGrid']
    npzopendata.close()
    
    # ======================================================================= #
    # get middle index of pSpaceGridSize, higher one if no exact middle
    paramIndicies = np.arange(0,pSpaceGridSize)
    try:
        middle_idx = np.median(paramIndicies)
        extra = middle_idx%1.0
        assert middle_idx%1.0 == 0.0
    except AssertionError:
        middle_idx = int(middle_idx + extra)
    
    # ======================================================================= #
    # enter a : by hand at the values of p1_idx and p2_idx
   # pSpace = sedGrid[:,middle_idx,middle_idx,:]
    chi_sq_grid = np.zeros([pSpaceGridSize,pSpaceGridSize])
    chi_sq_rand = np.zeros([nRandSamples])
    rand_ax1 = np.random.uniform(low=vHouse.lowerL[p1_idx],
                                 high=vHouse.higherL[p1_idx],
                                 size=nRandSamples)
    rand_ax2 = np.random.uniform(low=vHouse.lowerL[p2_idx],
                                high=vHouse.higherL[p2_idx],
                                size=nRandSamples)
                                
    # ======================================================================= #
    # load data for an object
    col_names = np.loadtxt(home+'data/active/filter_names.dat', dtype='str')
    nFilters = len(col_names)
    superFluxes = np.loadtxt(home+'data/active/superFluxes.dat')
    superErrors = np.loadtxt(home+'data/active/superErrors.dat')
    #superInfo = np.loadtxt(home+'data/active/superInfo.dat')
    superRedshift = np.loadtxt(home+'data/active/superRedshift.dat')
    nObjects, nSuperCols = np.shape(superFluxes)
    
    #for run_i in range(0, nObjects):
    for run_i in range(70, 71):
        
        dataAvgFlux = np.zeros([nFilters])
        dataAvgFluxErr = np.zeros([nFilters])
        
        z = superRedshift[run_i] #redshift for this object
        assert z >= 0.0, "Redshift must be positive...Right? Math breaks when <0?"
        
        dataAvgFlux = np.zeros([nFilters])
        dataAvgFluxErr = np.zeros([nFilters])
        
        # part of renaming input data...
        for filt_i in range(0, nFilters):
            # pull data from superFluxes/Errors into working data
            dataAvgFlux[filt_i] = superFluxes[run_i, filt_i  ]
            dataAvgFluxErr[filt_i] = superErrors[run_i, filt_i  ]
        
        # run the errors through error machine
        dataAvgFluxErr = error_machine(dataAvgFlux, dataAvgFluxErr)
                      
        # finally, make sure can we can actually do a walk, we want 2+ observations
        try:
            assert len(np.where(dataAvgFlux != 0.0)[0]) > 2 # more than 2 data points
        except AssertionError:
            print 'There is not enough data for this object. Skipping ',run_i
            continue    
        
    # ======================================================================= #
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
        p0 = np.zeros([vHouse.nWalkers, vHouse.nDim])
        for pi in range(0,vHouse.nDim):
            p0[:, pi] = np.random.uniform(vHouse.lowerL[pi],
                                          vHouse.higherL[pi],
                                          [vHouse.nWalkers])
        other1val = vHouse.paramRanges[1,middle_idx]
        other2val = vHouse.paramRanges[2,middle_idx]
        for a in range(0, pSpaceGridSize):
            for b in range(0, pSpaceGridSize):
                p1val = vHouse.paramRanges[p1_idx,a]
                p2val = vHouse.paramRanges[p2_idx,b]
                x = [p1val, other1val, other2val, p2val]
                chi_sq_grid[a,b] = lnprob(x,flux_interp_f,lambda_e,nFilters,
                                           z,NU_obs,D_nu,R_nu,dataAvgFlux,dataAvgFluxErr)                          
        
        for a in range(0, np.shape(rand_ax1)[0]):
            p1val = rand_ax1[a]
            p2val = rand_ax2[b]
            x = [p1val, other1val, other2val, p2val]
            chi_sq_rand[a] = lnprob(x,flux_interp_f,lambda_e,nFilters,
                                       z,NU_obs,D_nu,R_nu,dataAvgFlux,dataAvgFluxErr)
                
    # ======================================================================= #
        # plot
        cm = plt.cm.get_cmap('Reds')
        p1_grid, p2_grid = np.meshgrid(vHouse.paramRanges[p1_idx,:],
                                       vHouse.paramRanges[p2_idx,:])
        sc2 = plt.scatter(rand_ax1, rand_ax2,
                         c=chi_sq_rand, 
                         vmin=np.min(chi_sq_rand), vmax=np.max(chi_sq_rand), 
                         s=25, cmap=cm)       
        sc = plt.scatter(p1_grid, p2_grid,
                         c=chi_sq_grid, 
                         vmin=np.min(chi_sq_grid), vmax=np.max(chi_sq_grid), 
                         s=75, cmap=cm)
        plt.title(p1_name+' vs '+p2_name+' colored by Chi Squared')
        plt.xlabel(p1_name)
        plt.ylabel(p2_name)
        plt.colorbar(sc)
        plt.colorbar(sc2)
        plt.show()
    
    
# ======================================================================== #
# define a function to handle observational errors
def error_machine(obs_in, errs_in):
    ''' This function only permits observations with signal/noise > 3.
    Additionally, it adds a 5% error floor in quadrature.  If an error is 
    more than 5% of the observed flux, it is kept.'''
    output_errs = np.zeros(np.shape(errs_in))
    for filt_i in range(0, len(obs_in)):
        # signal/noise > 3:
        if obs_in[filt_i]/errs_in[filt_i] >= 3.0:
            if obs_in[filt_i]*0.05 > errs_in[filt_i]:
                output_errs[filt_i] = obs_in[filt_i]*0.05
            else:
                output_errs[filt_i] = errs_in[filt_i]   
        else:
            pass
    return output_errs
    
# =========================================================================== #
# check to make sure parameters are in bounds
def inBoundsCheck(x):
    import variable_house as vHouse
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
            dataAvgFluxErr):
                
    import variable_house as vHouse
    import factory
    import mcsed_utils as myTools
    import math_functions as mathF
     
    # check to make sure parameters are in bounds
    if inBoundsCheck(x) == False:
        return -1.0*np.inf
    else:
        pass
    if vHouse.walker_SED_type == 'Direct':
        assert False==True,  "Don't use this yet."
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
    # return the lnprob
    chi_sq = myTools.fracsum((dataAvgFlux-modelAvgFlux)**2,(dataAvgFluxErr)**2)
    return (-1.0/2.0)*chi_sq
    
# =========================================================================== #
''' ======================================================================= '''
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
''' ======================================================================= '''
# =========================================================================== #
def testInterp( plow, phigh, ntests, z,
               
                nFilters,
                NU_obs,
                D_nu, 
                R_nu,
               
                dataAvgFlux,
                dataAvgFluxErr):
                    
    pvals = np.linspace(plow, phigh, ntests)
    lnprobvals = np.zeros(ntests)
    seds = np.zeros([ntests, 1963])
    nu_obs_m = np.zeros([ntests, 1963])
    as_l = np.zeros(ntests)
    
    import variable_house as vHouse
    sp = vHouse.defSP()
    import math_functions as mathF
    import factory
    import mcsed_utils as myTools
    
    
    for i in range(0, ntests):
        all_pvals = [0.36, 0.74, 0.35, pvals[i]]
        sp,age = vHouse.updateFSPS(sp,vHouse.nameL,all_pvals)
        waves, lum_e = sp.get_spectrum(tage=age)
        s_nu, nu_obs = mathF.getFlux(lum_e, waves, z)

        model_S_nu = factory.S_nu(nFilters, s_nu)
        modelAvgFlux = mathF.getAvgFlux(     nFilters,
                                             model_S_nu,
                                             NU_obs,
                                             D_nu,
                                             R_nu,
                                             True)
        A = mathF.getA(modelAvgFlux, dataAvgFlux, dataAvgFluxErr)
        # return the lnprob
        chi_sq = myTools.fracsum((dataAvgFlux-A*modelAvgFlux)**2,(dataAvgFluxErr)**2)
        lnprobvals[i] = (-1.0/2.0)*chi_sq
        as_l[i] = A
        seds[i,:] = A*s_nu*10**(29) #account for uJy
        nu_obs_m[i,:] = nu_obs
    
    
    import matplotlib.pylab as plt
    #plt.plot(nu_obs_m[:,:100],seds[:,:100])
    plt.scatter(pvals, as_l, color='r', marker='^')

    plt.xlabel('log10(tage)')
    plt.ylabel('A')
    #plt.xlabel('wavelength (observed)')
    #plt.ylabel('uJy flux density per Hz')
    #plt.title('34-4770 (77) w/ [0.36, 0.74, 0.35, tage]')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.show()
    
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def histogram(predictions, pidx, title):
    
    okidx = np.where(predictions[:,pidx] != 0)[0]
    data = np.zeros(len(okidx))

    for i in range(0, len(okidx)):
        idx = okidx[i]
        data[i] = predictions[ idx, pidx ]
        
    import matplotlib.pylab as plt
    plt.hist(data,bins=40)
    
    print np.mean(data)
    print np.median(data)
    
    plt.title(r'Histogram of '+title+' , N = '+str(len(data)))
    plt.xlabel('parameter value')
    plt.ylabel('frequency in bin')
    
    plt.show()
    
# =========================================================================== #
''' ======================================================================= '''
# =============================================================================    
def plotSEDandResidsREST(aInput,
                         filterLambdas, dataAvgFlux, dataAvgFluxErr, nu_obs,
                         z, OBJi, field, objid,
                         save_opt=None):
                         
    ''' Rest Frame. Lambda vs. Flux_nu. '''

    ax1_legendNames=[]
    ax2_legendNames=[]
    
    import variable_house as vHouse
    paramNames = vHouse.nameL
    nDim = vHouse.nDim
    
    import math_functions as mathF
    import matplotlib.pylab as plt
    
    from matplotlib.ticker import FormatStrFormatter
    plt.figure(figsize=(8*3, 6*3), dpi=80)
    ax1 = plt.subplot(211)
    
    # cut indicies around only where we want
    lambda_obs = mathF.c_ang/nu_obs
    a = np.abs(lambda_obs-np.min(filterLambdas))
    b = np.abs(lambda_obs-np.max(filterLambdas))
    loi = int(np.where(a==np.min(a))[0])
    hii = int(np.where(b==np.min(b))[0])
    
    
    # make these tick labels invisible
    plt.setp( ax1.get_xticklabels(), visible=False)
             
    ## share x only
    ax2 = plt.subplot(212, sharex=ax1)
    
    
     # real values
    dataAvgFlux_uJy = dataAvgFlux * 10**(29) # put int uJy
    dataAvgFluxErr_uJy = dataAvgFluxErr * 10**(29) #put into uJy
    rest_filterLambdas = filterLambdas/(1.0+z) #unredden observationas
    ax1.errorbar(rest_filterLambdas, dataAvgFlux_uJy, 
             yerr=dataAvgFluxErr_uJy, 
             fmt='o', color='k', ecolor='k', capsize=0)
    ax1_legendNames.append('Observations')
             
    
    aAvgFlux_uJy_m = np.zeros([np.shape(aInput)[0],len(dataAvgFlux)]) #matrix meaning holds all
    for ai in range(0,np.shape(aInput)[0]):
        # first SED
        nu_obs, a_s_nu, aAvgFlux = arbitraryTruths(aInput[ai,:], z, dataAvgFlux, 
                                                   dataAvgFluxErr, paramNames)    
        a_s_nu_uJy = 10**(29) * a_s_nu #into uJy
        aAvgFlux_uJy_m[ai,:] = 10**(29) * aAvgFlux #into uJy
        lambda_rest = lambda_obs/(1.0+z)
        ax1.plot(lambda_rest[loi:hii], a_s_nu_uJy[loi:hii])
        
        # first residuals
        ok_IDXs = np.where(dataAvgFluxErr != 0)[0]
        resids = np.zeros(len(ok_IDXs))
        ok_filterLambdas = np.zeros(len(ok_IDXs))
        ok_aAvgFlux_uJy = np.zeros(len(ok_IDXs))
        
        for ia in range(0, len(ok_IDXs)):
            okIDX = ok_IDXs[ia]
            resids[ia] = (dataAvgFlux_uJy[okIDX] - aAvgFlux_uJy_m[ai,:][okIDX]) / dataAvgFluxErr_uJy[okIDX]
            ok_filterLambdas[ia] = filterLambdas[okIDX]
            ok_aAvgFlux_uJy[ia] = aAvgFlux_uJy_m[ai,:][okIDX]
            
        X2 = getChiSquared(aAvgFlux_uJy_m[ai,:], dataAvgFlux_uJy, dataAvgFluxErr_uJy)
        rX2 = X2/(len(ok_IDXs)-nDim+1)
        rX2s = str(rX2)[:4]
        ax2.scatter(ok_filterLambdas, resids)
        ax1_legendNames.append('ChiSq/Nu_DoF= '+str(rX2s))
        
        ax1.scatter(ok_filterLambdas, ok_aAvgFlux_uJy, marker='s', s=75)
        
    
    
    # plot resid = 0, and plus/minus 1 error
    ax2.plot(ok_filterLambdas, np.zeros(len(ok_filterLambdas)), 'k')
    ax2.plot(ok_filterLambdas, np.ones(len(ok_filterLambdas)), 'k', alpha=0.25)
    ax2.plot(ok_filterLambdas, -(np.ones(len(ok_filterLambdas))), 'k', alpha=0.25)
    
    # make graph nice
    ax1.set_xscale('log')
    #ax2.set_xscale('log')
    ax1.set_yscale('log')
    plt.xlabel('$\lambda$, Observed Frame   [ $\AA$ ]')
    ax1.set_ylabel(r'log( $F_{\nu}$ )   [ uJy ]')
    #ax2.set_ylabel('(Data - Model) / (Error)')
    ax1.set_title(field+'-'+objid+', Index: '+str(OBJi))
    ax1.legend(ax1_legendNames, loc=0)
    #ax2.legend(ax2_legendNames, loc=0)
    plt.subplots_adjust(hspace=0.0)
    # auto size plot
    minx = np.min(ok_filterLambdas)
    maxx = np.max(ok_filterLambdas)
    ax1.set_xlim([minx-0.10*minx, maxx+0.10*maxx])
    
    ax1.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    ax1.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    #ax1.set_xticks(np.arange(minx,maxx,1000))
    
    if save_opt != None:
        plt.savefig(save_opt)
    else:
        plt.show()
        
# =========================================================================== #
''' ======================================================================= '''   
# =============================================================================    
def plotGlobalAny2(p1_data, 
                   p1_errors, 
                   p1_idx, 
                   
                   p2_data,
                   p2_errors,
                   p2_idx,
                   
                   title,
                   xlabel,
                   ylabel,
                   
                   pColor_idx=None,
                   p1_mine=False,
                   p2_mine=False):
    
    import matplotlib.pylab as plt
    
    if pColor_idx != None:
        cm = plt.cm.get_cmap('Reds')
        coloring = []
        
    p1 = []
    p1_l_e = []
    p1_u_e = []
    p1_e = []
    
    p2 = []
    p2_l_e = []
    p2_u_e = []
    p2_e = []
    
    
    

    nonzero = np.where(p2_data[:,3] != 0)[0]
        

    for i in range(0, len(nonzero)):
        
        idx = nonzero[i]
        
        p1.append(float(p1_data[idx,p1_idx]))
        p2.append(float(p2_data[idx,p2_idx]))
        
        if p1_mine == True:
            p1_l_e.append(float(p1_errors[idx,p1_idx,0]))
            p1_u_e.append(float(p1_errors[idx,p1_idx,1]))
            
        if p1_mine == False:
            p1_e.append(float(p1_errors[idx]))
            
        if p2_mine == True:
            p2_l_e.append(float(p2_errors[idx,p2_idx,0]))
            p2_u_e.append(float(p2_errors[idx,p2_idx,1]))
            
        if p2_mine == False:
            p2_e.append(float(p2_errors[idx,p2_idx]))
    
        if pColor_idx != None:
            coloring.append(float(p1_data[idx, pColor_idx]))

    # ADJUST FOR ERRORS
    if p1_mine == True:
        p1_l_e = np.asarray(p1)-np.asarray(p1_l_e)
        p1_u_e = np.asarray(p1_u_e)-np.asarray(p1)
        
    if p2_mine == True:
        p2_l_e = np.asarray(p2)-np.asarray(p2_l_e)
        p2_u_e = np.asarray(p2_u_e)-np.asarray(p2)
        
        
    # DEFINE ERRORS
    if p1_mine == True:
        xerr = (p1_l_e, p1_u_e)
    if p1_mine == False:
        xerr = (p1_e)
    if p2_mine == True:
        yerr = (p2_l_e, p2_u_e)
    if p2_mine == False:
        yerr = (p2_e)
    
    if pColor_idx != None:
        sc = plt.scatter(p1, p2, 
                         c=coloring, vmin=np.min(coloring), vmax=np.max(coloring), 
                         s=65, cmap=cm, alpha = 0.5)
        plt.errorbar(p1, p2, xerr=xerr , yerr=yerr, 
                     linestyle="None", color='k', capsize=0, barsabove=False, alpha=0.15)   
                 
    if pColor_idx == None:
        
        
        
        sc = plt.scatter(p1, p2,  
                         s=65, alpha = 0.5)
        plt.errorbar(p1, p2, xerr=xerr, yerr=yerr, 
                     linestyle="None", color='k', capsize=0, barsabove=False, alpha=0.15)   
                     
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if pColor_idx != None:
        plt.colorbar(sc)
    plt.show()
