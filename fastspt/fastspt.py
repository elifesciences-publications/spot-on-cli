# A packaged version of the fastSPT code by Anders Sejr Hansen, Feb. 2016
# Python rewriting by MW, March 2017
#
# In this script we will compute jump length distributions, fit fastSPT data
# to the simple BOUND-UNBOUND 2-state model and determine the best-fit parameters
# 
# History: For the history of the script see the related CHANGELOG file.

## TODO MW (March 2017)
#- implement a meaningful logging system so that a log can be presented to the user
#- properly split this into several functions.
#- make a standalone version

## ==== Imports
import time
import lmfit
import numpy as np
from scipy.special import erfc

##
## ==== Auxiliary functions
##
def pdist(m):
    """Function returns the 2D distance between two points. Matrix form has:
       x1 y1
       x2 y2"""
    return ( (m[0,0]-m[1,0])**2+(m[0,1]-m[1,1])**2)**0.5

##
## ==== Main functions
##
def compute_jump_length_distribution(trackedPar,
                                     CDF=False, useAllTraj=False, TimePoints=8,
                                     GapsAllowed=1, JumpsToConsider=4,
                                     MaxJump=1.25, BinWidth=0.010):
    """Function that takes a series of translocations and computes an histogram of
    jump lengths. Returns both

    Arguments:
    - trackedPar: an object containing the trajectories
    - CDF (bool): compute the cumulative distribution function (CDF) instead of the probability distribution function (PDF)
    - useAllTraj (bool): True if we should use all trajectories to compute the histogram. This can lead to an overestimate of the bound fraction (see paper), but useful for troubleshooting
    - TimePoints (int): how many jump lengths to use for the fitting: 3 timepoints, yields 2 jumps
    - GapsAllowed (int): number of missing frames that are allowed in a single trajectory
    - JumpsToConsider (int): if `UseAllTraj` is False, then use no more than 3 jumps. 
    - TimeGap (float): time between frames in milliseconds;
    - MaxJump (float): for PDF fitting and plotting
    - BinWidth (float): for PDF fitting and plotting



    Returns:
    - An histogram at various \Delta t values.
    """

    PDF = not CDF
    tic = time.time() # Start the timer


    # Find total frames using a slight ad-hoc way: find the last frame with
    # a localization and round it. This is not an elegant solution, but it
    # works for your particle density:

    ## /!\ TODO MW: check this critical part of the code
    TempLastFrame = np.max([np.max(i[2]) for i in trackedPar]) #TempLastFrame = max(trackedPar(1,end).Frame)
    CellFrames = 100*round(TempLastFrame/100)
    CellLocs = sum([i[2].shape[1] for i in trackedPar]) # for counting the total number of localizations

    print "Number of frames: {}, number of localizations: {}".format(CellFrames, CellLocs)

    ##
    ## ==== Compile histograms for each jump lengths
    ##
    Min3Traj = 0   #for counting number of min3 trajectories;
    CellJumps = 0  # for counting the total number of jumps
    TransFrames = TimePoints+GapsAllowed*(TimePoints-1)
    TransLengths = []

    for i in range(TransFrames): # Initialize TransLengths
        TransLengths.append({"Step":[]}) # each iteration is a different number of timepoints
        
    if useAllTraj: ## Use all of the trajectory
        for i in range(len(trackedPar)): #1:length(trackedPar)
            CurrTrajLength = trackedPar[i][0].shape[0] #size(trackedPar(i).xy,1);

            if CurrTrajLength >= 3: # save lengths
                Min3Traj += 1

            #Now loop through the trajectory. Keep in mind that there are 
            #missing timepoints in the trajectory, so some gaps may be for 
            #multiple timepoints.

            #Figure out what the max jump to consider is:
            HowManyFrames = min(TimePoints-1, CurrTrajLength)
            if CurrTrajLength > 1:
                CellJumps = CellJumps+CurrTrajLength-1 # for counting all the jumps
                for n in range(1, HowManyFrames+1): #1:HowManyFrames
                    for k in range(CurrTrajLength-(n+1)): #=1:CurrTrajLength-n
                        # Find the current XY coordinate and frames between
                        # timepoints
                        CurrXY_points = np.vstack((trackedPar[i][0][k,:],
                                                   trackedPar[i][0][k+n,:]))
                        # vertcat(trackedPar(i).xy[k,:], trackedPar(i).xy[k+n,:])
                        CurrFrameJump = trackedPar[i][2][0][k+n] - \
                                        trackedPar[i][2][0][k]
                        #trackedPar(i).Frame(k+n) - trackedPar(i).Frame(k);

                        # Calculate the distance between the pair of points
                        TransLengths[CurrFrameJump-1]["Step"].append(pdist(CurrXY_points))

    elif not useAllTraj: ## Use only the first JumpsToConsider timepoints
        for i in range(len(trackedPar)): #1:length(trackedPar)
            CurrTrajLength =trackedPar[i][0].shape[0]#size(trackedPar(i).xy,1);
            if CurrTrajLength >= 3:
                Min3Traj += 1

            # Loop through the trajectory. If it is a short trajectory, you 
            # need to make sure that you do not overshoot. So first figure out 
            # how many jumps you can consider.

            # Figure out what the max jump to consider is:
            HowManyFrames = min([TimePoints-1, CurrTrajLength])
            if CurrTrajLength > 1:
                CellJumps = CellJumps+CurrTrajLength-1 #for counting all the jumps
                for n in range(1,HowManyFrames+1): #1:HowManyFrames
                    FrameToStop = min([CurrTrajLength, n+JumpsToConsider])
                    for k in range(FrameToStop-n): #=1:FrameToStop-n
                        # Find the current XY coordinate and frames between
                        # timepoints
                        CurrXY_points = np.vstack((trackedPar[i][0][k,:],
                                                   trackedPar[i][0][k+n,:]))

                        CurrFrameJump = trackedPar[i][2][0][k+n] - \
                                        trackedPar[i][2][0][k]
                        #trackedPar(i).Frame(k+n) - trackedPar(i).Frame(k);
                        # Compute the distance between the pair of points
                        TransLengths[CurrFrameJump-1]["Step"].append(pdist(CurrXY_points))

    ## Calculate the PDF histograms (required for CDF)
    HistVecJumps = np.arange(0, MaxJump+BinWidth, BinWidth)  # jump lengths in micrometers
    JumpProb = np.zeros((TimePoints-1, len(HistVecJumps))) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
    # TODO MW: investigate those differences
    for i in range(JumpProb.shape[0]):
        JumpProb[i,:] = np.float_(
            np.histogram(TransLengths[i]["Step"],
                         bins=np.hstack((HistVecJumps, HistVecJumps[-1])) )[0])/len(TransLengths[i]["Step"])
        
    if CDF: ## CALCULATE THE CDF HISTOGRAMS:
        BinWidthCDF = 0.001
        HistVecJumpsCDF = np.arange(0, MaxJump+BinWidthCDF, BinWidthCDF)  # jump lengths in micrometers
        JumpProbFine = np.zeros((TimePoints-1, len(HistVecJumpsCDF)))
        for i in range(JumpProb.shape[0]):
            JumpProbFine[i,:] = np.float_(
                np.histogram(TransLengths[i]["Step"],
                             bins=np.hstack((HistVecJumpsCDF, HistVecJumpsCDF[-1]) ))[0])/len(TransLengths[i]["Step"])
        
        JumpProbCDF = np.zeros((TimePoints-1, len(HistVecJumpsCDF))) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
        for i in range(JumpProbCDF.shape[0]): #1:size(JumpProbCDF,1)
            for j in range(2,JumpProbCDF.shape[1]+1): #=2:size(JumpProbCDF,2)
                JumpProbCDF[i,j-1] = sum(JumpProbFine[i,:j])

    toc = time.time()

    
    if PDF:
        return [HistVecJumps, JumpProb, {'time': toc-tic}]
    elif CDF:
        return [HistVecJumpsCDF,JumpProbCDF, HistVecJumps, JumpProb,
                {'time': toc-tic}]

def fit_jump_length_distribution(JumpProb, JumpProbCDF,
                                 HistVecJumps, HistVecJumpsCDF,
                                 LB, UB,
                                 LocError, iterations, dT, dZ, ModelFit, a, b,
                                 fit2states=True, fitSigma=False,
                                 verbose=True, init=None, useZcorr=True,
                                 solverparams = {'ftol':1e-20,
                                                 'xtol': 1e-20,
                                                 'maxfev': 100000,
                                             }):
    """Fits a kinetic model to an empirical jump length distribution.
    This applies a non-linear least squared fitting procedure.
    """
    ## Lower and Upper parameter bounds
    diff = np.array(UB) - np.array(LB)   # difference: used for initial parameters guess
    best_ssq2 = 5e10 # initial error

    # Need to ensure that the x-input is the same size as y-output
    if ModelFit == 1:
        ModelHistVecJumps = np.zeros((JumpProb.shape[0], len(HistVecJumps)))
        for i in range(JumpProb.shape[0]): # 1:size(JumpProb,1)
            ModelHistVecJumps[i,:] = HistVecJumps
        y = JumpProb
        x = np.repeat(HistVecJumps[:-1], JumpProb.shape[1])
            
    elif ModelFit == 2:
        ModelHistVecJumps = np.zeros((JumpProb.shape[0], len(HistVecJumpsCDF)))
        for i in range(JumpProb.shape[0]): # 1:size(JumpProb,1)
            ModelHistVecJumps[i,:] = HistVecJumpsCDF
        y = JumpProbCDF
        x = np.repeat(HistVecJumpsCDF[:-1], JumpProbCDF.shape[1]) 

    params = {"dT": dT,
              "dZ": dZ,
              "HistVecJumps": HistVecJumps,
              "HistVecJumpsCDF": HistVecJumpsCDF,
              "PDF_or_CDF": ModelFit,
              "JumpProb": JumpProb,
              "JumpProbCDF": JumpProbCDF,
              "LB": LB,
              "UB": UB,
              "a": a,
              "b": b,
              "useZcorr": useZcorr}

    ## ==== Get ready for the fitting
    def wrapped_jump_length_2states(x, D_free, D_bound, F_bound, sigma):
        """Wrapper for the main fit function (assuming global variables)"""
        if params['PDF_or_CDF'] == 1: # PDF fit
            dis = simulate_jump_length_distribution(
                (D_free, D_bound, F_bound), 
                params['JumpProb'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                sigma, #params['LocError'],
                params['PDF_or_CDF'],
                params['a'],
                params['b'], fit2states=True, useZcorr=params['useZcorr'], verbose=False)
        elif params['PDF_or_CDF'] == 2: # CDF fit
            dis = simulate_jump_length_distribution(
                (D_free, D_bound, F_bound), 
                params['JumpProbCDF'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                sigma, #params['LocError'],
                params['PDF_or_CDF'],
                params['a'],
                params['b'], fit2states=True, useZcorr=params['useZcorr'], verbose=False)
        return dis.flatten()
    
    def wrapped_jump_length_3states(x, D_fast, D_med, D_bound,
                                    F_fast, F_bound, sigma):
        """Wrapper for the main fit function (assuming global variables)"""
        if params['PDF_or_CDF'] == 1: # PDF fit
            dis = simulate_jump_length_distribution(
                (D_fast, D_med, D_bound, F_fast, F_bound), 
                params['JumpProb'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                sigma, #params['LocError'],
                params['PDF_or_CDF'],
                params['a'],
                params['b'], fit2states=False, useZcorr=params['useZcorr'], verbose=False)
        elif params['PDF_or_CDF'] == 2: # CDF fit
            dis = simulate_jump_length_distribution(
                (D_fast, D_med, D_bound, F_fast, F_bound), 
                params['JumpProbCDF'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                sigma, #params['LocError'],
                params['PDF_or_CDF'],
                params['a'],
                params['b'], fit2states=False, useZcorr=params['useZcorr'], verbose=False)
        if F_fast+F_bound<1:
            return np.hstack((dis.flatten(), 0))
        else:
            return np.hstack((dis.flatten(), 10000*(1-F_fast-F_bound)))
    
    if fit2states: # Instantiate model
        jumplengthmodel = lmfit.Model(wrapped_jump_length_2states) 
        pars = jumplengthmodel.make_params() ## Create an empty set of params
        y_init = y.flatten()
    else:
        jumplengthmodel = lmfit.Model(wrapped_jump_length_3states) 
        pars = jumplengthmodel.make_params() ## Create an empty set of params
        y_init = np.hstack((y.flatten(), 0))

    if init == None:
        init = []
        init1 = np.random.rand(len(LB))*diff+LB
        while not fit2states and len(init)<iterations:
            init1 = np.random.rand(len(LB))*diff+LB
            if init1[3]+init1[4]<1:
                init.append(init1)
        if fit2states:
            init = [np.random.rand(len(LB))*diff+LB for i in range(iterations)]
    elif len(init) != iterations:
        print "'iterations' variable ignored because 'init' is provided and has length {}".format(len(init))

    for (i,guess) in enumerate(init):
        eps = 1e-8
        if fit2states:
            if (guess.shape[0] != 3 and not fitSigma) or (guess.shape[0] != 4 and fitSigma): 
                print "init value has a wrong number of elements"
            pars['D_free'] .set(min=LB[0], max=UB[0], value=guess[0])
            pars['D_bound'].set(min=LB[1], max=UB[1], value=guess[1])
            pars['F_bound'].set(min=LB[2], max=UB[2], value=guess[2])
            if fitSigma:
                pars['sigma'].set(min=LB[3], max=UB[3], value=guess[3])
            else:
                pars['sigma'].set(value=LocError, vary=False)
            if abs(LB[0]-UB[0])<eps:
                pars['D_free'].set(value=LB[0], vary=False)
            if abs(LB[1]-UB[1])<eps:
                pars['D_bound'].set(value=LB[1], vary=False)
        else:
            if (guess.shape[0] != 5 and not fitSigma) or (guess.shape[0] != 6 and fitSigma):
                print "init value has a wrong number of elements"

            pars['D_fast'] .set(min=LB[0], max=UB[0], value=guess[0])
            pars['D_med']  .set(min=LB[1], max=UB[1], value=guess[1])
            pars['D_bound'].set(min=LB[2], max=UB[2], value=guess[2])            
            pars['F_bound'].set(min=LB[4], max=UB[4], value=guess[4])
            pars['F_fast'] .set(min=LB[3], max=UB[3], value=guess[3])
            if fitSigma:
                pars['sigma'].set(min=LB[5], max=UB[5], value=guess[5])
            else:
                pars['sigma'].set(value=LocError, vary=False)
            if abs(LB[0]-UB[0])<eps:
                pars['D_fast'].set(value=LB[0], vary=False)
            if abs(LB[1]-UB[1])<eps:
                pars['D_med'].set(value=LB[1], vary=False)
            if abs(LB[2]-UB[2])<eps:
                pars['D_bound'].set(value=LB[2], vary=False)

        out = jumplengthmodel.fit(y_init, x=x, params=pars, fit_kws=solverparams)
        ssq2 = (out.residual[:-1]**2).sum()/(out.residual.shape[0]-1)
        out.params.ssq2 = ssq2
        
        ## See if the current fit is an improvement:
        if ssq2 < best_ssq2:
            best_vals = out.params
            best_ssq2 = ssq2
            if verbose:
                print('==================================================')
                print('Improved fit on iteration {}'.format(i+1))
                print('Improved error is {}'.format(ssq2))
                print out.params.pretty_print(columns=['value', 'min', 'max', 'stderr'])
                print('==================================================')
        else:
            print('Iteration {} did not yield an improved fit'.format(i+1))
    return out

def _compute_2states(D_FREE, D_BOUND, F_BOUND, curr_dT, r, DeltaZ_use, LocError, useZcorr):
    """Subroutine for simulate_jump_distribution"""
    ## ==== Compute the integral
    HalfDeltaZ_use = DeltaZ_use/2.
    stp = DeltaZ_use/200.
    xint = np.linspace(-HalfDeltaZ_use, HalfDeltaZ_use, 200)
    yint = [C_AbsorBoundAUTO(i, curr_dT, D_FREE, HalfDeltaZ_use)*stp for i in xint]
    if useZcorr:
        Z_corr = 1/DeltaZ_use * np.array(yint).sum() # see below
    else:
        Z_corr=1
    
    # update the function output
    y1 = Z_corr * (1-F_BOUND) * (r/(2*(D_FREE*curr_dT+LocError**2)))
    y2 = np.exp( -(r**2)/(4*(D_FREE*curr_dT+LocError**2)))
    y3 = F_BOUND*(r /(2*(D_BOUND*curr_dT+LocError**2)))
    y4 = np.exp(-(r**2)/(4*(D_BOUND*curr_dT+LocError**2)))

    return y1*y2 + y3*y4
    
def _compute_3states(D_FAST, D_MED, D_BOUND, F_FAST, F_BOUND,
                     curr_dT, r, DeltaZ_useFAST, DeltaZ_useMED, LocError, useZcorr):
    """Subroutine for simulate_jump_distribution"""
    ## ==== Compute the integral
    HalfDeltaZ_useFAST = DeltaZ_useFAST/2.
    xintFAST = np.arange(-HalfDeltaZ_useFAST, HalfDeltaZ_useFAST, 4e-2)
    yintFAST = [C_AbsorBoundAUTO(i, curr_dT, D_FAST, HalfDeltaZ_useFAST)*4e-2 for i in xintFAST]
    if useZcorr:
        Z_corrFAST = 1/DeltaZ_useFAST * np.array(yintFAST).sum()
    else:
        Z_corrFAST = 1

    HalfDeltaZ_useMED = DeltaZ_useMED/2.
    xintMED = np.arange(-HalfDeltaZ_useMED, HalfDeltaZ_useMED, 4e-2)
    yintMED = [C_AbsorBoundAUTO(i, curr_dT, D_MED, HalfDeltaZ_useMED)*4e-2 for i in xintMED]
    if useZcorr:
        Z_corrMED = 1/DeltaZ_useMED * np.array(yintMED).sum()
    else:
        Z_corrMED = 1
    
    # update the function output
    y1 = F_BOUND*(r /(2*(D_BOUND*curr_dT+LocError**2)))
    y2 = np.exp(-(r**2)/(4*(D_BOUND*curr_dT+LocError**2)))
    y3 = Z_corrFAST * F_FAST * (r/(2*(D_FAST*curr_dT+LocError**2)))
    y4 = np.exp( -(r**2)/(4*(D_FAST*curr_dT+LocError**2)))
    y5 = Z_corrMED * (1-F_FAST-F_BOUND) * (r/(2*(D_MED*curr_dT+LocError**2)))
    y6 = np.exp( -(r**2)/(4*(D_MED*curr_dT+LocError**2)))
    
    return y1*y2 + y3*y4 + y5*y6

def simulate_jump_length_distribution(parameter_guess, JumpProb,
                                      HistVecJumpsCDF, HistVecJump,
                                      dT, dZ, LocError, PDF_or_CDF, a, b,
                                      fit2states = True, useZcorr=True, verbose=True):
    """Function 'SS_2State_model_Z_corr_v4' actually returns a distribution
    given the parameter_guess input. This function is to be used inside a
    least square fitting method, such as Matlab's `lsqcurvefit` or 
    Python's `lmfit`.
    
    Note that this function assumes some *global variables* that are provided
    by the main script: LocError dT HistVecJumps dZ HistVecJumpsCDF PDF_or_CDF
    """

    # ==== Initialize stuff
    HistVecJumps = HistVecJump.copy()
    HistVecJumps += HistVecJumps[1]/2.
    r = HistVecJumpsCDF.copy()
    r += r[1]/2.
    y = np.zeros((JumpProb.shape[0], len(r)))
    Binned_y_PDF = np.zeros((JumpProb.shape[0], JumpProb.shape[1]))

    if fit2states:
        D_FREE  = parameter_guess[0]
        D_BOUND = parameter_guess[1]
        F_BOUND = parameter_guess[2]
    else:
        D_FAST  = parameter_guess[0]
        D_MED   = parameter_guess[1]
        D_BOUND = parameter_guess[2]
        F_FAST  = parameter_guess[3]
        F_BOUND = parameter_guess[4]

    # ==== Precompute stuff
    # Calculate the axial Z-correction
    # First calculate the corrected DeltaZ:
    ##DeltaZ_use = dZ + 0.15716  * D_FREE**.5 + 0.20811 # See CHANGELOG_fit
    ##DeltaZ_use = dZ + 0.24472 * D_FREE**.5 + 0.19789
    if fit2states:
        DeltaZ_use = dZ + a * D_FREE**.5 + b #HalfDeltaZ_use = DeltaZ_use/2
    else:
        DeltaZ_useFAST = dZ + a * D_FAST**.5 + b #HalfDeltaZ_use = DeltaZ_use/2
        DeltaZ_useMED = dZ + a * D_MED**.5 + b #HalfDeltaZ_use = DeltaZ_use/2

    for iterator in range(JumpProb.shape[0]):
        # Calculate the jump length distribution of the parameters for each
        # time-jump
        curr_dT = (iterator+1)*dT
        if verbose:
            print "-- computing dT = {} ({}/{})".format(curr_dT, iterator+1, JumpProb.shape[0])
        if fit2states:
            y[iterator,:] = _compute_2states(D_FREE, D_BOUND, F_BOUND,
                                             curr_dT, r, DeltaZ_use, LocError, useZcorr)
        else:
            y[iterator,:] = _compute_3states(D_FAST, D_MED, D_BOUND,
                                             F_FAST, F_BOUND,
                                             curr_dT, r, DeltaZ_useFAST,
                                             DeltaZ_useMED, LocError, useZcorr)

    if PDF_or_CDF == 1:
        # norm_y = np.zeros_like(y)
        # for i in range(y.shape[0]): # Normalize y as a PDF
        #     norm_y[i,:] = y[i,:]/y[i,:].sum()        
	#Now bin the output y so that it matches the JumpProb variable:
	for i in range(JumpProb.shape[0]): #1:size(JumpProb,1)
	    for j in range(JumpProb.shape[1]): #=1:size(JumpProb,2)
		if j == (JumpProb.shape[1]-1):
		    Binned_y_PDF[i,j] = y[i,maxIndex:].mean()
		else:
                    minIndex = np.argmin(np.abs(r-HistVecJumps[j]))
                    maxIndex = np.argmin(np.abs(r-HistVecJumps[j+1]))
		    Binned_y_PDF[i,j] = y[i,minIndex:maxIndex].mean()
	for i in range(JumpProb.shape[0]): #1:size(JumpProb,1) ## Normalize
	    Binned_y_PDF[i,:] = Binned_y_PDF[i,:]/sum(Binned_y_PDF[i,:]);
        Binned_y = Binned_y_PDF #You want to fit to a histogram, so no need to calculate the CDF
        return Binned_y

    elif PDF_or_CDF == 2:
        # You want to fit to a CDF function, so first we must calculate the CDF
        # from the finely binned PDF
        Binned_y_CDF = np.zeros((JumpProb.shape[0], JumpProb.shape[1]))

        ## Normalize the PDF
        for i in range(Binned_y_CDF.shape[0]):
            Binned_y_PDF[i,:] = y[i,:]/y[i,:].sum()
    
        ## calculate the CDF
        for i in range(Binned_y_CDF.shape[0]): #1:size(Binned_y_CDF,1):
            Binned_y_CDF[i,:] = np.cumsum(Binned_y_PDF[i,:])
            #for j in range(1, Binned_y_CDF.shape[1]): #=2:size(Binned_y_CDF,2):
            #    Binned_y_CDF[i,j] = Binned_y_PDF[i,:j].sum()
        Binned_y = Binned_y_CDF ##Output the final variable

    return Binned_y

def generate_jump_length_distribution(fitparams, JumpProb, r,
                                      LocError, dT, dZ, a, b, fit2states=True,
                                      norm=False, useZcorr=True):
    """
    This function has no docstring. This is bad
    """
    if fit2states:
        D_free = fitparams['D_free']
        D_bound = fitparams['D_bound']
        F_bound = fitparams['F_bound']
    else:
        D_fast = fitparams['D_fast']
        D_med = fitparams['D_med']
        D_bound = fitparams['D_bound']
        F_fast = fitparams['F_fast']
        F_bound = fitparams['F_bound']

    y = np.zeros((JumpProb.shape[0], r.shape[0]))
    #Z_corr = np.zeros(JumpProb.shape[0]) # Assume ABSORBING BOUNDARIES

    # Calculate the axial Z-correction
    if fit2states:
        DeltaZ_use = dZ + a * D_free**.5 + b #HalfDeltaZ_use = DeltaZ_use/2
    else:
        DeltaZ_useFAST = dZ + a * D_fast**.5 + b #HalfDeltaZ_use = DeltaZ_use/2
        DeltaZ_useMED = dZ + a * D_med**.5 + b #HalfDeltaZ_use = DeltaZ_use/2

    for iterator in range(JumpProb.shape[0]):
        # Calculate the jump length distribution of the parameters for each
        # time-jump
        curr_dT = (iterator+1)*dT

        if fit2states:
            y[iterator,:] = _compute_2states(D_free, D_bound, F_bound,
                                             curr_dT, r, DeltaZ_use, LocError, useZcorr)
        else:
            y[iterator,:] = _compute_3states(D_fast, D_med, D_bound,
                                             F_fast, F_bound,
                                             curr_dT, r, DeltaZ_useFAST,
                                             DeltaZ_useMED, LocError, useZcorr)

        if norm:
            norm_y = np.zeros_like(y)
            for i in range(y.shape[0]): # Normalize y as a PDF
                norm_y[i,:] = y[i,:]/y[i,:].sum()
            
            #scaled_y = (float(len(HistVecJumpsCDF))/len(HistVecJumps))*norm_y #scale y for plotting next to histograms
            y = norm_y
    return y


def C_AbsorBoundAUTO(z, CurrTime, D, halfZ):
    """
    This is a corrected version of equation 16 in Kues and Kubitscheck, Single
    Molecules, 2002 and a corrected version of equation Suppl 5.7 in Mazza et
    al, Nucleic Acids Research, 2012. Both equations are wrong, but they are
    wrong in different ways and the one below is correct. 
    Moreover, this implementation automatically stops the sum when the error
    is negligble. 

    Original Matlab code in SS_2State_model_Z_corr_v4.m
    """

    WhenToStop = 1e-10
    f = np.inf
    n = 0 # iterator
    h = 1
    
    while np.abs(f) > WhenToStop:
        if CurrTime != 0:
            z1 =  ((2*n+1)*halfZ-z)/(4*D*CurrTime)**.5
            z2 =  ((2*n+1)*halfZ+z)/(4*D*CurrTime)**.5
        elif (2*n+1)*halfZ-z<0:
            z1 = -np.inf
            z2 = np.inf
        else:
            z1 = np.inf
            z2 = np.inf
            
        f = ((-1)**n) * ( erfc(z1) +  erfc(z2) )
        h -= f
        n += 1
    return h
