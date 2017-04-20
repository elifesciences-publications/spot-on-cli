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
                                     TimeGap=4.477, MaxJump=1.25, BinWidth=0.010):
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
    # CellNumb = 0
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
    #print('2. compiling histogram of jump lengths...') # DBG
        
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
                        # Calculate the distance between the pair of points
                        TransLengths[CurrFrameJump-1]["Step"].append(pdist(CurrXY_points))

    ## Calculate the PDF histograms (required for CDF)
    HistVecJumps = np.arange(0, MaxJump+BinWidth, BinWidth)  # jump lengths in micrometers
    JumpProb = np.zeros((TimePoints-1, len(HistVecJumps))) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
    # TODO MW: investigate those differences
    for i in range(JumpProb.shape[0]):
        JumpProb[i,:] = np.float_(
            np.histogram(TransLengths[i]["Step"],
                         bins=np.hstack((HistVecJumps, MaxJump)))[0])/len(TransLengths[i]["Step"])
        
    if CDF: ## CALCULATE THE CDF HISTOGRAMS:
        BinWidthCDF = 0.001
        HistVecJumpsCDF = np.arange(0, MaxJump+BinWidthCDF, BinWidthCDF)  # jump lengths in micrometers
        JumpProbFine = np.zeros((TimePoints-1, len(HistVecJumpsCDF)))
        for i in range(JumpProb.shape[0]):
            JumpProbFine[i,:] = np.float_(
                np.histogram(TransLengths[i]["Step"],
                             bins=np.hstack((HistVecJumpsCDF, MaxJump)))[0])/len(TransLengths[i]["Step"])
        
        JumpProbCDF = np.zeros((TimePoints-1, len(HistVecJumpsCDF))) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
        for i in range(JumpProbCDF.shape[0]): #1:size(JumpProbCDF,1)
            for j in range(2,JumpProbCDF.shape[1]): #=2:size(JumpProbCDF,2)
                JumpProbCDF[i,j] = sum(JumpProbFine[i,:j])

    toc = time.time()

    
    if PDF:
        return [HistVecJumps, JumpProb, {'time': toc-tic}]
    elif CDF:
        return [HistVecJumpsCDF,JumpProbCDF, HistVecJumps, JumpProb,
                {'time': toc-tic}]


def fit_jump_length_distribution(JumpProb, JumpProbCDF,
                                 HistVecJumps, HistVecJumpsCDF,
                                 D_Free, D_Bound, Frac_Bound,
                                 LocError, iterations, dT, dZ, ModelFit,
                                 verbose=True):
    """Fits a kinetic model to an empirical jump length distribution.

    This applies a non-linear least squared fitting procedure.
    """
    ## Lower and Upper parameter bounds
    ## /!\ TODO MW: these are specific to the Matlab least square solver
    LB = np.array([D_Free[0], D_Bound[0], Frac_Bound[0]])
    UB = np.array([D_Free[1], D_Bound[1], Frac_Bound[1]])
    diff = UB - LB   # difference: used for initial parameters guess
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

    params = {"LocError": LocError,
              "dT": dT,
              "dZ": dZ,
              "HistVecJumps": HistVecJumps,
              "HistVecJumpsCDF": HistVecJumpsCDF,
              "PDF_or_CDF": ModelFit,
              "JumpProb": JumpProb,
              "JumpProbCDF": JumpProbCDF,
              "LB": LB,
              "UB": UB}

    ## ==== Get ready for the fitting
    def wrapped_jump_length(x, D_free, D_bound, F_bound):
        """Wrapper for the main fit function (assuming global variables)"""
        if params['PDF_or_CDF'] == 1: # PDF fit
            dis = simulate_jump_length_distribution(
                (D_free, D_bound, F_bound), 
                params['JumpProb'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                params['LocError'],
                params['PDF_or_CDF'], verbose=False)
        elif params['PDF_or_CDF'] == 2: # CDF fit
            dis = simulate_jump_length_distribution(
                (D_free, D_bound, F_bound), 
                params['JumpProbCDF'],
                params['HistVecJumpsCDF'],
                params['HistVecJumps'],
                params['dT'],
                params['dZ'],
                params['LocError'],
                params['PDF_or_CDF'], verbose=False)
        return dis.flatten()
    
    jumplengthmodel = lmfit.Model(wrapped_jump_length) ## Instantiate model
    pars = jumplengthmodel.make_params() ## Create an empty set of params

    for i in range(iterations):
        guess = np.random.rand(len(LB))*diff+LB ## Guess a random set of parameters
        pars['D_free'] .set(min=LB[0], max=UB[0], value=guess[0])
        pars['D_bound'].set(min=LB[1], max=UB[1], value=guess[1])
        pars['F_bound'].set(min=LB[2], max=UB[2], value=guess[2])

        out = jumplengthmodel.fit(y.flatten(), x=x, params=pars) ## FIT
        ssq2 = (out.residual**2).sum()/out.residual.shape[0]

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
                print('Cell {}: Iteration {} did not yield an improved fit'.format(CellNumb+1, i+1))
    return out

def simulate_jump_length_distribution(parameter_guess, JumpProb,
                              HistVecJumpsCDF, HistVecJumps,
                              dT, dZ, LocError, PDF_or_CDF,
                              verbose=True):
    """Function 'SS_2State_model_Z_corr_v4' actually returns a distribution
    given the parameter_guess input. This function is to be used inside a
    least square fitting method, such as Matlab's `lsqcurvefit`.
    
    Note that this function assumes some *global variables* that are provided
    by the main script: LocError dT HistVecJumps dZ HistVecJumpsCDF PDF_or_CDF
    """

    # ==== Initialize stuff    
    r = HistVecJumpsCDF
    y = np.zeros((JumpProb.shape[0], len(r)))
    Binned_y_PDF = np.zeros((JumpProb.shape[0], JumpProb.shape[1]))

    D_FREE  = parameter_guess[0]
    D_BOUND = parameter_guess[1]
    F_BOUND = parameter_guess[2]

    # Assume ABSORBING BOUNDARIES
    Z_corr = np.zeros(JumpProb.shape[0])

    # ==== Precompute stuff
    # Calculate the axial Z-correction
    # First calculate the corrected DeltaZ:
    DeltaZ_use = dZ + 0.15716  * D_FREE**.5 + 0.20811 # See CHANGELOG_fit
    HalfDeltaZ_use = DeltaZ_use/2    

    for iterator in range(JumpProb.shape[0]): #=1:size(JumpProb,1):
        # Calculate the jump length distribution of the parameters for each
        # time-jump
        curr_dT = (iterator+1)*dT
        if verbose:
            print "-- computing dT = {} ({}/{})".format(curr_dT, iterator+1, JumpProb.shape[0])
        
        ## ==== Compute the integral
        # Getting the theoretical bounds of the integral so that I can properly
        # compute in the first place the erfc
        # - z \in [-HalfDeltaZ_use, HalfDeltaZ_use]
        # - currTime \in [0, JumpProb.shape[0]*dT -> [0, 0.03]
        # - D \in [0.15, 25]
        # - halfZ \in [0.3, 1.7]
        #print JumpProb.shape[0]*dT
        #print dZ + 0.15716 * 0.15**.5 + 0.20811, dZ + 0.15716 * 25**.5 + 0.20811

        xint = np.arange(-HalfDeltaZ_use, HalfDeltaZ_use, 1e-3)
        yint = [C_AbsorBoundAUTO(i, curr_dT, D_FREE, HalfDeltaZ_use)*1e-3 for i in xint]
        
        Z_corr[iterator] = 1/DeltaZ_use * np.array(yint).sum() # see below
        #1/DeltaZ_use * integral(@(z)C_AbsorBoundAUTO(z,curr_dT, D_FREE, HalfDeltaZ_use),-HalfDeltaZ_use,HalfDeltaZ_use);
    
        # update the function output
        #print Z_corr[iterator], (1-F_BOUND), D_FREE, curr_dT, LocError, D_BOUND
        y1 = Z_corr[iterator] * (1-F_BOUND) * (r/(2*(D_FREE*curr_dT+LocError**2)))
        y2 = np.exp( -(r**2)/(4*(D_FREE*curr_dT+LocError**2)))
        y3 = F_BOUND*(r /(2*(D_BOUND*curr_dT+LocError**2)))
        y4 = np.exp(-(r**2)/(4*(D_BOUND*curr_dT+LocError**2)))
        
        y[iterator,:] = y1*y2 + y3*y4
        #Z_corr[iterator] * (1-F_BOUND) * (r/(2*(D_FREE*curr_dT+LocError^2))) * exp(-(r^2.)/(4*(D_FREE*curr_dT+LocError^2.))) + F_BOUND*(r/(2*(D_BOUND*curr_dT+LocError^2.)))*exp(-(r^2.)/(4*(D_BOUND*curr_dT+LocError^2.)))

    if PDF_or_CDF == 1:
	## Now bin the output y so that it matches the JumpProb variable: 
	for i in range(JumpProb.shape[0]): #1:size(JumpProb,1)
	    for j in range(JumpProb.shape[1]): #=1:size(JumpProb,2)
		if j == JumpProb.shape[1]:
		    Binned_y_PDF[i,j] = y[i,maxIndex:].mean()
		else:
		    #[~, minIndex] = min(abs(r-HistVecJumps(j)));
		    #[~, maxIndex] = min(abs(r-HistVecJumps(j+1)));
                    minIndex = np.argmin(np.abs(r-HistVecJumps[j]))
                    maxIndex = np.argmin(np.abs(r-HistVecJumps[j+1]))
		    Binned_y_PDF[i,j] = y[i,minIndex:maxIndex].mean()

	## Normalize:
	for i in range(JumpProb.shape[0]): #1:size(JumpProb,1)
	    Binned_y_PDF[i,:] = Binned_y_PDF[i,:]/sum(Binned_y_PDF[i,:]);

        Binned_y = Binned_y_PDF #You want to fit to a histogram, so no need to calculate the CDF

    elif PDF_or_CDF == 2:
        # You want to fit to a CDF function, so first we must calculate the CDF
        # from the finely binned PDF
        Binned_y_CDF = np.zeros((JumpProb.shape[0], JumpProb.shape[1]))

        ## Normalize the PDF
        for i in range(Binned_y_CDF.shape[0]):
            Binned_y_PDF[i,:] = y[i,:-1]/y[i,:-1].sum()
    
        ## calculate the CDF
        for i in range(Binned_y_CDF.shape[0]): #1:size(Binned_y_CDF,1):
            for j in range(1, Binned_y_CDF.shape[1]): #=2:size(Binned_y_CDF,2):
                Binned_y_CDF[i,j] = Binned_y_PDF[i,:j].sum()

        Binned_y = Binned_y_CDF ##Output the final variable

    return Binned_y

def generate_jump_length_distribution(D_free, D_bound, F_bound, JumpProb, r,
                                LocError, dT, dZ):
    """
    Anybody really interested in what is really
    happening would notice that this function is only generating a distribution
    and thus this code is in a way redundant with the one of the fitting procedure.
    /!\ TODO MW: remove this redundancy
    """

    y = np.zeros((JumpProb.shape[0], r.shape[0]))
    Z_corr = np.zeros(JumpProb.shape[0]) # Assume ABSORBING BOUNDARIES

    # Calculate the axial Z-correction
    DeltaZ_use = dZ + 0.15716  * D_free**.5 + 0.20811 # See CHANGELOG_fit
    HalfDeltaZ_use = DeltaZ_use/2

    for iterator in range(JumpProb.shape[0]):
        # Calculate the jump length distribution of the parameters for each
        # time-jump
        curr_dT = (iterator+1)*dT

        ## ==== Compute the integral
        xint = np.arange(-HalfDeltaZ_use, HalfDeltaZ_use, 1e-3)
        yint = [C_AbsorBoundAUTO(i, curr_dT, D_free, HalfDeltaZ_use)*1e-3 for i in xint]
        Z_corr[iterator] = 1/DeltaZ_use * np.array(yint).sum()
        #Z_corr(1,iterator) =  1/DeltaZ_use * integral(@(z)C_AbsorBoundAUTO(z,curr_dT, D_free, HalfDeltaZ_use),-HalfDeltaZ_use,HalfDeltaZ_use);
    
        ## ==== update the function output
        y1 = Z_corr[iterator] * (1-F_bound) * (r/(2*(D_free*curr_dT+LocError**2)))
        y2 = np.exp( -(r**2)/(4*(D_free*curr_dT+LocError**2)))
        y3 = F_bound*(r /(2*(D_bound*curr_dT+LocError**2)))
        y4 = np.exp(-(r**2)/(4*(D_bound*curr_dT+LocError**2)))
        
        y[iterator,:] = y1*y2 + y3*y4        
        #y(iterator,:) = Z_corr(1,iterator) .* (1-F_bound).*(r./(2*(D_free*curr_dT+LocError^2))).*exp(-r.^2./(4*(D_free*curr_dT+LocError^2))) + F_bound.*(r./(2*(D_bound*curr_dT+LocError^2))).*exp(-r.^2./(4*(D_bound*curr_dT+LocError^2))) ;
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

    WhenToStop = 1e-6 # 1e-10
    f = np.inf
    n = 0 # iterator
    h = 1
    
    while np.abs(f) > WhenToStop:
        #print ((2*n+1)*halfZ-z)/(4*D*CurrTime)**.5, ((2*n+1)*halfZ+z)/(4*D*CurrTime)**.5 ## DBG
        f = ((-1)^n) * ( erfc( ((2*n+1)*halfZ-z)/(4*D*CurrTime)**.5) +  erfc( ((2*n+1)*halfZ+z)/(4*D*CurrTime)**.5) )
        h -= f
        n += 1
    return h
