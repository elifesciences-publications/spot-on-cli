# A packaged version of the fastSPT code by Anders Sejr Hansen, Feb. 2016
# Python rewriting by MW, March 2017
#
# In this script we will fit fastSPT data to the simple BOUND-UNBOUND
# 2-state model and determine the best-fit parameters
# 
# History: For the history of the script see the related CHANGELOG file.

## TODO MW (March 2017)
#- implement a meaningful logging system so that a log can be presented to the user
#- properly split this into several functions.
#- make a standalone version

## ==== Imports
import time
import numpy as np

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
    HistVecJumps = np.arange(0, MaxJump, BinWidth)  # jump lengths in micrometers
    JumpProb = np.zeros((TimePoints-1, len(HistVecJumps)-1)) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
    # TODO MW: investigate those differences
    for i in range(JumpProb.shape[0]):
        JumpProb[i,:] = np.float_(np.histogram(TransLengths[i]["Step"], bins=HistVecJumps)[0])/len(TransLengths[i]["Step"]) ## /!\ TODO MW: check that Matlab histc and np.histogram are equivalent
        
    if CDF: ## CALCULATE THE CDF HISTOGRAMS:
        HistVecJumpsCDF = np.arange(0, MaxJump, 0.001)  # jump lengths in micrometers
        JumpProbFine = np.zeros((TimePoints-1, len(HistVecJumpsCDF)-1))
        
        for i in range(JumpProb.shape[0]):
            JumpProbFine[i,:] = np.float_(np.histogram(TransLengths[i]["Step"], bins=HistVecJumpsCDF)[0])/len(TransLengths[i]["Step"]) ## /!\ TODO MW: check that Matlab histc and np.histogram are equivalent
        
        JumpProbCDF = np.zeros((TimePoints-1, len(HistVecJumpsCDF)-1)) ## second -1 added when converting to Python due to inconsistencies between the histc and histogram function
        for i in range(JumpProbCDF.shape[0]): #1:size(JumpProbCDF,1)
            for j in range(2,JumpProbCDF.shape[1]): #=2:size(JumpProbCDF,2)
                JumpProbCDF[i,j] = sum(JumpProbFine[i,:j])

    toc = time.time()

    
    if PDF:
        return [JumpProb, {'time': toc-tic}]
    elif CDF:
        return [JumpProbCDF, JumpProb, {'time': toc-tic}]
