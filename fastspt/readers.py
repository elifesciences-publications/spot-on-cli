## readers.py is part of the fastspt library
## By MW, GPLv3+, Oct. 2017
## readers.py imports many file formats widespread in SPT analysis
## Imported from Spot-On

import scipy.io, os, json, xmltodict
import numpy as np
import pandas as pd

## ==== evalSPT
def read_evalspt(fn, framerate, pixelsize):
    return read_arbitrary_csv(fn, col_traj=3, col_x=0, col_y=1, col_frame=2,
                              framerate=framerate/1000., pixelsize=pixelsize/1000.,
                              sep="\t", header=None)

## ==== MOSAIC suite file format
def read_mosaic(fn, framerate, pixelsize):
    return read_arbitrary_csv(fn, col_traj="Trajectory", col_x="x", col_y="y", col_frame="Frame", framerate=framerate/1000., pixelsize=pixelsize/1000.)

## ==== TrackMate file format
def read_trackmate_csv(fn, framerate):
    """Do not call directly, wrapped into `read_trackmate`"""    
    def cb(da):
        return da[da.TRACK_ID!="None"]
    return read_arbitrary_csv(fn, col_traj="TRACK_ID", col_x="POSITION_X", col_y="POSITION_Y", col_frame="FRAME", framerate=framerate/1000., cb=cb)

def read_trackmate_xml(fn):
    """Do not call directly, wrapped into `read_trackmate`."""
    x=xmltodict.parse(open(fn, 'r').read())
    # Checks
    spaceunit = x['Tracks']['@spaceUnits']
    if spaceunit not in ('micron', 'um'):
        raise IOError("Spatial unit not recognized: {}".format(spaceunit))
    if x['Tracks']['@timeUnits'] != 'ms':
        raise IOError("Time unit not recognized")
    
    # parameters
    framerate = float(x['Tracks']['@frameInterval'])/1000. # framerate in ms
    traces = []
    
    for particle in x['Tracks']['particle']:  
        traces.append([(float(d['@x']), float(d['@y']), float(d['@t'])*framerate, float(d['@t'])) for d in particle['detection']])
    return traces

## ==== CSV file format
def read_csv(fn):
    return read_arbitrary_csv(fn, col_traj="trajectory", col_x="x", col_y="y", col_frame="frame", col_t="t")

## ==== Anders' file format
def read_anders(fn, new_format=True):
    """The file format sent by Anders. I don't really know where it 
    comes from.
    new_format tells whether we should perform weird column manipulations
    to get it working again..."""
    
    def new_format(cel):
        """Converts between the old and the new Matlab format. To do so, it 
        swaps columns 1 and 2 of the detections and transposes the matrices"""
        cell = cel.copy()
        for i in range(len(cell)):
            f = cell[i][2].copy()
            cell[i][2] = cell[i][1].T.copy()
            cell[i][1] = f.T
        return cell
    
    ## Sanity checks
    if not os.path.isfile(fn):
        raise IOError("File not found: {}".format(fn))

    try:
        mat = scipy.io.loadmat(fn)
        m=np.asarray(mat['trackedPar'])
    except:
        raise IOError("The file does not seem to be a .mat file ({})".format(fn))

    if new_format:
        m[0] = new_format(m[0])
    
    ## Make the conversion
    traces_header = ('x','y','t','f')
    traces = []
    for tr in m[0]:
        x = [float(i) for i in tr[0][:,0]]
        y = [float(i) for i in tr[0][:,1]]
        t = [float(i) for i in tr[1][0]]
        f = [int(i) for i in tr[2][0]]
        traces.append(zip(x,y,t,f))
    return traces

## ==== Format for fastSPT
def to_fastSPT(f):
    """Returns an object formatted to be used with fastSPT from a parsed dataset
    (in the internal representation of the GUI). f is a file descriptor (thus the
    function assumes that the file exists).

    Actually, the fastSPT code is a little bit picky about what it likes and what
    it doesn't. It cares strictly about the file format, that is a nested numpy
    object, and of the data types. I expect many bugs to arise from improper 
    converters that do not fully comply with the file format."""

    da = json.loads(f.read()) ## Load data

    ## Create the object
    dt = np.dtype([('xy', 'O'), ('TimeStamp', 'O'), ('Frame', 'O')]) # dtype
    DT = np.dtype('<f8', '<f8', 'uint16')
    trackedPar = []
    for i in da:
        xy = []
        TimeStamp = []
        Frame = []
        for p in i:
            xy.append([p[0],p[1]])
            TimeStamp.append(p[2])
            Frame.append(p[3])
        trackedPar.append((np.array(xy, dtype='<f8'),
                           np.array([TimeStamp], dtype='<f8'),
                           np.array([Frame], dtype='uint16')))
    return np.asarray(trackedPar, dtype=dt)

##
## ==== This are some helper functions
##

def traces_to_csv(traces):
    """Returns a CSV file with the format 
    trajectory,x,y,t,frame
    """
    csv = "trajectory,x,y,t,frame\n"
    for (tr_n, tr) in enumerate(traces):
        for pt in tr:
            csv +="{},{},{},{},{}\n".format(tr_n, pt[0],pt[1],pt[2],pt[3])
    return csv

def read_arbitrary_csv(fn, col_x="", col_y="", col_frame="", col_t="t",
                       col_traj="", framerate=None, pixelsize=None, cb=None,
                       sep=",", header='infer'):
    """This function takes the file name of a CSV file as input and parses it to
    the list of list format required by Spot-On. This function is called by various
    CSV importers and it is advised not to call it directly."""
    
    da = pd.read_csv(fn, sep=sep, header=header) # Read file
    
    # Check that all the columns are present:
    cols = da.columns
    if (not (col_traj in cols and col_x in cols and col_y in cols and col_frame in cols)) or (not (col_t in cols) and framerate==None):
        raise IOError("Missing columns in the file, or wrong header")
        
    # Correct units if needed
    if framerate != None:
        da[col_t]=da[col_frame]*framerate
    if pixelsize != None:
        da[col_x]*=pixelsize
        da[col_y]*=pixelsize
        
    # Apply potential callback
    if cb != None:
        da = cb(da)
        
    # Split by traj
    out = []
    for (idx,t) in da.sort_values(col_traj).groupby(col_traj):
        tr = [(tt[1][col_x], tt[1][col_y], tt[1][col_t], int(tt[1][col_frame])) for tt in t.sort_values(col_frame).iterrows()] # Order by trace, then by frame
        out.append(tr)
    
    return out





def lol():
    return "Why so serious?"

