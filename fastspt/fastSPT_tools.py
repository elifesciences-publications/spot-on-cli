## fastSPT_tools
## Some tools for the fastSPT package
## By MW, GPLv3+
## March 2017

## ==== Imports
import pickle, sys, scipy.io, os
import numpy as np

## ==== Sample dataset-related functions
def list_sample_datasets(path):
    """Simple relay function that allows to list datasets from a datasets.py file"""
    sys.path.append(path)
    import datasets
    reload(datasets) # Important I think
    return datasets.list(path, string=True)

def load_dataset(path, datasetID, cellID):
    """Simple helper function to load one or several cells from a dataset"""
    ## Get the information about the datasets
    sys.path.append(path)
    import datasets
    reload(datasets) # Important I think
    li = datasets.list(path, string=False)

    if type(cellID) == int:
        cellID = [cellID]
    
    try: ## Check if our dataset(s) is/are available
        for cid in cellID:
            if not li[1][datasetID][cid].lower() == "found":
                raise IOError("This dataset does not seem to be available. Either it couldn't be found or it doesn't exist in the database.")
    except:
        raise IOError("This dataset does not seem to be available. Either it couldn't be found or it doesn't exist in the database or there is a problem with the database.")

    da_info = li[0][datasetID]

    ## Load the datasets
    AllData = []
    for ci in cellID:
        mat = scipy.io.loadmat(os.path.join(path,
                                            da_info['path'],
                                            da_info['workspaces'][ci]))
        AllData.append(np.asarray(mat['trackedPar'][0]))
    return np.hstack(AllData) ## Concatenate them before returning

def load_dataset_from_path(path):
    """Returns a dataset object from a Matlab file"""
    mat = scipy.io.loadmat(path)
    return np.asarray(mat['trackedPar'][0])
