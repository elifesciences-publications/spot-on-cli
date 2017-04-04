import sys

sys.path.append("/home/maxime/code/fastSPT/")
import fastSPT, fastSPT_tools

sys.path.append("/home/maxime/thesis/9_SPT/fastSPT/SPTGUI/")
import parsers

def gui_to_fastSPT_conversion(fn):
    """A test to make sure that the files we convert from the GUI
    format are usable for fastSPT"""
    with open(fn, 'r') as f:
        cell = parsers.to_fastSPT(f)
        fastSPT.compute_jump_length_distribution(cell, CDF=True, useAllTraj=True)
        
if __name__ == "__main__":
    fn = "/home/maxime/Bureau/Thesis/9_SPT/fastSPT/uploads/uploads/20160526_mESC_C87_Halo-mCTCF_25nM_PA-JF646_1ms633_3-405_4msCam_cell1._agbfilG.parsed"
    
    gui_to_fastSPT_conversion(fn)
