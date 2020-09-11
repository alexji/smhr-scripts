import numpy as np
from smh import Session

import optparse
import sys, os

if __name__=="__main__":
    assert len(sys.argv) == 4, "Usage: python transfer_spectral_models.py smh_file_with_masks.smh smh_file_needing_masks.smh output_smh_file.smh"
    origfname = sys.argv[1]
    newfname = sys.argv[2]
    outfname = sys.argv[3]
    assert os.path.exists(fname), "input {} does not exist".format(fname)
    assert not os.path.exists(outfname), "output {} already exists".format(outfname)

    session = Session.load(fname)
    
    
