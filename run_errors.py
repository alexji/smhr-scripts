#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import glob, os, sys, time
from smh import Session

if __name__=="__main__":
    #smh_fnames = glob.glob("*.smh")
    #
    fname = sys.argv[1]
    newfname = sys.argv[2]
    assert not os.path.exists(newfname), f"{newfname} already exists! Stopping..."
    #eTeff = float(sys.argv[2])
    #elogg = float(sys.argv[3])
    #evt = float(sys.argv[4])
    #eMH = float(sys.argv[5])
    assert os.path.exists(fname)
    smh_fnames = [fname]
    for fname in smh_fnames:
        print("Processing {}".format(fname))
        start = time.time()
        #newfname = fname.replace("syntheses","errors")
        if os.path.exists(newfname):
            print("Already done: {}".format(newfname))
        session = Session.load(fname)
        session.stellar_parameter_uncertainty_analysis(systematic_errors=[100,0.2,0.2,0.2])
        #session.stellar_parameter_uncertainty_analysis(systematic_errors=[50, 0.15, 0.10, 0.05])
        #eTeff, elogg, evt, eMH = session.stellar_parameters_staterr
        #session.set_stellar_parameters_errors("sys",0,0,0,0)
        #session.set_stellar_parameters_errors("stat",eTeff, elogg, evt, eMH)
        session.compute_all_abundance_uncertainties()
        try:
            session.save(newfname, overwrite=True)
        except OSError:
            pass
        print("Finished {} in {:.1f}s".format(newfname, time.time()-start))
