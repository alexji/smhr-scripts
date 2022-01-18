from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import glob, time, os, sys
import numpy as np
from smh import Session
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)

def update_abundance_table(session, model):
    assert isinstance(model, SpectralSynthesisModel)
    summary_dict = session.summarize_spectral_models(organize_by_element=True)
    for elem in model.metadata["rt_abundances"]:
        try:
            model.metadata["rt_abundances"][elem] = summary_dict[elem][1]
        except KeyError:
            model.metadata["rt_abundances"][elem] = np.nan
    
    # Fill in fixed abundances
    try:
        fitted_result = model.metadata["fitted_result"]
    except KeyError:
        #logger.info("Run at least one fit before setting abundances of "
        #      "fitted element {}!".format(elem))
        pass
    else:
        for i,elem in enumerate(model.elements):
            try:
                abund = summary_dict[elem][1]
            except KeyError:
                #logger.warn("No abundance found for {}, using nan".format(elem))
                abund = np.nan
            key = "log_eps({})".format(elem)
            fitted_result[0][key] = abund
            fitted_result[2]["abundances"][i] = abund

if __name__=="__main__":
    """
    ipython run_syntheses.py file_in.smh file_out.smh
    """
    ## Edit these master list paths to match where your files are
    ## You'll have to edit the master list files themselves to point to the right place too
    master_list_path = "/Users/alexji/Dropbox/S5/linelists/s5_sorted_master_list.txt"
    master_list_extra_path = "/Users/alexji/Dropbox/S5/linelists/s5_sorted_master_list.txt"
    numiter = 3
    
    assert len(sys.argv) == 3, "Format: python run_syntheses.py file_in.smh file_out.smh"
    fname_in = sys.argv[1]
    assert os.path.exists(fname_in)
    assert fname_in.endswith(".smh")

    fname_out = sys.argv[2]
    assert not os.path.exists(fname_out), "{} already exists!".format(fname_out)
    assert fname_out.endswith(".smh")
    
    ## Rename variables to match below
    fname, new_fname = fname_in, fname_out
    
    start = time.time()
    if os.path.exists(new_fname):
        print(new_fname,"already exists! Skipping",fname)
        exit
    print("===========================================")
    print("===========================================")
    print("===========================================")
    print(fname, new_fname)
    print("===========================================")
    print("===========================================")
    print("===========================================")
    
    session = Session.load(fname)
    feh_approx = session.stellar_parameters[3]
    session.import_master_list(master_list_path)
    
    # Iterate fitting all syntheses
    for it in range(numiter):
        print("===========================================")
        print("======Running synthesis iter {}=============".format(it+1))
        print("===========================================")
        for model in session.spectral_models:
            if isinstance(model, ProfileFittingModel): continue
            print("Fitting {} {}".format(model.wavelength, model.species))
            update_abundance_table(session, model)
            try:
                model.iterfit(maxiter=5)
            except:
                print("Failed on this one")
                model.is_acceptable = False
                model.user_flag = True
                continue
            # An approximate upper limit flag
            if model.abundances_to_solar[0] - feh_approx < -3:
                model.is_acceptable = False
                model.user_flag = True
    
    session.import_master_list(master_list_extra_path)
    session.save(new_fname)
    print("===========================================")
    print("===========================================")
    print("===========================================")
    print(fname, new_fname, "{:.1f}".format(time.time()-start))
    print("===========================================")
    print("===========================================")
    print("===========================================")
    #break
