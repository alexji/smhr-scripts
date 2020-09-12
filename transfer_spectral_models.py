import numpy as np
from smh import Session
from smh.spectral_models import ProfileFittingModel, SpectralSynthesisModel

from optparse import OptionParser
import sys, os

if __name__=="__main__":
    parser = OptionParser(usage="%prog file_with_models.smh file_needing_models.smh output_file.smh")
    parser.add_option("--force", dest="force", action="store_true", default=False)
    parser.add_option("--refit", dest="refit", action="store_false", default=True)
    options, args = parser.parse_args()
    
    origfname, datafname, outfname = args

    tmpfname = "./.tmp.pkl"
    while os.path.exists(tmpfname):
        tmpfname = tmpfname[:-4]+"p.pkl"
    
    assert os.path.exists(origfname), "input 1 {} does not exist".format(origfname)
    assert os.path.exists(datafname), "input 2 {} does not exist".format(datafname)
    assert options.force or (not os.path.exists(outfname)), "output {} already exists".format(outfname)
    
    session1 = Session.load(origfname)
    session2 = Session.load(datafname)
    
    profile_keys = ["profile","central_weighting","window","continuum_order","detection_sigma",
                    "detection_pixels","max_iterations","velocity_tolerance","mask","antimask_flag",
                    "elements","species"]
    synthesis_keys = ["mask","window","continuum_order","velocity_tolerance","smoothing_kernel",
                      "initial_abundance_bounds","elements","species",
                      "manual_continuum","manual_sigma_smooth","manual_rv"]
    
    new_spectral_models = []
    for model in session1.spectral_models:
        if isinstance(model, ProfileFittingModel):
            newmodel = ProfileFittingModel(session2, model.transitions)
            for key in profile_keys:
                newmodel.metadata[key] = model.metadata[key]
            newmodel._update_parameter_names()
            newmodel._verify_transitions()
            newmodel._verify_metadata()
            if options.refit and model.is_acceptable:
                newmodel.fit()
        elif instance(model, SpectralSynthesisModel):
            newmodel = ProfileFittingModel(session2, model.transitions, model.metadata["elements"],
                                           what_species=model.metadata["species"][0],
                                           what_wavelength=model.wavelength,
                                           what_expot=model.expot,
                                           what_loggf=model.loggf)
            for key in synthesis_keys:
                newmodel.metadata[key] = model.metadata[key]
            newmodel._update_parameter_names()
            newmodel._verify_transitions()
            pass
        new_spectral_models.append(newmodel)
    session2.metadata["spectral_models"] = new_spectral_models
    session2.save(outfname, overwrite=options.force)
