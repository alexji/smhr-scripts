import os, sys, time
from smh import Session
from smh.spectral_models import SpectralSynthesisModel
from smh.specutils import Spectrum1D
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

Nrow, Ncol = 7,5

if __name__=="__main__":
    ## Hack for now, just use the same index
    #imodel = int(sys.argv[1])
    
    tab = ascii.read("stellar_params.txt", delimiter='&')
    tab.sort("Teff")
    #tab = tab[0:10]
    
    stars = np.array(tab["star"])
    smhr_fnames = [f"{star.upper()}.smh" for star in stars]
    to_keep = np.array([os.path.exists(fn) for (star,fn) in zip(stars,smhr_fnames)])
    print(f"Currently have {to_keep.sum()} Stars")
    
    stars = stars[to_keep]
    Teffs, MHs = tab["Teff"][to_keep], tab["MH"][to_keep]
    smhr_fnames = [fn for (fn, k) in zip(smhr_fnames, to_keep) if k]
    start = time.time()
    sessions = [Session.load(fn) for fn in smhr_fnames]
    print(f"Time to read {len(sessions)} smhr files {time.time()-start:.1f}")
    
    Ns = [len(session.spectral_models) for session in sessions]
    for star, N in zip(stars, Ns):
        print(star,N,"lines")
    print(f"(If not all the same number of lines, will do first {N} lines and assume they're sorted")
    
    N = len(stars)
    assert N <= Nrow*Ncol, (N, Nrow, Ncol)
    
    for imodel in range(len(sessions[-1].spectral_models)):
        start = time.time()
        plotmodels = [session.spectral_models[imodel] for session in sessions]
        fig, axes = plt.subplots(Nrow, Ncol, figsize=(Ncol*3,Nrow*3),sharex=True,sharey=True)
        for ax, star, Teff, MH, plotmodel, session in zip(axes.flat, stars, Teffs, MHs, plotmodels, sessions):
            #ii = (spec.dispersion >= w1-1) & (spec.dispersion <= w2+1)
            #ax.plot(spec.dispersion[ii], spec.flux[ii], 'k-')
            spectrum = session.normalized_spectrum
            if isinstance(plotmodel, SpectralSynthesisModel):
                transition1 = plotmodel.transitions[0]
                transition2 = plotmodel.transitions[-1]
                w1, w2 = transition1["wavelength"], transition2["wavelength"]
                species = plotmodel.species[0][0]
                wave = plotmodel.wavelength
            else:
                transition = plotmodel.transitions[0]
                w1, w2 = transition["wavelength"]-1, transition["wavelength"]+1
                species = transition["species"]
                wave = transition["wavelength"]
            ii = (spectrum.dispersion >= w1-0.2) & (spectrum.dispersion <= w2+0.2)
            ax.plot(spectrum.dispersion[ii], spectrum.flux[ii], 'k-', lw=3)
            try:
                named_p_opt, cov, meta = plotmodel.metadata["fitted_result"]
                ax.set(ylim=(0,1.1), xlim=(w1,w2))
                ax.text(.01,.99,f"{star} Teff={Teff:.0f} [Fe/H]={MH:.2f}",ha='left',va='top',transform=ax.transAxes,fontsize=8)
                color = "r" if plotmodel.is_acceptable else "b"
                ax.plot(meta["model_x"], meta["model_y"], '-', lw=2, color=color, alpha=1)
            except KeyError:
                pass
        fig.tight_layout()
        fig.savefig(f"stamp_plots/stampfit_{species:.1f}_{wave:.1f}.png",dpi=200)
        plt.close(fig)
        print(f"stamp_plots/stampfit_{species:.1f}_{wave:.1f}.png took {time.time()-start:.1f}s")
