from pathlib import Path
from astropy.table import Table
import numpy as np
import scipy.interpolate as spinter
import matplotlib.pyplot as plt

#for now I am locally defining this across the board but I should make this its own thing in a file.io function library
def searcher(start_path, dirname):
    # Start from a high-level but not-too-huge root directory
    search_root = Path(str(start_path))

    # Recursively look for the directory
    matches = list(search_root.rglob(str(dirname)))

    if matches:
        jwst_dir = matches[0]  # Use the first match (or loop over all if multiple found)
        print(f"Found directory: {jwst_dir}")

        # List all files in its subdirectories
        all_files = [f for f in jwst_dir.glob("**/*") if f.is_file()]
        print(f"Found {len(all_files)} files:")
        disp_files = [str(f) for f in all_files if f.suffix == ".fits"]
        file_dict = {}
        for f in disp_files:
            df = Table.read(f)
            dex = (f.split("/")[-1]).split("_")[2]
            file_dict[dex] = df
        return file_dict
    else:
        print(f"No directory named {dirname} found.")



#Uncomment to see disperser options
#print(JWST_disp_dict.keys())

#Uncomment to view specific fits files


###TODO: Add the corresponding data format for the data files so that we can instantiate the disperser class once for each disperser, 
# then call the binner function for multiple spectra.

class JWST_disperser:
    def __init__(self, disperser_dict): 
        self.respow = np.array(disperser_dict["R"])
        self.dlds = np.array(disperser_dict["DLDS"])
        self.wavelength = np.array(disperser_dict["WAVELENGTH"])
        self.minwave = min(self.wavelength)
        self.maxwave = max(self.wavelength)
    
    def edge_finder(self, wavelengths = None, bin_widths = None, respower = None, start = None, end = None): # I should probably be asking for a list here as only taking part of this is probably not accurate
        # Default to instance values if no arguments provided
        if wavelengths is None:
            wavelengths = self.wavelength
        if bin_widths is None:
            bin_widths = self.dlds
        if respower is None:
            respower = self.respow
        if start is None:
            start = self.minwave
        if end is None:
            end = self.maxwave
        
        # Interpolate bin width as a function of wavelength
        interp_respow = spinter.CubicSpline(wavelengths, respower) #if I need this outside, need to assign as an attribute
        interp_bin_width = spinter.CubicSpline(wavelengths, bin_widths * 2.2)
        # Create adaptive bins
        bin_edges = [start]
        current = start
        while current < end:
            width = current / interp_respow(current)
            next_edge = current + width
            if next_edge > end:
                break
            bin_edges.append(next_edge)
            current = next_edge
        # Ensure last bin edge reaches the end
        if bin_edges[-1] < end:
            bin_edges.append(end)
        bincents = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
        xerrbars = np.diff(np.array(bin_edges)) / 2
        self.bincs = bincents
        self.xerbs = xerrbars
        # Convert bin edges to numpy array
        bin_edges = np.array(bin_edges)

        self.bin_edges = bin_edges
        self.interp_bin_width = interp_bin_width
        return bin_edges, bincents, xerrbars, interp_bin_width
    
    def binner(self, lcents, sigmas, lams, fluxes, bin_edges = None): 
        if bin_edges is None:
            bin_edges = self.bin_edges
        ### The below will need to be appended to get arrays of these values for each individual spectrum
        self.lcents = lcents
        self.sigmas = sigmas
        #self.lams = lams #the wavelength array (x values)
        #self.fluxes = fluxes #y-values
        #TODO: Breakdown
        
        ### Goal of this loop is to 
        if len(sigmas) == 1:
            sigmas = np.ones(len(lcents)) * sigmas
        elif len(sigmas) == len(lcents):
            return
        else:
            print(f"!!!Length mismatch!!! \nLine centers length: {len(lcents)} \nSigmas length: {len((sigmas))} ")

        #set limits for valid bin ranges
        lolims = lcents - 3 * sigmas
        uplims = lcents + 3 * sigmas
        for i in range(lcents):
            #is bin edge within +- 3 sigma of lcent
            boolmsk = lolims[i] <= bin_edges <= uplims[i]
            

            print("not done")
            #loop over lcents (lcents and sigmas should be same dimensional array, len(sigmas) = 1 or len(lcents))
        
        #TODO
        #integrate overfluxes in each bin(how to do coherently without getting binning artifacts for small bins?)
        #make flux plot (histogram?, scatter? line? je ne sais quoi!)
        return 

    def pix_scale(self, wavelengths = None, bin_widths = None, respower = None): #see edge finder comment
        if wavelengths is None:
            wavelengths = self.wavelength
        if bin_widths is None:
            bin_widths = self.dlds
        if respower is None:
            respower = self.respow

        
        FWHM_l = wavelengths/(respower)
        dl_dpx = FWHM_l / 2.2
        return dl_dpx
        
        
#TODO: Incorporate dsiperser output into filter to obtain actual filtered spectrum.
class JWST_filter:
    def __init__(self, filter_dict):
        self.throughput = np.array(filter_dict["THROUGHPUT"])
        self.wavelength = np.array(filter_dict["WAVELENGTH"])
        self.minwave = min(self.wavelength)
        self.maxwave = max(self.wavelength)
        
    
    def inter(self):
        smth_thr = spinter.interp1d(self.wavelength, self.throughput, kind = 'cubic')
        return smth_thr

#the first instantiation of this class as an example
'''
f110w = JWST_filter(JWST_fil_dict["f110w"])
interpolated_test = JWST_filter.inter(f110w)

g140h = JWST_disperser(JWST_disp_dict["g140h"])

g395h = JWST_disperser(JWST_disp_dict["g395h"])

'''


if __name__ == "__main__":

    JWST_disp_dict = searcher("/Users/lamoreau/python/ASpec", "NIRSpecdis")
    JWST_fil_dict = searcher("/Users/lamoreau/python/ASpec", "NIRSpecfil")

    #print(JWST_disp_dict["g395h"])

    test_disp = JWST_disperser(JWST_disp_dict["prism"])

    #print(test_disp.pix_scale())
    #print(test_disp.dlds)
    #we are officially twinning!!!!




    import numpy as np
    from scipy import signal
    from scipy.interpolate import interp1d

    # wavelength grid that matches your dispersion table
    wave_hr = test_disp.wavelength       # make sure your class has this
    dlds    = test_disp.dlds       # Δλ per pixel (micron/pixel)
    R_hr    = test_disp.respow     # resolving power array

    # --- Construct high-res input spectrum ---
    flux_hr = 0.2 + 0.02*np.sin(10*wave_hr)  # baseline continuum wiggles

    # add random emission lines
    rng = np.random.default_rng(42)
    line_positions = rng.uniform(0.7, 5.0, 20)
    for lam0 in line_positions:
        flux_hr += 1.0 * np.exp(-0.5*((wave_hr-lam0)/1e-4)**2)

    # FWHM and sigma for LSF from R(λ)
    FWHM_hr  = wave_hr / R_hr
    sigma_hr = FWHM_hr / (2 * np.sqrt(2 * np.log(2)))

    def convolve_locally(wave, flux, sigma_lambda):
        """Gaussian convolution with wavelength-dependent sigma."""
        n = len(wave)
        out = np.zeros_like(flux)
        chunk = 2000
        for i0 in range(0, n, chunk):
            i1 = min(n, i0 + chunk)
            ww = wave[i0:i1]
            ff = flux[i0:i1]
            sig = np.median(sigma_lambda[i0:i1])
            dl  = np.mean(np.diff(ww))
            half_npix = int(np.ceil(8 * sig / dl))
            kx = np.arange(-half_npix, half_npix + 1) * dl
            kernel = np.exp(-0.5 * (kx / sig) ** 2)
            kernel /= np.sum(kernel)
            padded = np.pad(ff, half_npix, mode='constant', constant_values=0.0)
            conv   = signal.fftconvolve(padded, kernel, mode='same')
            out[i0:i1] = conv[half_npix:-half_npix]
        return out

    conv_flux_hr = convolve_locally(wave_hr, flux_hr, sigma_hr)

    # Build pixel grid directly from dlds (flux-conserving bins)
    pix_edges = [wave_hr[0]]
    lam = wave_hr[0]
    while lam < wave_hr[-1]:
        idx = np.searchsorted(wave_hr, lam)
        if idx >= len(dlds):
            break
        step = dlds[idx]
        lam += step
        pix_edges.append(lam)
    pix_edges = np.array(pix_edges)
    pix_centers = 0.5 * (pix_edges[:-1] + pix_edges[1:])

    # Bin convolved spectrum into pixels
    interp_flux = interp1d(wave_hr, conv_flux_hr, bounds_error=False, fill_value=0.0)
    binned_flux = np.array([
        np.trapz(interp_flux([pix_edges[i], pix_edges[i+1]]),
                 [pix_edges[i], pix_edges[i+1]]) /
        (pix_edges[i+1] - pix_edges[i])
        for i in range(len(pix_centers))
    ])

    print(f"Simulated {len(pix_centers)} NIRSpec pixels "
          f"from {pix_centers[0]:.2f} to {pix_centers[-1]:.2f} μm")
    
    # --- Build pixel grid from dlds ---
    pix_edges = [wave_hr[0]]
    lam = wave_hr[0]
    while lam < wave_hr[-1]:
        idx = np.searchsorted(wave_hr, lam)
        if idx >= len(test_disp.dlds):
            break
        step = test_disp.dlds[idx]
        lam += step
        pix_edges.append(lam)
    pix_edges = np.array(pix_edges)
    pix_centers = 0.5*(pix_edges[:-1] + pix_edges[1:])

    # Bin convolved flux into pixels
    interp_flux = interp1d(wave_hr, conv_flux_hr, bounds_error=False, fill_value=0.0)
    binned_flux = np.array([
        np.trapz(interp_flux([pix_edges[i], pix_edges[i+1]]),
                [pix_edges[i], pix_edges[i+1]]) /
        (pix_edges[i+1]-pix_edges[i])
        for i in range(len(pix_centers))
    ])

    # --- Plot ---
    plt.figure(figsize=(10,5))
    plt.plot(wave_hr, flux_hr, color="gray", alpha=0.6, label="High-res input")
    plt.plot(wave_hr, conv_flux_hr, color="blue", lw=1, label="After LSF convolution")
    plt.step(pix_centers, binned_flux, where="mid", color="red", label="Binned to NIRSpec pixels")
    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Flux (arb. units)")
    plt.title("JWST NIRSpec PRISM Resolution & Binning")
    plt.legend()
    plt.tight_layout()
    plt.show()