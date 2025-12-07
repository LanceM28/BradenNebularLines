from pathlib import Path
from astropy.table import Table
import numpy as np
import scipy.interpolate as spinter
import matplotlib.pyplot as plt
import os
from scipy import signal
from scipy.interpolate import interp1d
import astropy.units as u

#ham
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



figstore = "/Users/lamoreau/Documents/ASpecfigs/"
os.makedirs(figstore[0:-1], exist_ok=True)
#Uncomment to see disperser options
JWST_disp_dict = searcher("/Users/lamoreau/python/ASpec", "NIRSpecdis")

#Uncomment to view specific fits files


###TODO: Add the corresponding data format for the data files so that we can instantiate the disperser class once for each disperser, 
# then call the binner function for multiple spectra.

class JWST_disperser:
    JWST_disp_dict = searcher("/Users/lamoreau/python/ASpec", "NIRSpecdis")
    JWST_fil_dict = searcher("/Users/lamoreau/python/ASpec", "NIRSpecfil")
    def __init__(self, disp_name, disperser_dict = JWST_disp_dict):
        self.disp_name = disp_name 
        self.respow = np.array(disperser_dict[disp_name]["R"])
        self.dlds = np.array(disperser_dict[disp_name]["DLDS"])
        self.wavelength = np.array(disperser_dict[disp_name]["WAVELENGTH"])
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
    
    def binner(self, flux_in = None, wave_in = None): 
        dispname = self.disp_name #maybe change this for consistent naming at some point
        # --- Input or generate flux ---
        if flux_in is None:
            # synthetic spectrum if none provided
            wave_hr = np.linspace(self.minwave,self.maxwave,1000) * (u.micron)
            flux_hr = (0.2 + 0.02*np.sin(10*wave_hr.value)) * (u.erg / (u.s * u.cm**2 * u.micron))
            rng = np.random.default_rng(42)
            for lam0 in rng.uniform(self.minwave, self.maxwave, 20):
                flux_hr += 1.0 * np.exp(-0.5*((wave_hr.value - lam0)/1e-4)**2) * (u.erg / (u.s * u.cm**2 * u.micron))
        else:
            wave_hr = wave_in * u.microns
            flux_hr = flux_in * (u.erg / (u.s * u.cm**2 * u.micron))

        # FWHM and sigma for LSF from R(λ)
        FWHM_hr  = self.wavelength / self.respow * (u.micron)#DINGO: FWHM from wavelength and R, where does this equaition come from
        sigma_hr = FWHM_hr / (2 * np.sqrt(2 * np.log(2))) #DINGO: might as well look up this one too
        #frogged up to here
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
                half_npix = int(np.ceil(8 * sig / dl).to_value())
                kx = np.arange(-half_npix, half_npix + 1) * dl
                kernel = np.exp(-0.5 * (kx / sig).to_value() ** 2)
                kernel /= np.sum(kernel)
                ff_val = ff.to_value()
                padded = np.pad(ff_val, half_npix, mode='edge')
                conv   = signal.fftconvolve(padded, kernel, mode='same')
                conv_q = conv * flux.unit
                out[i0:i1] = conv_q[half_npix:-half_npix]
            return out
        
        # --- Extend wavelength and flux before convolution to avoid edge losses ---
        # JWST LSF convolution smears flux across a few resolution elements.
        # If we stop the wavelength grid abruptly, the convolution kernel "falls off"
        # the edge and dilutes the first/last bins. To avoid this, we pad the
        # spectrum with a flat continuum at both ends before convolution, then trim.

        pad_width = 5  # in units of resolution elements (FWHM) to extend
        dl = np.median(np.diff(wave_hr)).to(u.micron)

        # how far to pad (pick ~5 FWHM at each end)
        pad_left  = pad_width * np.median(FWHM_hr[:100])
        pad_right = pad_width * np.median(FWHM_hr[-100:])#problems if len < 100
        
        print(wave_hr[0].to_value(), pad_left.to_value(), dl.to_value())
        # build extended wavelength grid
        wave_ext = np.concatenate([
            np.arange(wave_hr[0].to_value()-pad_left.to_value(), wave_hr[0].to_value(), dl.to_value()),
            wave_hr.to_value(),
            np.arange(wave_hr[-1].to_value()+dl.to_value(), wave_hr[-1].to_value()+pad_right.to_value(), dl.to_value())
        ]) #what is wave_ext and wave_hr

        # extend flux with flat extrapolation (baseline = edge values)
        flux_ext = np.concatenate([
            np.full(len(wave_ext) - len(wave_hr) - (len(wave_ext)-len(wave_hr))//2, flux_hr[0]) * flux_hr.unit,
            flux_hr,
            np.full((len(wave_ext)-len(wave_hr))//2, flux_hr[-1]) * flux_hr.unit
        ])
        #print(f"Flux_hr = {flux_hr} \nFlux_ext = {flux_ext}")
        # extend sigma array (approximate edges with nearest values)
        sigma_ext = np.concatenate([
            np.full(len(wave_ext) - len(wave_hr) - (len(wave_ext)-len(wave_hr))//2, sigma_hr[0]) * sigma_hr.unit,
            sigma_hr,
            np.full((len(wave_ext)-len(wave_hr))//2, sigma_hr[-1]) * sigma_hr.unit
        ])

        # now convolve on the extended grid
        conv_flux_ext = convolve_locally(wave_ext, flux_ext, sigma_ext)

        # trim back to original wavelength coverage
        mask = (wave_ext >= wave_hr[0].to_value()) & (wave_ext <= wave_hr[-1].to_value())
        conv_flux_hr = conv_flux_ext[mask]
        

        conv_flux_hr = convolve_locally(wave_hr, flux_hr, sigma_hr)

        # Build pixel grid directly from dlds (flux-conserving bins)
        pix_edges = [wave_hr[0].to_value()]
        lam = wave_hr[0] #in microns
        while lam < wave_hr[-1]:
            idx = np.searchsorted(wave_hr, lam)
            if idx >= len(self.dlds):
                break
            step = self.dlds[idx] * wave_hr.unit
            lam += step
            pix_edges.append(lam.to_value())
        pix_edges = np.array(pix_edges) *wave_hr.unit
        pix_edges[-1] = wave_hr[-1]
        pix_centers = 0.5 * (pix_edges[:-1] + pix_edges[1:])

        # Bin convolved spectrum into pixels
        interp_flux = interp1d(wave_hr, conv_flux_hr, bounds_error=False, fill_value=0.0)
        
        ######
        # trying to utelize binning that is based on the subgrid scale that I already have.
        
        # (1) sanity checks #do I need this and what is the slowdown?
        if not np.all(np.diff(wave_hr) > 0):
            raise ValueError("wave_hr must be strictly increasing")

        # (2) build cumulative integral using trapezoidal rule
        # area of each small segment between wave_hr[i] and wave_hr[i+1]
        dl = np.diff(wave_hr)                               # length N-1
        seg_area = 0.5 * (conv_flux_hr[:-1] + conv_flux_hr[1:]) * dl  # length N-1
        cumI = np.concatenate(([0.0], np.cumsum(seg_area))) # length N ; cumI[i] = integral from wave_hr[0] to wave_hr[i]

        # (3) make an interpolator for the cumulative integral vs wavelength
        cumI_interp = interp1d(wave_hr, cumI, kind='linear',
                            bounds_error=False,
                            fill_value=(0.0, cumI[-1]), assume_sorted=True)

        # (4) evaluate cumulative integral at pixel edges, then take differences
        I_edges = cumI_interp(pix_edges)            # length M+1 if pix_edges length M+1
        pixel_areas = np.diff(I_edges)              # integral within each pixel
        pixel_widths = np.diff(pix_edges)
        binned_flux = pixel_areas / pixel_widths    # average flux density inside each pixel

        #######

        # binned_flux corresponds to pixels with centers
        pix_centers = 0.5*(pix_edges[:-1] + pix_edges[1:])

        print(f"Simulated {len(pix_centers)} NIRSpec pixels "
            f"from {pix_centers[0]:.2f} to {pix_centers[-1]:.2f} μm")

        # --- Plot ---
        plt.figure(figsize=(10,5))
        plt.plot(wave_hr, flux_hr, color="gray", alpha=0.6, label="High-res input")
        plt.plot(wave_hr, conv_flux_hr, color="blue", lw=1, label="After LSF convolution")
        plt.step(pix_centers, binned_flux, where="mid", color="red", label="Binned to NIRSpec pixels")
        plt.xlabel("Wavelength (μm)")
        plt.ylabel("Flux (arb. units)")
        plt.title(f"JWST NIRSpec Resolution & Binning ({dispname})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{figstore}JWST NIRSpec Resolution & Binning ({dispname})", transparent = True)
        
        #TODO
        #integrate overfluxes in each bin(how to do coherently without getting binning artifacts for small bins?)
        #make flux plot (histogram?, scatter? line? je ne sais quoi!)
        return 

if __name__ == "__main__":

    g395h = JWST_disperser("g395m")
    g395h.binner()

