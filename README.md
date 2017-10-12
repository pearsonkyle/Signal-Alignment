# Signal-Alignment
Algorithms to align 1D signals

These algorithms were used in "Ground-Based Spectroscopy of XO-2b using a Systematic Wavelength Calibration" to correct for temporal wavelength shifts in spectra.

*LINK*
*Image*


I have send two different cross-correlation algorithms in the script that is attached. The first one is labeled "align_spectra" this function shifts spectra by linearly interpolated a reference spectrum to a template within a region of interest, it will then try to minimize the chi-squared between the two spectra as a function of the shift amount. Each spectra is normalized prior to shifting so that the algorithm focuses maining on comparing the shape of the data. 

The second alignment routine is called "phase_spectra" this will perform a phase correlation in fourier space (https://en.wikipedia.org/wiki/Phase_correlation). It applies the same technique as that link except in 1D instead of 2D. This algorithm should do the same thing as the chi-squared alignment technique but perhaps at a lower precision. The only caveat about the phase correlation is if you cross-correlate your data at the native resolution you can only achieve integer precision. I've added an option to upscale the data using the "res" parameter in the function. res will upscale the data N times so that you can achieve a precision of 1/N as far as the shift value goes. That means if N=100, then the smallest shift value you can get is 0.01. 

I have also included a function to fit a gaussian plus a line that will derive the center of the gaussian. This function is called fit_psf .

I've also included a small example at the bottom of the script for how to use most of these functions. 
