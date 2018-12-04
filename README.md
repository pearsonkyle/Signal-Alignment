# Signal-Alignment
Algorithms to align 1D signals by deriving the offset using different cross-correlation methods. 


## Alignment via chi-squared minimization
The first alignment technique is labeled "align_spectra". This function shifts a signal by resampling it via linear interpolation, the edges of the signal are ignored because the algorithm focuses on aligning two signals within a region of interest. The shifted signal is then compared to a "template" signal where a chi-squared minimization occurs to derive the optimal offset between the two. Each signal (template and shifted) is normalized (mean=1) prior to shifting so that the algorithm focuses mainly on comparing the shape of the data. 


## Alignment via Fourier space
The second alignment routine is called "phase_spectra" and uses a phase correlation in fourier space (https://en.wikipedia.org/wiki/Phase_correlation). It applies the same technique as that link except in 1D instead of 2D. After independent testing this algorithm behaves in a similar manner as the chi-squared minimization except at a lower precision. The only caveat about the phase correlation is if you cross-correlate your data at the native resolution you can only achieve integer precision. To achieve sub-pixel/ sub-integer precision, I've added an option to upscale the data prior to the cross-correlation. By upscaling the data N times you can achieve a precision of 1/N as far as the shift value goes. That means if N=100, then the smallest shift value you can get is 0.01. 

![](https://github.com/pearsonkyle/Signal-Alignment/blob/master/images/cross-correlation.png)
Here is an image of the alignment algorithm in use for real data. I was interested in comparing the shift between two signals that happen to be astrophysical spectra at two different wavelength regions. 

## Example
Using each algorithm is fairly straight forward, pass each function two 1D arrays as well as bounds 
``` python
    NPTS = 100
    SHIFTVAL = 4
    NOISE = 1e-3

    # generate some noisy data and simulate a shift
    x = np.linspace(0,4*np.pi,NPTS)
    y = signal.gaussian(NPTS, std=4) * np.random.normal(1,NOISE,NPTS)
    shifted = np.roll( signal.gaussian(NPTS, std=4) ,SHIFTVAL) * np.random.normal(1,NOISE,NPTS)
    # np roll can only do integer shifts

    # align the shifted spectrum back to the real
    s = phase_spectra(y, shifted, [10,190])
    print('phase shift value to align is',s)

    # chi squared alignment at native resolution
    s = align_spectra(y, shifted, [10,190],init=-4,b=1)
    print('chi square alignment',s)

    plt.plot(x,y,'k-',label='original data')
    plt.plot(x,shifted,'r-',label='shifted data')
    plt.plot(x,shift(shifted,s),'o--',label='aligned data') # use shift function to linearly interp data
    plt.legend(loc='best')
    plt.show()
```

## Use cases
* These algorithms are used in [Ground-Based Spectroscopy of XO-2b using a Systematic Wavelength Calibration](https://www.overleaf.com/read/nprsyzqjtqqb) to correct for time-dependent wavelength shifts as a result of atmospheric differential refraction. 

* Additionally, these techniques are used as part of the on-chip guiding system for the [NESSI](https://en.wikipedia.org/wiki/New_Mexico_Exoplanet_Spectroscopic_Survey_Instrument) instrument at Palomar Observatory (See [SPIE proceeding: The New NESSI – Refurbishment of an NIR MOS for characterizing exoplanets using the Hale Telescope](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/10702/107023K/The-new-NESSI--refurbishment-of-a-NIR-MOS-for/10.1117/12.2314242.short?SSO=1))
