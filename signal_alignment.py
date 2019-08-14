#
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift
from statsmodels.tsa.stattools import ccovf

def equalize_array_size(array1,array2):
    '''
    reduce the size of one sample to make them equal size. 
    The sides of the biggest signal are truncated

    Args:
        array1 (1d array/list): signal for example the reference
        array2 (1d array/list): signal for example the target

    Returns:
        array1 (1d array/list): middle of the signal if truncated
        array2 (1d array/list): middle of the initial signal if there is a size difference between the array 1 and 2
        dif_length (int): size diffence between the two original arrays 
    '''
    len1, len2 = len(array1), len(array2)
    dif_length = len1-len2
    if dif_length<0:
        array2 = array2[int(np.floor(-dif_length/2)):len2-int(np.ceil(-dif_length/2))]
    elif dif_length>0:
        array1 = array1[int(np.floor(dif_length/2)):len1-int(np.ceil(dif_length/2))]
    return array1,array2, dif_length

def chisqr_align(reference, target, roi=None, order=1, init=0.1, bound=1):
    '''
    Align a target signal to a reference signal within a region of interest (ROI)
    by minimizing the chi-squared between the two signals. Depending on the shape
    of your signals providing a highly constrained prior is necessary when using a
    gradient based optimization technique in order to avoid local solutions.

    Args:
        reference (1d array/list): signal that won't be shifted
        target (1d array/list): signal to be shifted to reference
        roi (tuple): region of interest to compute chi-squared
        order (int): order of spline interpolation for shifting target signal
        init (int):  initial guess to offset between the two signals
        bound (int): symmetric bounds for constraining the shift search around initial guess

    Returns:
        shift (float): offset between target and reference signal 
    
    Todo:
        * include uncertainties on spectra
        * update chi-squared metric for uncertainties
        * include loss function on chi-sqr

    '''
    reference, target, dif_length = equalize_array_size(reference,target)
    if roi==None: roi = [0,len(reference)-1] 
  
    # convert to int to avoid indexing issues
    ROI = slice(int(roi[0]), int(roi[1]), 1)

    # normalize ref within ROI
    reference = reference/np.mean(reference[ROI])

    # define objective function: returns the array to be minimized
    def fcn2min(x):
        shifted = shift(target,x,order=order)
        shifted = shifted/np.mean(shifted[ROI])
        return np.sum( ((reference - shifted)**2 )[ROI] )

    # set up bounds for pos/neg shifts
    minb = min( [(init-bound),(init+bound)] )
    maxb = max( [(init-bound),(init+bound)] )

    # minimize chi-squared between the two signals 
    result = minimize(fcn2min,init,method='L-BFGS-B',bounds=[ (minb,maxb) ])

    return result.x[0] + int(np.floor(dif_length/2))


def phase_align(reference, target, roi, res=100):
    '''
    Cross-correlate data within region of interest at a precision of 1./res
    if data is cross-correlated at native resolution (i.e. res=1) this function
    can only achieve integer precision 

    Args:
        reference (1d array/list): signal that won't be shifted
        target (1d array/list): signal to be shifted to reference
        roi (tuple): region of interest to compute chi-squared
        res (int): factor to increase resolution of data via linear interpolation
    
    Returns:
        shift (float): offset between target and reference signal 
    '''
    # convert to int to avoid indexing issues
    ROI = slice(int(roi[0]), int(roi[1]), 1)

    # interpolate data onto a higher resolution grid 
    x,r1 = highres(reference[ROI],kind='linear',res=res)
    x,r2 = highres(target[ROI],kind='linear',res=res)

    # subtract mean
    r1 -= r1.mean()
    r2 -= r2.mean()

    # compute cross covariance 
    cc = ccovf(r1,r2,demean=False,unbiased=False)

    # determine if shift if positive/negative 
    if np.argmax(cc) == 0:
        cc = ccovf(r2,r1,demean=False,unbiased=False)
        mod = -1
    else:
        mod = 1

    # often found this method to be more accurate then the way below
    return np.argmax(cc)*mod*(1./res)

    # interpolate data onto a higher resolution grid 
    x,r1 = highres(reference[ROI],kind='linear',res=res)
    x,r2 = highres(target[ROI],kind='linear',res=res)

    # subtract off mean 
    r1 -= r1.mean()
    r1 -= r2.mean()

    # compute the phase-only correlation function
    product = np.fft.fft(r1) * np.fft.fft(r2).conj()
    cc = np.fft.fftshift(np.fft.ifft(product))

    # manipulate the output from np.fft
    l = reference[ROI].shape[0]
    shifts = np.linspace(-0.5*l,0.5*l,l*res)

    # plt.plot(shifts,cc,'k-'); plt.show()
    return shifts[np.argmax(cc.real)]


def highres(y,kind='cubic',res=100):
    '''
    Interpolate data onto a higher resolution grid by a factor of *res*

    Args:
        y (1d array/list): signal to be interpolated
        kind (str): order of interpolation (see docs for scipy.interpolate.interp1d)
        res (int): factor to increase resolution of data via linear interpolation
    
    Returns:
        shift (float): offset between target and reference signal 
    '''
    y = np.array(y)
    x = np.arange(0, y.shape[0])
    f = interp1d(x, y,kind='cubic')
    xnew = np.linspace(0, x.shape[0]-1, x.shape[0]*res)
    ynew = f(xnew)
    return xnew,ynew


if __name__ == "__main__":
    from scipy import signal
    import matplotlib.pyplot as plt

    NPTS = 100
    SHIFTVAL = 4
    NOISE = 1e-2 # can perturb offset retrieval from true
    print('true signal offset:',SHIFTVAL)

    # generate some noisy data and simulate a shift
    y = signal.gaussian(NPTS, std=4) + np.random.normal(1,NOISE,NPTS)
    shifted = shift( signal.gaussian(NPTS, std=4) ,SHIFTVAL) + np.random.normal(1,NOISE,NPTS)

    # align the shifted spectrum back to the real
    s = phase_align(y, shifted, [10,90])
    print('phase shift value to align is',s)

    # chi squared alignment at native resolution
    s = chisqr_align(y, shifted, [10,90], init=-3.5,bound=2)
    print('chi square alignment',s)

    # make some diagnostic plots
    plt.plot(y,label='original data')
    plt.plot(shifted,label='shifted data')
    plt.plot(shift(shifted,s,mode='nearest'),ls='--',label='aligned data') 
    plt.legend(loc='best')
    plt.show()
