import colorsys
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift
from statsmodels.tsa.stattools import ccovf, ccf
from scipy import signal
import matplotlib.pyplot as plt

def align_spectra(reference, target, ROI, order=1,init=0.1,res=1,b=1):
    '''
        NH[0], NH[i]

        Aligns the target spectrum with in the region of interest (ROI) to the reference spectrum's ROI

        res - resolution of the data, only used if passing in higher resolution data and the initial value
        is given in native pixel coordinates not the high res coordinates

        b - symmetric bounds for constraining the shift search around the initial guess
          -
    '''
    ROI[0] = int(ROI[0]*res)
    ROI[1] = int(ROI[1]*res)

    # ROI - region of interest to focus on computing the residuals for
    # LIMS - shifting limits
    reference = reference/np.mean(reference[ROI[0]:ROI[1]])

    # define objective function: returns the array to be minimized
    def fcn2min(x):
        # x = shift length
        shifted = shift(target,x,order=order)
        shifted = shifted/np.mean(shifted[ROI[0]:ROI[1]])
        return np.sum( ((reference - shifted)**2 )[ROI[0]:ROI[1]] )


    #result = minimize(fcn2min,init,method='Nelder-Mead')
    minb = min( [(init-b)*res,(init+b)*res] )
    maxb = max( [(init-b)*res,(init+b)*res] )
    result = minimize(fcn2min,init,method='L-BFGS-B',bounds=[ (minb,maxb) ])
    return result.x[0]/res

# add a nested sampler
# try alignment with just |chi| (doesn't seem to make much difference)



def phase_spectra(ref,tar,ROI,res=100):
    '''
        Cross-Correlate data within ROI with a precision of 1./res
        interpolate data onto higher resolution grid and
        align target to reference
    '''
    x,r1 = highres(ref[ROI[0]:ROI[1]],kind='linear',res=res)
    x,r2 = highres(tar[ROI[0]:ROI[1]],kind='linear',res=res)

    r1 -= r1.mean()
    r2 -= r2.mean()

    cc = ccovf(r1,r2,demean=False,unbiased=False)
    if np.argmax(cc) == 0:
        cc = ccovf(r2,r1,demean=False,unbiased=False)
        mod = -1
    else:
        mod = 1

    s1 = np.argmax(cc)*mod*(1./res)
    return s1

    x,r1 = highres(ref[ROI[0]:ROI[1]],kind='linear',res=res)
    x,r2 = highres(tar[ROI[0]:ROI[1]],kind='linear',res=res)

    r1 -= r1.mean()
    r1 -= r2.mean()

    # compute the POC function
    product = np.fft.fft(r1) * np.fft.fft(r2).conj()
    cc = np.fft.fftshift(np.fft.ifft(product))

    l = ref[ROI[0]:ROI[1]].shape[0]
    shifts = np.linspace(-0.5*l,0.5*l,l*res)

    return shifts[np.argmax(cc.real)]

    # applying filter removes spectral info at the edge of bin
    #f1 = np.fft.fft(r1)
    #f1[20:-20] = 0
    #df = np.fft.ifft(f1)


def highres(y,kind='cubic',res=100):
    # interpolate onto higher resolution grid with res* more data points than original input
    # from scipy import interpolate
    y = np.array(y)
    x = np.arange(0, y.shape[0])
    f = interp1d(x, y,kind='cubic')
    xnew = np.linspace(0, x.shape[0]-1, x.shape[0]*res)
    ynew = f(xnew)
    return xnew,ynew

def error(x,y):
    # basic uncertainty on poisson quantities of x and y for f(x,y) = x/y
    sigx = np.sqrt(x)
    sigy = np.sqrt(y)
    dfdx = 1./y
    dfdy = x/(y*y)

    er = np.sqrt( dfdx**2 * sigx**2 + dfdy**2 * sigy**2 )
    return er


def fit_psf(flux, ci=0, npts=10, debug=False):
    # ci is the initial guess to gaussian center 

    # pull out flux region around spectral line
    Y = flux[int(np.sign(ci)*ci-0.5*npts):int(np.sign(ci)*ci+0.5*npts)]
    X = np.arange(int(np.sign(ci)*ci-0.5*npts),int(np.sign(ci)*ci+0.5*npts))

    # just find the lowest point in the spectrum and center around there
    if ci > 0: # fit trough
        cp = np.argmin(Y)
        ci = X[cp]
    else: # fit peak
        cp = np.argmax(Y)
        ci = X[cp]

    # pull out flux region around spectral line
    Y = flux[int(ci-0.5*npts):int(ci+0.5*npts)]
    X = np.arange(int(ci-0.5*npts),int(ci+0.5*npts))


    # assume all spectral lines have negative amplitude
    ai = np.min(Y) - np.max(Y)
    bi = np.max(Y)

    line = lambda x,m,b : m*x + b
    def star_psf(x,a,x0,sigx,m,b):
        gaus = a * np.exp(-(x-x0)**2 / (2*sigx**2) )
        gaus += line(x,m,b)
        return gaus

    def fcn2min(p):
        return np.sum(  (Y - star_psf(X,*p))**2 )

    #return X[np.argmin(Y)

    init = [ai,ci,2*(npts/10.),0,bi]
    result = minimize(fcn2min,init,method='Nelder-Mead')

    if debug:
        import pdb; pdb.set_trace()

    return result.x[1]



if __name__ == "__main__":

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

    plt.plot(x,y,label='original data')
    plt.plot(x,shifted,label='shifted data')
    plt.plot(x,shift(shifted,s),label='aligned data') # use shift function to linearly interp data
    plt.legend(loc='best')
    plt.show()
