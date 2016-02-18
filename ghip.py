import numpy as np
import scipy
import scipy.misc
import scipy.integrate
import ipdb

def ghfunc(x,par,ipro=None,offset=None,plot_key=None):

    nx = len(x)

    amp = fan([1.,par[1:10],par[15:18]],nx)
    if par[0] != oldwidth:
        width = par[0]
        p = set_param(x,n,nx,param_set=True)
        oldwidth = par[0]
    #S which axis is this summed over? I believe the first due to idl docs on 
    #S 'total'
    ipro = np.sum(amp*p,axis=0)
    return ipro/scipy.integrate.simps(ipro, x=x)

def make_herm(n,transpose=False):

    """
    Create an array of the first n hermite polynomial coefficients. 
    Returns an array of (n+1)X(n+1) of the coeffs, and can be transposed.
    Emulating code from johnjohn's make_herm.pro

    For reference, the polynomials H_N for order N are
    H_0 = 1
    H_1 = 2x
    H_2 = 4x**2 -2
    H_3 = 8x**3 -12x
    H_4 = 16x**4 - 48x**2 + 12
    """
    #S we need special cases to handle the coefficients less than two, as the
    #S recursion formula works only for n>2. These cases aren't hard though!

    #S make the first element equal to 1
    h = np.zeros([n+1,n+1],dtype=np.float64)
    h[0,0] = 1.
    
    #S if the array is large enough, make element_2,2 equal to 2
    if n > 0:
        h[1,1] = 2.
    #S formula seems to work, found a different one on wikipedia. this one from
    #S make_herm.pro, maybe just the same result? need to work them out to 
    #S equivalence. this returns correct array up to H_10
    if n > 1:
        for ind in range(2,n+1):
            h[ind,:] = np.roll(h[ind-1,:],1)*2. -\
                2.*float(ind-1)*h[ind-2,:]
    #S if we want the transpose
    if transpose:
        return h.T

    #S otherwise just send out the h array
    else:
        return h

#do ghfunc except param, just ignore. make all attributes their own variables
#

def set_param(x,n,nx,param_set=True):
    
    #S number of terms
    n = 15
    #S number of oversampled points
    nx = len(x) #121
    herm_arr = make_herm(n-1,transpose=True)
    if param_set:
        pow = np.arange(n)
        pow_arr = fan(pow, nx)
        xarr = fan(x,n,transpose=True)
        zarr = xarr**powarr
        norm = 1./np.sqrt(2.**pow * np.sqrt(np.pi) * scipy.misc.factorial(pow))
        normarr = fan(norm,nx)
    beta = np.zeros(nx, dtype=np.float64) + width
    betarr = fan(beta, n)**powarr
    gauss = np.exp(-(x * beta)**2/2.)
    gaussarr = fan(gauss,n,transpose=True)
    t = np.dot(herm_arr,zarr * betarr)
    final_p = normarr * gaussarr * t
    
def fan(array,nfan,transpose=False):
    #S function to emulate fan.pro
    temp_list = []
    for i in range(nfan):
        temp_list.append(array)
    #S let's return a numpy array for ease
    fan_array = np.array(temp_list)
    if transpose:
        fan_array = fan_array.T
    return fan_array
    
if __name__ == '__main__':
    ipdb.set_trace()
    z = make_herm(5,transpose=True)
    x = np.arange(121)
    w = set_param(x,15,len(x))
