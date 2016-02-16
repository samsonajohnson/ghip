import numpy as np
import scipy
import scipy.integrate
import ipdb

def ghfunc(x,par,ipro=None,offset=None,plot_key=None):
    amp = fan([1.,par[1:10],par[15:18]],121)
    if par[0] != param.oldwid:
        param.wid = par[0]
        set_param
        param.oldwid = par[0]
    #S which axis is this summed over? I believe the first due to idl docs on 
    #S 'total'
    ipro = np.sum(amp*param.p,axis=0)
    
    return ipro/scipy.integrate.simps(ipro, x=param.x)

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
    
    #S if we want a transpose of the herm array
    if transpose:
        if n > 1:
            #S for 2 to n+1
            for ind in range(2,n+1):
                h[:,ind] = np.roll(h[ind-1,:],1)*2. -\
                    2.*float(ind-1)*h[:,ind-2]
    elif not transpose:
        if n > 1:
            for ind in range(2,n+1):
                h[ind,:] = np.roll(h[ind-1,:],1)*2. -\
                    2.*float(ind-1)*h[ind-2,:]

    #S return the array of the coeeficients
    return h


def set_param(param):
    
    #S number of terms
    n = 15
    #S number of oversampled points(?)
    nx = 121

    if not param.set:
        make_herm(n-1,c)
        param.coeff = c
        pow = np.arange(n)
        param.powarr = fan(pow,n)
        xarr = fan(param.x,n,transpose=True)
        param.zarr = xarr**pow

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
    z = make_herm(5)
