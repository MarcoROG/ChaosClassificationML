import numpy as np

def SALI(tangA, tangB):
    tangA = tangA / np.linalg.norm(tangA, axis=1)[:, None]
    tangB = tangB / np.linalg.norm(tangB, axis=1)[:, None]

    diff = tangA - tangB
    sum = tangA + tangB

    # Row by row norm of diff
    diffNorm = np.linalg.norm(diff, axis=1)
    # Row by row norm of sum
    sumNorm = np.linalg.norm(sum, axis=1)

    return np.minimum(diffNorm, sumNorm)

def MEGNO_nonorm(tang, t_max):
    N = tang.shape[0] - 1
    h = t_max / N
    times = np.arange(N+1) * h
    times = times[1:]

    incr = np.diff(tang, axis=0)
    tang = 0.5 * (tang[0:-1,:] + tang[1:,:])

    tangNorm = np.power(np.linalg.norm(tang, axis=1), 2)
    incrNorm = np.einsum('ij,ij->i', incr, tang)
    fraction = np.divide(incrNorm, tangNorm)

    integrand = np.multiply(fraction, times) * h # integrand with the dt
   
    
    invtimes = 1.0 / times
    int_cumsum = np.cumsum(integrand)

    Y = 2 * np.multiply(invtimes, int_cumsum)
    Ybar = np.multiply(invtimes, np.cumsum(Y*h))

    return Ybar

def MEGNO(tang, t_max):
    N = tang.shape[0] - 1
    h = t_max / N
    times = np.arange(N+1) * h
    times = times[1:]
    
    norms = np.linalg.norm(tang, axis=1)
    norms_plus = norms[1:]
    norms = norms[0:-1]
    
    tang_plus = tang[1:,:]
    tang = tang[0:-1,:]
    
    scp = 2 * np.einsum('ij,ij->i', tang, tang_plus)
        
    numerator = np.multiply(norms, np.power(norms_plus, 2) - 1.0)
    denominator = np.multiply(norms, np.power(norms_plus, 2)) + norms + scp
    
    fraction = 2 * np.divide(numerator, denominator)
        
    integrand = np.multiply(fraction, times) # integrand without the dt because it cancels
   
    invtimes = 1.0 / times
    int_cumsum = np.cumsum(integrand)

    Y = 2 * np.multiply(invtimes, int_cumsum)
    Ybar = np.multiply(invtimes, np.cumsum(Y*h))

    return Ybar


def FLI(tangs):
    normsteps = np.stack([np.cumprod(np.linalg.norm(tang, axis=1)) for tang in tangs])
    
    return normsteps.max(axis=0)
