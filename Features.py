import numpy as np
from scipy.fft import fft
from Utils import *

def ValidFeature(exp):
    orbit = ValidRatio(exp["orbit"])
    tangents = min([ValidRatio(t) for t in exp["tangents"]])
    exp["feature_valid"] = min(orbit, tangents)
    return exp

def DeltaEFeature(exp):
    valid = OnlyValid(exp["orbit"])
    Es = HH_Hamiltonian(valid, exp["w0"], exp["w1"])    
    exp["feature_deltaE"] = np.log(np.max(Es) - np.min(Es))
    return exp

def MegnoLinearFitFeature(exp, weight=None):
    megno = exp["MEGNO"]
    times = np.linspace(0, exp["Tf"], megno.shape[0]) # Get all times
    megno = OnlyValid(megno)
    times = times[0:megno.shape[0]] # Only the valid ones

    tF = times[-1]
    
    if weight is None:
        weights = np.ones_like(times)
    else:
        weights = 1.0 / (1.0 + np.power(tF - times, weight))
    
    slope, intercept, score = LinearFit(times, megno, weights)

    exp["feature_megno_slope"] = slope
    exp["feature_megno_intercept"] = intercept
    exp["feature_megno_score"] = score

    return exp

def SALIExpFitFeature(exp, weight=None):
    l_sali = np.log(np.clip(exp["SALI"], 1e-15, a_max=None))
    times = np.linspace(0, exp["Tf"], l_sali.shape[0]) # Get all times

    l_sali = OnlyValid(l_sali)
    times = times[0:l_sali.shape[0]] # Only the valid ones
    tF = times[-1]

    if weight is None:
        weights = np.ones_like(times)
    else:
        weights = 1.0 / (1.0 + np.power(tF - times, weight))
    
    slope, intercept, score = LinearFit(times, l_sali, weights)

    exp["feature_sali_slope"] = slope
    exp["feature_sali_intercept"] = intercept
    exp["feature_sali_score"] = score

    return exp

def SALIFFTFeature(exp):
    salis = exp["SALI"]
    sali = OnlyValid(salis)
    ff = np.abs(fft(sali))
    
    simdiff = np.max(np.abs(ff - ff[::-1]))
    max_v = np.max(ff)
    
    peak = np.sum(ff[0:10]) / np.sum(ff)

    steps = len(ff)
    start = int(steps * 0.2)
    stop = int(steps * 0.8)

    tv = np.sum(np.abs(np.diff(ff[start:stop])))
    
    exp["feature_sali_max_freq"] = max_v
    exp["feature_sali_peak_ratio"] = peak
    exp["feature_sali_symmetric"] = simdiff / max_v
    exp["feature_sali_tv_freq"] = tv / max_v

    return exp

def FLIFFTFeature(exp):
    flis = exp["FLI"]
    fli = np.clip(OnlyValid(flis), 0, 1e100)
    ff = np.abs(fft(fli))
    
    simdiff = np.max(np.abs(ff - ff[::-1]))
    max_v = np.max(ff)
    
    peak = np.sum(ff[0:10]) / np.sum(ff)

    steps = len(ff)
    start = int(steps * 0.2)
    stop = int(steps * 0.8)

    tv = np.sum(np.abs(np.diff(ff[start:stop])))
    
    exp["feature_fli_max_freq"] = np.log(max_v)
    exp["feature_fli_peak_ratio"] = peak
    exp["feature_fli_symmetric"] = simdiff / max_v
    exp["feature_fli_tv_freq"] = tv / max_v

    return exp

def FLIExpFitFeature(exp, weight=None):
    l_fli = np.log(exp["FLI"])
    times = np.linspace(0, exp["Tf"], l_fli.shape[0]) # Get all times

    l_fli = OnlyValid(l_fli)
    times = times[0:l_fli.shape[0]] # Only the valid ones
    tF = times[-1]

    if weight is None:
        weights = np.ones_like(times)
    else:
        weights = 1.0 / (1.0 + np.power(tF - times, weight))
    
    slope, intercept, score = LinearFit(times, l_fli, weights)

    exp["feature_fli_slope"] = slope
    exp["feature_fli_intercept"] = intercept
    exp["feature_fli_score"] = score

    return exp

