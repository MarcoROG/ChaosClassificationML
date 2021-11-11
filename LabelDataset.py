import numpy as np
from Utils import OnlyValid

def MEGNO_Label(exp, tresh):
    return np.abs(exp["MEGNO"][-1]) > tresh

def SALI_Label(exp, tresh):
    return np.abs(exp["SALI"][-1]) < tresh

def FLI_Label(exp, tresh):
    fli = OnlyValid(exp["FLI"])
    last = fli[-1] 
    res = last > tresh * exp["Tf"]
    return res
