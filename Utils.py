import numpy as np
import json
import matplotlib.pyplot as plt
import os
from ImportData import ProcessExperiment
from sklearn.linear_model import LinearRegression
from PoincareSection import *

def RunExperimentDict(dic, NAME, dim=4):
    return RunExperiment(dic["w0"],dic["w1"],dic["q0"],dic["E"],dic["q1"],dic["p1"],dic["Tstep"],dic["Tf"],dic["Tprint"],dic["solver"],NAME,dim=dim)

def RunExperiment(w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, NAME, dim=4):
    p0square = 2.*E/w0 - (w1/w0)*(p1*p1+q1*q1) - q0*q0 -(2./w0)*q0*q0*q1 + (2./(3.*w0))*q1*q1*q1

    if p0square < 0:
        return None
        
    ID = NAME + "_live"
    
    # w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, id
    params = [w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, ID]
    params = " ".join([str(p) for p in params])
    
    dic = {}
    dic["w0"] = w0
    dic["w1"] = w1
    dic["q0"] = q0
    dic["E"] = E
    dic["q1"] = q1
    dic["p1"] = p1
    dic["Tstep"] = Tstep
    dic["Tf"] = Tf
    dic["Tprint"] = Tprint
    dic["solver"] = solver
    dic["ID"] = ID

    with open("Runs/" + ID + '.json', 'w') as f:
        json.dump(dic, f)
        
    command = "./GenOrbits2.o " + params
    os.system(command)
    
    exp = ProcessExperiment(ID, "Runs", dim)
    
    return exp


def HH_Hamiltonian(orbit, w0, w1):
    q0 = orbit[:, 0]
    q1 = orbit[:, 1]

    p0 = orbit[:, 2]
    p1 = orbit[:, 3]

    osc_harm_0 = (w0 / 2.0) * (np.power(q0, 2) + np.power(p0, 2))
    osc_harm_1 = (w1 / 2.0) * (np.power(q1, 2) + np.power(p1, 2))

    coupling_0 = np.multiply( np.power(q0, 2), q1 )
    coupling_1 = -(1.0/3.0) * np.power(q1, 3)

    return osc_harm_0 + osc_harm_1 + coupling_0 + coupling_1

def OnlyValid(arr):
    if arr.ndim == 1:
        nans = np.isnan(arr)
        infs = np.isinf(arr)
        return arr[np.invert(np.logical_or(nans, infs))]
    else:
        nans = np.any(np.isnan(arr), axis=1)
        infs = np.any(np.isinf(arr), axis=1)
        return arr[np.invert(np.logical_or(nans, infs)), :]

def ValidRatio(arr):
    return OnlyValid(arr).size / arr.size

def LinearFit(X, y, weights=None):
    X = X.reshape(-1,1)

    reg = LinearRegression().fit(X, y, weights)
    score = reg.score(X, y, weights)

    intercept = reg.intercept_
    slope = reg.coef_[0]

    return slope, intercept, score


def PresentExperiment(exp, time=0, poinc=None):
    orb = exp["orbit"]
    times = np.linspace(0, exp["Tf"], orb.shape[0])
    
    plane = Plane(np.array([1.0, 0, 0, 0]), 0)
    sec = PoincareSection(plane)
    
    stop = int(time*orb.shape[0]/exp["Tf"])

    if poinc is None:
        inters_pre = np.array(sec.ComputeFromOrbitAxis(orb[0:stop,:]))
        inters_post = np.array(sec.ComputeFromOrbitAxis(orb[stop:,:]))
    else:        
        inters_pre = np.array(sec.ComputeFromOrbitAxis(poinc[0:stop,:]))
        inters_post = np.array(sec.ComputeFromOrbitAxis(poinc[stop:,:]))
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,9))
    
    title = "E=" + str(exp["E"]) + " q1=" + str(exp["q1"]) + " p1=" + str(exp["p1"])
    
    fig.suptitle(title, fontsize=14)
    
    
    ax[0,0].set_title("Poincare section q0=0")
    ax[0,0].set_xlabel("q1")
    ax[0,0].set_ylabel("p1")
    ax[0,0].scatter(inters_pre[:,1], inters_pre[:,3], color='r', s=0.5)
    ax[0,0].scatter(inters_post[:,1], inters_post[:,3], color='b', s=0.5)
    
    
    ax[0,1].set_title("SALI")
    ax[0,1].set_xlabel("t")
    ax[0,1].set_ylabel("SALI")
    ax[0,1].loglog(times[0:stop], exp["SALI"][0:stop], color='r', linewidth=0.5)
    ax[0,1].loglog(times[stop+1:], exp["SALI"][stop+1:], color='b', linewidth=0.5)
    
    ax[1,0].set_title("MEGNO")
    ax[1,0].set_xlabel("t")
    ax[1,0].set_ylabel("MEGNO")
    ax[1,0].plot(times[0:stop], exp["MEGNO"][0:stop], color='r', linewidth=0.5)
    ax[1,0].plot(times[stop+1:-1], exp["MEGNO"][stop+1:], color='b', linewidth=0.5)
    
    ax[1,1].set_title("FLI")
    ax[1,1].set_xlabel("t")
    ax[1,1].set_ylabel("FLI")
    ax[1,1].loglog(times[0:stop], exp["FLI"][0:stop], color='r', linewidth=0.5)
    ax[1,1].loglog(times[stop+1:], exp["FLI"][stop+1:], color='b', linewidth=0.5)
    
    fig.tight_layout()
    plt.show()
