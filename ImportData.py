import struct
import numpy as np
import json
import array
import os
from Indicators import *

def ImportData(fileName):
    with open(fileName, 'rb') as f:
        Nrows = struct.unpack('i', f.read(4))[0]
        Ncols = struct.unpack('i', f.read(4))[0]
    
    data = array.array('d')
    with open(fileName, 'rb') as fin:
        n = os.fstat(fin.fileno()).st_size // 8
        data.fromfile(fin, n)
        
    data = data[1:]
    Nrows = int(len(data) / Ncols)
    data = np.array(data).reshape((Nrows,Ncols))
        
    return data

def ImportExperiment(name, folder):
    with open(folder + "/" + name + '.json', 'r') as fp:
        experiment = json.load(fp)
    
    orbit = ImportData(folder + "/RUNHIST" + name + ".bin")
    tangent = ImportData(folder + "/RUNTANG" + name + ".bin")

    return experiment, orbit, tangent


def ProcessExperiment(name, folder, dim):
    exp, orbit, tang = ImportExperiment(name, folder)

    tangents = []
    for i in range(int(tang.shape[1]/dim)):
        tangents.append(tang[:,i*dim:(i+1)*dim])

    sali = SALI(tangents[0], tangents[1])
    fli = FLI(tangents)
    megno = MEGNO(tangents[0], exp["Tf"])

    exp["SALI"] = sali
    exp["FLI"] = fli
    exp["MEGNO"] = megno
    exp["orbit"] = orbit
    exp["tangents"] = tangents

    return exp
