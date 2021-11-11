import os
import string
import random
import numpy as np
import time
import json
from multiprocessing import Pool
from ImportData import ProcessExperiment
from PoincareSection import *
import pandas as pd
import tqdm

random.seed(1234)
NGEN = 512
CHUNKSIZE = 32
NPROC = 4
Tstep = 1e-3
Tprint = 1e-1
Tf = 10000.0
solver = 1

#w1 = np.sqrt(2)
#Es = [0.1565, 0.1575, 0.16] #for sqrt(2)

#w1 = 1.0
#Es = [ 0.11, 0.125, 0.13 ] #for 1

#w1 = (np.sqrt(5.0)+1.0) / 2.0
#Es = [0.17, 0.19, 0.21, 0.23]


w1 = (np.sqrt(5.0)+1.0) / 2.0
Es = [0.2, 0.22]

w0 = 1.0

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def RunJob(command):
    os.system(command)
    return command
    


FIXID = ''.join(random.choices(string.ascii_uppercase, k=6))

generated = 0

runs = []

q0 = 0.0

for E in Es:
    xmin = -np.sqrt((2.*E)/w1)
    xmax = np.sqrt(E/((w1/2.)-(1./3.)*(-xmin)))
    ymin = xmin
    ymax = -ymin

    q1s = np.linspace(xmin, xmax, NGEN)
    p1 = 0.0

    for q1 in q1s:
        if generated % 10 == 0:
            print(generated)

        p0square = 2.*E/w0 - (w1/w0)*(p1*p1+q1*q1) - q0*q0 -(2./w0)*q0*q0*q1 + (2./(3.*w0))*q1*q1*q1

        if p0square < 0:
            continue

        ID = FIXID + str(int(time.time())) + "_" + str(generated)   


        # w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, id
        params = [w0, w1, q0, E, q1, p1, Tstep, Tf, Tprint, solver, ID]
        params = " ".join([str(p) for p in params])
        print(params)
        #os.system("./GenOrbits2.o " + params)

        runs.append("./GenOrbits2.o " + params)

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

        with open('Runs/'+ ID + '.json', 'w') as f:
            json.dump(dic, f)

        generated += 1

print("Running")

pool = Pool(processes=NPROC)
for idx, chunk in enumerate(chunks(runs, CHUNKSIZE)):
    print("Chunk " + str(idx+1) + " / " + str(int(generated / CHUNKSIZE)))

    results = list(tqdm.tqdm(pool.imap(RunJob, chunk), total=len(chunk)))
