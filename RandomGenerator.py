import os
import string
import random
import numpy as np
import time
import json
from multiprocessing import Pool

random.seed(1234)
NGEN = 10  
Tstep = 1e-3
Tprint = 1e2
Tf = 100000.0
solver = 1


def RunJob(command):
    os.system(command)

FIXID = ''.join(random.choices(string.ascii_uppercase, k=6))

generated = 0

runs = []

while generated < NGEN:
    if generated % 10 == 0:
        print(generated)

    w0 = 1.0
    w1 = random.choice([1.0, (np.sqrt(5.0)+1.0) / 2.0]) #random.uniform(1.0, 2.0)

    q0 = 0.0
    #Eratio = 0.85 #random.uniform(0.001, 0.5)    
    Ecrit = min(w0*w0*w0/24. + w0*w0*w1/8., w1*w1*w1/6.)
    E = 0.17#Eratio * Ecrit

    xmin = -np.sqrt((2.*E)/w1)
    xmax = np.sqrt(E/((w1/2.)-(1./3.)*(-xmin)))
    ymin = xmin
    ymax = -ymin

    q1 = random.uniform(xmin, xmax) 
    p1 = random.uniform(0, xmax) 

    p0square = 2.*E/w0 - (w1/w0)*(p1*p1+q1*q1) - q0*q0 -(2./w0)*q0*q0*q1 + (2./(3.*w0))*q1*q1*q1

    if p0square < 0:
        continue

    ID = FIXID + str(int(time.time())) + "_" + str(generated)   


    # w0, w1, q0, E, q1, p1, N, Tf, Nprint, solver, id
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

pool = Pool(processes=2)
pool.map(RunJob, runs)

