import sys
sys.path.append("../")
from ImportData import ProcessExperiment
from multiprocessing import Pool
import glob, os
from PoincareSection import *
import pandas as pd
from Utils import HH_Hamiltonian
import Features
from Features import *
from LabelDataset import *
from Utils import RunExperiment
import tqdm

NPROC = 4
FOLDER = "Runs_def_res"
treshes = [0.1, 0.05, 0.01, 0.005, 0.001]

def GetLen(row):
    return len(row["SALI"])

def GetSlice(row, ratio):
    Nold = GetLen(row)
    Nnew = int((Nold - 1) * ratio + 1)
    
    short = dict(row)
    
    short["SALI"] = row["SALI"][0:Nnew]
    short["FLI"] = row["FLI"][0:Nnew]
    short["MEGNO"] = row["MEGNO"][0:Nnew]
    
    return short

def GetFeatures(row):
    exps = []
    N = GetLen(row)

    pick = min(512, N)
    
    idx = np.round(np.linspace(0, N-1, pick)).astype(int)

    exp = dict(row)
    exp["feature_FLI"] = np.log(np.clip(row["FLI"], 1e-15, 1e100))[idx]
    exp["feature_SALI"] = row["SALI"][idx]
    exp["feature_MEGNO"] = row["MEGNO"][idx]

    return exp

def Finalise(row):
    todrop = ["FLI", "SALI", "MEGNO", "orbit", "tangents", "inters", "allsame",
              "label_FLI", "label_MEGNO", "label_SALI"]
    
    smaller = dict(row)

    for col in todrop:
        del smaller[col]

    return smaller

def Process(data):
    sali_tresh = 1e-4
    megno_tresh = 10
    fli_tresh = 1.0

    file, FOLDER, treshes = data
    exp = ProcessExperiment(file, "../" + FOLDER, 4)
    
    exp["label_MEGNO"] = MEGNO_Label(exp, megno_tresh)
    exp["label_SALI"] = SALI_Label(exp, sali_tresh)
    fl = FLI_Label(exp,fli_tresh)
    exp["label_FLI"] = fl
       
    exp["allsame"] = (exp["label_MEGNO"] == exp["label_SALI"]) and (exp["label_MEGNO"] == exp["label_FLI"])
    exp["inters"] = None
        
    if not exp["allsame"]:
        return None

    exp["label"] = exp["label_MEGNO"]

    results = []

    for idx, tresh in enumerate(treshes):
        copy = dict(exp)
        copy = GetSlice(copy, tresh)
        N = GetLen(copy)

        copy = GetFeatures(copy)
        copy = Finalise(copy)

        results.append(copy)

    return results


# Begin
os.chdir("../")

# Load all experiments
count = 0
files = []

os.chdir(FOLDER)
for file in glob.glob("[!test_]*.json"):
    if count % 50 == 0:
        print(file)
        print(count)

    count+=1

    file = file.split(".")[0]
    files.append((file,FOLDER,treshes))

pool = Pool(processes=NPROC)
results = list(tqdm.tqdm(pool.imap(Process, files), total=len(files)))
results = [r for r in results if r is not None]

results = list(map(list, zip(*results)))

for idx,res in enumerate(results):
    print(idx)
    N = int((res[0]["Tf"] / res[0]["Tprint"]) * treshes[idx])
    df = pd.DataFrame(res)
    df.to_pickle("Series_512_res" + str(N))
