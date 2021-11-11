from abc import ABC, abstractmethod
from joblib import dump, load
import json
import random
import string
from Utils import RunExperimentDict
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial
from Features import *
from tensorflow import keras
import os

# Computes an Hashvalue for a dictory by hashing all the keys
def hash_dict(dic):
    hashval = 0
    for key in dic.keys():
        if isinstance(dic[key], float):
            hashval = hashval ^ hash(dic[key])

    return hashval

# Generates a random string
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

# Divides a list in chunks of size n
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

class Model(ABC):
    def __init__(self, name):
        with open(name + "_info.json", "r") as f:
            data = json.load(f)

        self.data = data
        self.Tf = data["Tf"]
        self.Tprint = data["Tprint"]
        self.Tstep = data["Tstep"]
        self.solver = data["solver"]

    @abstractmethod
    def GetFeatures(self, exp):
        pass

    def PredictFeatures(self, feats):
        return self.model.predict_proba(feats)[:,1]

    def PredictIntegration(self, exp):
        feats = self.GetFeatures(exp)
        return self.PredictFeatures(feats)

    def GetIntegration(self, dic, dim=4, cleanup=False):
        dic["Tstep"] = self.Tstep
        dic["Tprint"] = self.Tprint
        dic["Tf"] = self.Tf
        dic["solver"] = self.solver

        NAME = id_generator() + str(hash_dict(dic))
        exp = RunExperimentDict(dic, NAME, dim)
    
        if cleanup:
            os.remove(NAME + "_live.json")
            os.remove("RUNHIST" + NAME + "_live.json")
            os.remove("RUNTANG" + NAME + "_live.json")

        return exp

    def PredictInitialData(self, dic, dim=4, cleanup=False):
        exp = self.GetIntegration(dic, dim, cleanup)
        return self.PredictIntegration(exp)

    def PredictGrid(self, dics, dim=4, batchsize=64, concurrency=None):
        # Determine how many CPU cores to use
        if concurrency is None:
            concurrency = multiprocessing.cpu_count()

        # Make sure the function "GetIntegration" is correctly parallelizable
        integrate = partial(self.GetIntegration, cleanup=False, dim=dim)

        # Save the prediction
        preds = []

        # Split the list of simulation "dics" in chunks of size "batchsize"
        for chunk in chunks(dics, batchsize):
            # Work on each chunk separately
            with Pool(concurrency) as p:
                # Integrate the chunk in parallel (efficient because each integration is independent)
                integrations = p.map(integrate, chunk) 
                
                # Compute features from all the integrations in parallel (efficient for the same reason)
                features = p.map(self.GetFeatures, integrations)
                # Convert to numpy array containing the features of all the simulations in the chunk
                features = np.squeeze(np.array(features))
                
                # Predict a class for all the elements in chunk (already parallel thanks to Keras)
                pred_chunk = self.PredictFeatures(features)

                # Append the chunk predictions to the final array
                if len(preds) == 0:
                    preds = pred_chunk
                else:
                    preds = np.append(preds, pred_chunk)

        # Format the result
        preds = np.array(preds).flatten()
        return preds            

class FeatureModel(Model):
    def __init__(self, name):
        super(FeatureModel, self).__init__(name)        
        self.model = load(name + "_model")


class SequenceModel(Model):
    def __init__(self, name):
        super(SequenceModel, self).__init__(name)
        self.model = keras.models.load_model(name + "_model.h5")


class RFModel(FeatureModel):
    def GetFeatures(self, exp):     
        exp = ValidFeature(exp)
        exp = DeltaEFeature(exp)
        exp = MegnoLinearFitFeature(exp, 0.5)
        exp = SALIExpFitFeature(exp, 0.5)
        exp = SALIFFTFeature(exp)
        exp = FLIExpFitFeature(exp, 0.5)
        exp = FLIFFTFeature(exp)

        vals = np.array([exp[v] for v in exp.keys() if v.startswith("feature_")])

        return vals.reshape(1, -1)

class SVMModel(FeatureModel):
    def GetFeatures(self, exp):     
        exp = ValidFeature(exp)
        exp = DeltaEFeature(exp)
        exp = MegnoLinearFitFeature(exp, 0.5)
        exp = SALIExpFitFeature(exp, 0.5)
        exp = SALIFFTFeature(exp)
        exp = FLIExpFitFeature(exp, 0.5)
        exp = FLIFFTFeature(exp)

        vals = np.array([exp[v] for v in exp.keys() if v.startswith("feature_")])

        return vals.reshape(1, -1)

class CNNModel(SequenceModel):

    def PredictFeatures(self, feats):
        if len(feats.shape) == 3:
            feats = feats[..., np.newaxis]
        return self.model.predict(feats)[:,0]

    def GetFeatures(self, row):
        N = int(self.Tf / self.Tprint)

        pick = min(512, N+1)
        
        idx = np.round(np.linspace(0, N-1, pick)).astype(int)

        vals = np.zeros((1, 3, pick, 1))
        vals[0,0,:,0] = np.log(np.clip(row["FLI"], 1e-15, 1e100))[idx]
        vals[0,1,:,0] = row["SALI"][idx]
        vals[0,2,:,0] = row["MEGNO"][idx]

        return vals

