import numpy as np

class Problem:
    def __init__(self, FluxA, FluxB, params, eps, FluxCorr = None):
        self.FluxA = FluxA
        self.FluxB = FluxB
        self.FluxCorr = FluxCorr
        self.params = params
        self.epsilon = eps