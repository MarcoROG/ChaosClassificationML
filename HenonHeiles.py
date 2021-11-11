from Splitting import SplittingMethod, LeapFrog, SABA3, SABA3_corr
from Problem import Problem
from PoincareSection import Plane, PoincareSection
import matplotlib.pyplot as plt
import numpy as np

def T(h, state, params):
    w0 = params[0]
    w1 = params[1]

    state[0] += h * w0 * state[2]
    state[1] += h * w1 * state[3]

def V(h, state, params):
    w0 = params[0]
    w1 = params[1]

    state[2] -= h * state[0] * (2 * state[1] + w0)
    state[3] -= h * (state[0]*state[0] - state[1]*state[1] + w1*state[1])

def TV_V(h, state, params):
    w0 = params[0]
    w1 = params[1]

    F0 = -4.0*w1*state[0]*state[0]*state[0] + (4*w1-8*w0)*state[0]*state[1]*state[1] + (-8*w0*w0 -4*w1*w1)*state[0]*state[1] -2*w0*w0*w0*state[0]
    F1 = (-8*w0+4*w1)*state[0]*state[0]*state[1] + (-4*w0*w0 -2*w1*w1)*state[0]*state[0] -4*w1*state[1]*state[1]*state[1] + 6*w1*w1*state[1]*state[1] - 2*w1*w1*w1*state[1]

    state[2] += h * F0
    state[3] += h * F1
