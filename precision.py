import numpy as np

def hexize(a, prc=np.float32):
    return a.values.astype(prc)[:]

def hexize_val(a, prc=np.float32):
    return a.values.astype(prc)