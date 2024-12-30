import pandas as pd
from .cea import CEA
import numpy as np
def PCEA_test():
    ts = pd.read_csv(
        r"../data/demo.csv",
        index_col=0,
        header=0)
    ts = ts.dropna()
    ts = ts.iloc[:, [0]]
    ts = pd.concat([ts, ts], axis=1)

    ts.columns = ['var1', 'var2']
    cea = CEA(ts, is_binary_array=False, max_gap=1, max_gap_length=1, delta=1, direction="forward",
              threshold=[-np.inf, -0.5], tau_i=1)
    fig, axes = cea.plot_signal(relationship='interval')
    fig.save('../data/test.png')