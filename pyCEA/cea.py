"""
Created on 08 Jun 2023
Update on 08 Jun 2023
@author: Hanyu Jin
version: 0.1
"""
from __future__ import annotations

import warnings

import numpy as np
import proplot as pplt  # I strong suggest use high-level plot package to visualize
from pandas import DataFrame

from cea_core import *


class CEA(object):

    def __init__(self, si: DataFrame | ndarray, threshold: ndarray | list | tuple = None,
                 tau: int | list | tuple | ndarray = None,
                 delta: int | tuple | ndarray = None, operator: list | tuple = ('ge', 'le'),
                 max_gap_length: int | list | tuple | ndarray = None,
                 max_gap: int | list | tuple | ndarray = None,
                 is_binary_array=False, get_compound_event=True):
        """

        :param si: (m x n) numpy array, m denotes time and n denotes variable (such as SPEI, SRI, SMI)
        :param tau: (1 x n) numpy array, window size of the compound events (such as [1, 1, 1])
        :param threshold: (2 x n) numpy array, the threshold of events (such as [[-0.5, -0.5, -0.5], [-inf, -inf, -inf]])
        :param delta: (1 x n) numpy array, the minimum duration of signal
        :param operator:
        :param is_binary_array: if si is not bool array, and threshold must be provided
        """
        #

        self.n_sample = si.shape[0]
        self.n_var = si.shape[1]

        self.si = si
        self.threshold = threshold
        self.tau = tau if tau is not None else 1
        self.delta = delta if delta is not None else 1
        self.operator = operator  # the threshold, False means it is below the threshold
        self.max_gap_length = max_gap_length if max_gap_length is not None else 1
        self.max_gap = max_gap if max_gap is not None else 1
        self.is_binary_array = is_binary_array  # whether the input signal si is a bool matrix, True means it exceeds
        self.get_compound = get_compound_event

        self.var_name = None
        self.time = None
        self._validate_input()
        self.event_trip_info = None
        self.compound_event_trip_info = None

    def _validate_input(self):
        # the following is safety check for parameter si
        if (not isinstance(self.si, ndarray)) and (not isinstance(self.si, DataFrame)):
            raise TypeError(f'signal {type(self.si)} is not valid type, signal should be numpy ndarray or pandas '
                            f'DataFrame')
        else:
            if (self.si.dtypes != 'float').all() and (self.si.dtypes != 'bool'):
                try:
                    self.si = self.si.astype('float')
                except ValueError:
                    raise ValueError("sigal contains non-numeric value, please remove it")
        if isinstance(self.si, DataFrame):
            self.var_name = self.si.columns
            self.time = self.si.index
            self.si = self.si.values
        elif isinstance(self.si, ndarray):
            self.si = self.si
            self.var_name = [f"var{x}" for x in range(self.n_var)]
            self.time = list(range(self.n_var))
        if np.isnan(self.si).any():
            raise ValueError("signal contains NaN values, please remove them.")

        # the following is safety check for parameter threshold

        if self.threshold is not None:
            if (not isinstance(self.threshold, ndarray)) and (not isinstance(self.si, list)) and (
                    not isinstance(self.si, tuple)) and (not self.is_binary_array):
                raise TypeError(
                    f', and event threshold {type(self.threshold)} is not valid type, threshold should be numpy.ndarray, list or tuple')
            else:
                if isinstance(self.threshold, ndarray):
                    if (self.threshold.ndim <= 2) and (self.threshold[0] != 2):
                        if self.threshold.ndim == 1:
                            if self.threshold.size != 2:
                                raise ValueError(
                                    f"only provide bound {self.threshold}, at the least two bounds need to be provided")
                            else:
                                self.threshold = np.tile(self.threshold, (self.n_var, 1))
                        else:
                            if (self.threshold.size != self.n_var * 2):
                                raise ValueError(
                                    f"only provide {self.threshold[0] / 2} threshold {self.threshold}, but there is {self.n_var} "
                                    f"variables: {self.var_name}")


                    else:
                        raise ValueError(f"too manny thresholds {self.threshold}")

                    if self.threshold.dtype != 'float':
                        try:
                            self.threshold = self.threshold.astype('float')
                        except ValueError:
                            raise ValueError("threshold contains non-numeric value, please remove it")
                if isinstance(self.threshold, list) or isinstance(self.threshold, tuple):
                    if len(self.threshold) != 2:
                        raise ValueError(f"current provide thresholds {self.threshold} not eq 2")
                    else:
                        self.threshold = np.tile(np.array(self.threshold), (self.n_var, 1))
        else:
            if (not self.is_binary_array):
                self.threshold = np.tile([-np.inf, -0.5], (self.n_var, 1))
                warnings.warn(f"threshold is None, and current thresholds is {self.threshold}")

        valid_operators = ['g', 'ge', 'l', 'le']
        for i in self.operator:
            if i not in valid_operators:
                raise ValueError(f"'{i}' is not a valid relational operator. Valid operators are: g, ge, l, le.")

        # Safety check for parameter tau
        if self.tau is not None:
            self.tau = self._check_param(self.tau, "tau")

        # Safety check for parameter delta
        if self.delta is not None:
            self.delta = self._check_param(self.delta, "delta")

        # Safety check for parameter max_gap_length
        if self.max_gap_length is not None:
            self.max_gap_length = self._check_param(self.max_gap_length, "max_gap_length")

        # Safety check for parameter max_gap
        if self.max_gap is not None:
            self.max_gap = self._check_param(self.max_gap, "max_gap")

    def _check_param(self, param, name):
        if isinstance(param, (list, tuple)):
            if len(param) == 1:
                return np.tile(np.array(param, dtype=int), (self.n_var, 1))
            elif len(param) == self.n_var:
                return np.array(param, dtype=int).reshape(self.n_var, 1)
            else:
                raise ValueError(f"Length of {name} must be 1 or {self.n_var}, got {len(param)}.")
        elif isinstance(param, np.ndarray):
            if param.ndim > 2:
                raise ValueError(f"{name} dimensions must be <= 2.")
            if param.size == 1:
                return np.tile(param.flatten(), (self.n_var, 1))
            elif param.size == self.n_var:
                return param.reshape(self.n_var, 1)
            else:
                raise ValueError(f"Size of {name} must be 1 or {self.n_var}, got {param.size}.")
        elif isinstance(param, int):
            return np.full((self.n_var), param)
        else:
            raise TypeError(f"{name} must be int, list, tuple, or ndarray.")

    def run_cea(self):
        """
        run compound event analysis
        :return: None
        """
        self.__flag_event()
        self.__summary_once_events()
        if self.get_compound:
            self.__flag_compound_event()
        else:
            pass

    def plot_signal(self, axes=None, **kwargs):
        """
        plot each variable in a figure
        :param axes:
        :return:
        """
        if axes is None:
            fig, axes = pplt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, refaspect=3)

        n = 1
        for i in range(self.n_var):
            si_normalize = (self.si[:, i] - self.si[:, i].min()) / (self.si[:, i].max() - self.si[:, i].min())
            axes.plot(self.time, si_normalize + n)
            n += 1
        axes.format(xlocator=pplt.Locator('maxn', 5),
                    **kwargs)
    def summary(self):
        """

        :return:
        """

    def __flag_event(self):
        """
        this function flags the beginning, terminating location of events
        :return:
        """
        if not self.is_binary_array:
            event_binary = []
            for i in range(self.n_var):
                event_binary.append(apply_operator(self.si[:, i], self.threshold[i, :], self.operator))
            event_binary = np.array(event_binary).transpose()
        else:
            event_binary = self.si.copy()

        event_trip_info = []  # the location index of beginning and terminating of events
        for i in range(self.n_var):
            event_trip_info.append(
                find_consecutive(event_binary[:, i], self.delta[i], self.max_gap_length[i], self.max_gap[i]))
        self.event_trip_info = event_trip_info

    def __flag_compound_event(self):
        # 遍历所有相邻变量的组合 (例如 (ENSO, SPEI), (SPEI, WILDFIRE))

        self.compound_event_trip_info = find_compound_event(self.event_trip_info, self.tau)

    def __summary_once_events(self):
        results_ = []
        for i in range(self.n_var):
            idx = self.event_trip_info[i]
            results_container = pd.DataFrame(columns=["var", 'time', 'peak', 'duration', 'strength'])
            for jhat, j in enumerate(idx):
                if ((j[1] - j[0]) == 0) and self.delta[i] == 1:
                    event_sig = self.si[j[0], i]  # only one value
                else:
                    event_sig = self.si[slice(*j), i]

                peak = event_sig[np.argmax(np.abs(event_sig))]
                duration = event_sig.size
                strength = event_sig.mean()
                results_container.loc[jhat, 'var'] = self.var_name[i]
                results_container.loc[jhat, 'time'] = f"{self.time[j[0]]}-{self.time[j[1]]}"
                results_container.loc[jhat, 'peak'] = peak
                results_container.loc[jhat, 'duration'] = duration
                results_container.loc[jhat, 'strength'] = strength
            results_.append(results_container)
        results_ = pd.concat(results_)
        self.once_events_statistic = results_


if __name__ == '__main__':
    import pandas as pd

    ts = pd.read_csv(
        r"../data/demo.csv",
        index_col=0,
        header=0)
    ts = ts.dropna()
    ts = ts.iloc[:, [0, 1, 2]]
    CEA(ts, np.array([-np.inf, -0.5]), is_binary_array=False, max_gap=1, max_gap_length=1, delta=3).run_cea()
