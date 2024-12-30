"""
Created on 08 Jun 2023
Update on 08 Jun 2023
@author: Hanyu Jin
version: 0.1
"""
from __future__ import annotations

import warnings

import numpy as np
import proplot as pplt
from matplotlib.patches import Rectangle
from pandas import DataFrame
from numpy import ndarray
from cea_core import *


class CEA(object):

    def __init__(self, si: DataFrame | ndarray, threshold: ndarray | list | tuple = None,
                 delta: int | tuple | ndarray = None, operator: list | tuple = ('ge', 'le'),
                 max_gap_length: int | list | tuple | ndarray = None,
                 max_gap: int | list | tuple | ndarray = None,
                 is_binary_array=False, get_compound_event=True,
                 tau_i: int | ndarray = None,
                 tau_o: int | ndarray = None,
                 tau_p: int | ndarray = None,
                 direction=None
                 ):
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
        # parameters of signal variable event identification
        self.delta = delta if delta is not None else 1
        self.operator = operator  # the threshold, False means it is below the threshold
        self.max_gap_length = max_gap_length if max_gap_length is not None else 1
        self.max_gap = max_gap if max_gap is not None else 1

        self.var_name = None
        self.time = None

        # compound parameters
        self.is_binary_array = is_binary_array  # whether the input signal si is a bool matrix, True means it exceeds
        self.get_compound = get_compound_event
        self.tau_i = tau_i if tau_i is not None else 6
        self.tau_o = tau_o if tau_o is not None else 6
        self.tau_p = tau_p if tau_p is not None else 1
        self.direction = direction
        # initial variable to store cea results

        # the following is variable to store CEA results
        self.event_trip_info = None
        self.compound_event_trip_info = None
        self.coincident_rate = None

        self._validate_input()

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
            if (not isinstance(self.threshold, ndarray)) and (not isinstance(self.threshold, list)) and (
                    not isinstance(self.threshold, tuple)) and (not self.is_binary_array):
                raise TypeError(
                    f'event threshold {type(self.threshold)} is not valid type, threshold should be numpy.ndarray, list or tuple')
            else:
                if isinstance(self.threshold, ndarray):
                    if (self.threshold.ndim <= 2) and (self.threshold.shape[1] != 2):
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

                    elif self.threshold.ndim == 2:
                        pass
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

        # Safety check for compound event parameter tau
        if self.tau_i is not None:
            self.tau_i = self.__check_compound_param(self.tau_i, "tau_i")
        if self.tau_o is not None:
            self.tau_o = self.__check_compound_param(self.tau_o, "tau_o")
        if self.tau_p is not None:
            self.tau_p = self.__check_compound_param(self.tau_p, "tau_p")

        # Safety check for parameter delta
        if self.delta is not None:
            self.delta = self._check_param(self.delta, "delta")

        # Safety check for parameter max_gap_length
        if self.max_gap_length is not None:
            self.max_gap_length = self._check_param(self.max_gap_length, "max_gap_length")

        # Safety check for parameter max_gap
        if self.max_gap is not None:
            self.max_gap = self._check_param(self.max_gap, "max_gap")

    def __check_compound_param(self, param, name):
        if isinstance(param, int):
            return np.full([self.n_var, self.n_var], param)
        elif isinstance(param, np.ndarray):
            if param.ndim != 2:
                raise ValueError(f"{name} dimensions must be <= 2.")
            if param.ndim == 2:
                if param.shape != (self.n_var, self):
                    raise ValueError(f"{name} shape must ({self.n_var}, {self.n_var}).")
        else:
            raise TypeError(f"{name} must be int, list, tuple, or ndarray.")

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

    def run_cea(self,
                verbose=True,
                save_path=None):
        """
        run compound event analysis
        :return: None
        """

        self.__flag_event()
        self.__summary_one_variable_events()
        if self.get_compound:
            self.__flag_compound_event()
        else:
            pass

        self.__summary_compound_events()

        if save_path is not None:
            pointer = pd.ExcelWriter(save_path)
            self.one_variable_events_statistic.to_excel(pointer, sheet_name='one_variable')
            if self.direction in ["backward", "both"]:
                self.backward_compound_events_statistic.to_excel(pointer, sheet_name="backward")
            if self.direction in ["forward", "both"]:
                self.forward_compound_events_statistic.to_excel(pointer, sheet_name="forward")
            pointer.close()
        self.__calc_coincidence_rate()
        # significant test base surrogate datasets (routine of R coincidence package)
    def plot_signal(self, relationship='all', fig=None, axes=None, **kwargs):
        """
        plot each variable in a figure
        :param axes:
        :return:
        """
        if self.event_trip_info is None and self.compound_event_trip_info is None:  # make sure already have cea results
            self.run_cea()

        if axes is None:
            fig, axes = pplt.subplots(nrows=self.n_var, ncols=1, sharex=False, sharey=False, refaspect=3)

        variable_event_info = self.one_variable_events_statistic
        bounds = [-abs(variable_event_info['strength'].min()), abs(variable_event_info['strength'].min())]
        variable_event_info = variable_event_info.set_index('var')
        for i in range(self.n_var):
            ax = axes[i]
            one_variable_info = variable_event_info.loc[self.var_name[i], :]
            one_variable_info.index = range(one_variable_info.shape[0])

            idx_interval_total = []
            idx_intersect_total = []
            idx_containing_total = []

            if self.direction in ['both', 'backward']:
                idx_intersect = self.compound_event_trip_info['backward']['intersect']
                idx_containing = self.compound_event_trip_info['backward']['containing']
                idx_interval = self.compound_event_trip_info['backward']['idx_interval']
                idx_intersect = [x[i] for x in idx_intersect]
                idx_containing = [x[i] for x in idx_containing]
                idx_interval = [x[i] for x in idx_interval]

                idx_interval_total.extend(idx_interval)
                idx_intersect_total.extend(idx_intersect)
                idx_containing_total.extend(idx_containing)
            if self.direction in ['both', 'forward']:
                idx_intersect = self.compound_event_trip_info['forward']['intersect']
                idx_containing = self.compound_event_trip_info['forward']['containing']
                idx_interval = self.compound_event_trip_info['forward']['interval']
                idx_intersect = [x[i] for x in idx_intersect]
                idx_containing = [x[i] for x in idx_containing]
                idx_interval = [x[i] for x in idx_interval]
                idx_interval_total.extend(idx_interval)
                idx_intersect_total.extend(idx_intersect)
                idx_containing_total.extend(idx_containing)

            for j in one_variable_info.index:
                x0, x1 = one_variable_info.loc[j, 'index_start'], one_variable_info.loc[j, 'index_end']
                height = one_variable_info.loc[j, 'strength']
                if height >= 0:
                    fc = 'tab:orange'
                else:
                    fc = 'tab:blue'
                hatch = None
                if (x0, x1) in idx_containing_total:
                    hatch = '////'
                if (x0, x1) in idx_interval_total:
                    hatch =  'xxxx'
                if (x0, x1) in idx_intersect_total:
                    hatch = '++++'
                if hatch is not  None:
                    rect = Rectangle((x0, 0), x1 - x0, height, facecolor=fc, edgecolor='k', hatch=hatch)
                else:
                    rect = Rectangle((x0, 0), x1 - x0, height, facecolor=fc, edgecolor='k')

                ax.add_patch(rect)
            ax.format(ylabel=self.var_name[i],
                      xlim=(0, self.n_sample),
                      ylim=bounds)
            ax.hlines(0, 0, self.n_sample, lw=0.5, color='k')
        return fig, axes


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


        # The variable 'cea_result1' contains specific relationships.
        # The variable 'cea_results2' contains the overall index of compound events.
        cea_results = find_compound_event(self.event_trip_info, self.tau_i, self.tau_o, self.tau_p,
                                                         direction=self.direction)
        self.compound_event_trip_info = cea_results

    def __summary_one_variable_events(self):
        """

        :return:
        """
        results_ = []
        for i in range(self.n_var):
            idx = self.event_trip_info[i]
            results_container = pd.DataFrame(columns=['var', 'event_start', 'event_end', 'peak', 'duration', 'strength', "index_start", "index_end"])
            for jhat, j in enumerate(idx):
                if ((j[1] - j[0]) == 0) and self.delta[i] == 1:
                    event_sig = self.si[j[0], i]  # only one value
                else:
                    event_sig = self.si[slice(*j), i]
                if self.delta[i] == 1 and isinstance(event_sig, float):
                    event_sig = np.array([event_sig])
                peak = event_sig[np.argmax(np.abs(event_sig))]
                duration = event_sig.size
                strength = event_sig.mean()
                results_container.loc[jhat, 'var'] = self.var_name[i]
                results_container.loc[jhat, "evnet_start"] = self.time[j[0]]
                results_container.loc[jhat, "event_end"] = self.time[j[1]]

                results_container.loc[jhat, 'index_start'] = j[0]
                results_container.loc[jhat, "index_end"] =j[1]
                results_container.loc[jhat, 'peak'] = peak
                results_container.loc[jhat, 'duration'] = duration
                results_container.loc[jhat, 'strength'] = strength
            results_.append(results_container)
        results_ = pd.concat(results_)

        self.one_variable_events_statistic = results_

    def __summary_compound_events(self):
        """
        :save_path: the path save cea results
        This function convert compound results to DataFrame
        :return: None
        """
        def __core_dict2frame(arr, var_name):
            columns = [f"{x}_{i}" for x in var_name for i in ["start", "end"]]
            columns.extend(["relationship"])
            container = pd.DataFrame(columns=columns)
            n = 0
            for i in arr:
                arr_ = arr[i]
                for j in arr_:
                    arr__ = np.array(j)
                    arr__ = np.ravel(arr__)
                    time__ = self.time[arr__]
                    container.loc[n, 'relationship'] = i
                    container.loc[n, columns[:-1]] = time__
                    n += 1
            return container
        df_backward = __core_dict2frame(self.compound_event_trip_info['backward'], self.var_name)
        df_forward = __core_dict2frame(self.compound_event_trip_info['forward'], self.var_name)
        self.forward_compound_events_statistic = df_forward
        self.backward_compound_events_statistic = df_backward

    def __calc_coincidence_rate(self):

        coincidence_rate_backward = pd.DataFrame(columns=["interval", "intersect", "containing", "total"])
        coincidence_rate_forward = pd.DataFrame(columns=["interval", "intersect", "containing", "total"])
        for i in range(self.n_var):
            sample_var_length = len(self.event_trip_info[i])
            if self.direction in ["both", "backward"]:
                sample_intersect_length = len(self.compound_event_trip_info['backward']['intersect'])
                sample_containing_length = len(self.compound_event_trip_info["backward"]["containing"])
                sample_interval_length = len(self.compound_event_trip_info["backward"]['interval'])
                total_sample_length = sample_containing_length + sample_intersect_length + sample_interval_length
                coincidence_rate_backward.loc[self.var_name[i], :] = [x / sample_var_length for x in [sample_interval_length, sample_intersect_length, sample_containing_length, total_sample_length]]
            if self.direction in ["both", 'forward']:
                sample_intersect_length = len(self.compound_event_trip_info['forward']['intersect'])
                sample_containing_length = len(self.compound_event_trip_info["forward"]["containing"])
                sample_interval_length = len(self.compound_event_trip_info["forward"]['interval'])
                total_sample_length = sample_containing_length + sample_intersect_length + sample_interval_length
                coincidence_rate_forward.loc[self.var_name[i], :] = [x / sample_var_length for x in
                                                                    [sample_interval_length, sample_intersect_length,
                                                                     sample_containing_length, total_sample_length]]
        self.forward_coincidence_rate = coincidence_rate_forward
        self.backward_coincidence_rate = coincidence_rate_backward

if __name__ == '__main__':
    import pandas as pd

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