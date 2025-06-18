import numpy as np
import os
import re

'''
from: http://www.mirzatrokic.ca/FILES/codes/fracdiff.py
small modification: wrapped 2**np.ceil(...) around int()
https://github.com/SimonOuellette35/FractionalDiff/blob/master/question2.py
'''

def fast_fracdiff(x, d):
    import pylab as pl
    T = len(x)
    np2 = int(2 ** np.ceil(np.log2(2 * T - 1)))
    k = np.arange(1, T)
    b = (1,) + tuple(np.cumprod((k - d - 1) / k))
    z = (0,) * (np2 - T)
    z1 = b + z
    z2 = tuple(x) + z
    dx = pl.ifft(pl.fft(z1) * pl.fft(z2))
    return np.real(dx[0:T])


def get_weights(d, size):
    # thres>0 drops insignificant weights
    w = [1.]
    for k in range(1, size):
        w_ = -w[-1] / k * (d - k + 1)
        w.append(w_)
    w = np.array(w[::-1]).reshape(-1, 1)
    return w


def fracDiff_original_impl(series, d, thres=.01):
    # 1) Compute weights for the longest series
    w = get_weights(d, series.shape[0])
    # 2) Determine initial calcs to be skipped based on weight-loss threshold
    w_ = np.cumsum(abs(w))
    w_ /= w_[-1]
    skip = w_[w_ > thres].shape[0]
    # 3) Apply weights to values
    # df = {}
    output = {}
    for name in series.columns:
        seriesF = series[[name]].fillna(method='ffill').dropna()
        for iloc in range(skip, seriesF.shape[0]):
            loc = seriesF.index[iloc]
            if not np.isfinite(series.loc[loc, name]): continue  # exclude NAs
            output[loc] = np.dot(w[-(iloc + 1):, :].T, seriesF.loc[:loc])[0, 0]
        # df[name] = df_.copy(deep=True)
    # df = pd.concat(df, axis=1)
    return output


def get_weight_ffd(d, thres, lim):
    w, k = [1.], 1
    ctr = 0
    while True:
        w_ = -w[-1] / k * (d - k + 1)
        if abs(w_) < thres:
            break
        w.append(w_)
        k += 1
        ctr += 1
        if ctr == lim - 1:
            break
    w = np.array(w[::-1]).reshape(-1, 1)
    return w


def fracDiff_FFD_original_impl(series, d, thres=1e-5):
    import pandas as pd
    # 1) Compute weights for the longest series
    w = get_weight_ffd(d, thres, len(series))
    width = len(w) - 1
    # df = {}
    output = []
    for name in series.columns:
        seriesF, df_ = series[[name]].fillna(method='ffill').dropna(), pd.Series()
        output.extend([0] * width)
        for iloc1 in range(width, seriesF.shape[0]):
            loc0, loc1 = seriesF.index[iloc1 - width], seriesF.index[iloc1]
            if not np.isfinite(series.loc[loc1, name]):
                continue  # exclude NAs
            # df_[loc1] =
            output.append(np.dot(w.T, seriesF.loc[loc0:loc1])[0, 0])
        # df[name] = df_.copy(deep=True)
    # df = pd.concat(df, axis=1)
    return output


def frac_diff_ffd(x, d, thres=1e-5):
    w = get_weight_ffd(d, thres, len(x))
    width = len(w) - 1
    output = []
    output.extend([0] * width)
    for i in range(width, len(x)):
        output.append(np.dot(w.T, x[i - width:i + 1])[0])
    return np.array(output)

def plot_multi(data, cols=None, spacing=.1, **kwargs):
    from pandas import plotting
    
    # Get default color style from pandas - can be changed to any other color list
    if cols is None: cols = data.columns
    if len(cols) == 0: return
#    colors = getattr(getattr(plotting, 'style'), '_get_standard_colors')(num_colors=len(cols))

    # First axis
#    ax = data.loc[:, cols[0]].plot(label=cols[0], color=colors[0], **kwargs)
    ax = data.loc[:, cols[0]].plot(label=cols[0], **kwargs)
    ax.set_ylabel(ylabel=cols[0])
    lines, labels = ax.get_legend_handles_labels()

    for n in range(1, len(cols)):
        # Multiple y-axes
        ax_new = ax.twinx()
        ax_new.spines['right'].set_position(('axes', 1 + spacing * (n - 1)))
 #       data.loc[:, cols[n]].plot(ax=ax_new, label=cols[n], color=colors[n % len(colors)])
        data.loc[:, cols[n]].plot(ax=ax_new, label=cols[n])
        ax_new.set_ylabel(ylabel=cols[n])

        # Proper legend position
        line, label = ax_new.get_legend_handles_labels()
        lines += line
        labels += label

    ax.legend(lines, labels, loc=0)
    return ax

#------------------------------------------Main functions--------------------------------------------------------------------------------------------------------------


def loop_fracs_ADfuller_test(np_array, p_value, Max, step) :
    from statsmodels.tsa.stattools import adfuller
    fracCoef = 0
    fracs = fast_fracdiff(np_array, fracCoef)
    while (adfuller(fracs)[1] > p_value) :
        fracCoef += step
        fracs = fast_fracdiff(np_array, fracCoef)
        if fracCoef > Max :
            print('ERROR -- Data not stationary enough to perform Granger Causality Test')
            break
    return fracs, fracCoef, adfuller(fracs)[1]


def Granger_Test(np_array1, np_array2, max_lag, p_value): #np_array2 Granger causes np_array 1
    from statsmodels.tsa.stattools import grangercausalitytests
    N = np_array1.shape[0]
    if (np_array2.shape[0] != N):
        print('ERROR -- Wrong dimension in Granger Test -- Provided arrays must have same dim')    
    GrangerData = np.zeros((N,2))
    for i in range(N) :
        GrangerData[i][0] = np_array1[i]
        GrangerData[i][1] = np_array2[i]
    results = grangercausalitytests(GrangerData, maxlag = max_lag, verbose = False)
    
    # loop on lags to check if the null hypthesis is rejected for some of them
    for lag in range(1,max_lag) :
        if results[lag][0]['ssr_ftest'][1] < p_value :
            return [lag, results[lag][0]['ssr_ftest'][1]]
    return('test failed')
    

#---------------------------------------------------Main----------------------------------------------------------------------------------------------------------

Data1 = "Data from time series 1"
Data2 = "Data from time series 2"

Stationary_Data1, frac_coef1, pvalue1 = loop_fracs_ADfuller_test(Data1, 0.05, 2, 0.05)
Stationary_Data2, frac_coef2, pvalue2 = loop_fracs_ADfuller_test(Data2, 0.05, 2, 0.05)

# time series 1 causes time series 2 ?
result = Granger_Test(Stationary_Data1, Stationary_Data2, 10, 0.05)
