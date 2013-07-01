#!/usr/bin/env python
import csv
import pylab as pl
from pprint import pprint
from warnings import warn
import datetime as dt


def load_config(file_name="wfdb.conf"):
    """
    Loads config file.
    In:
        file_name : str, config file name
    Out:
        config : dict, config file contents
    """
    config = {}
    execfile(file_name, config)
    return config


def read_signal(file_name=None):
    """
    Read processing results signal.
    In:
        file_name : str, processing results signal file name, *.csv
    Out:
        result : dict, file contents
    """
    beta, std, cov, mean = [], [], [], []
    result = {'record':False}
    try:
        with open(file_name, 'rb') as f:
            reader = csv.reader(f)
            for row in reader:
                if result['record']: # Already read header
                    beta.append(float(row[2]))
                    std.append(float(row[3]))
                    cov.append(float(row[4]))
                    mean.append(float(row[5]))
                else:
                    result['record'] = row[0]
                    if row[1] in ['F','M']:
                        result['gender'] = row[1]
                    else:
                        result['gender'] = None
                    if row[2].isdigit():
                        result['age'] = row[2]
                    else:
                        result['age'] = None
                    result['diagnosis'] = row[3]
    except IOError:
        print 'Could not read %s' % (file_name,)
        return None
    if len(beta) > 0:
        result['beta'] = beta
        result['std'] = std
        result['cov'] = cov
        result['mean'] = mean
        return result
    else:
        print "Result is empty"
        return None


def filter2d(x, y, axes=['y'], algos=['2sigma']):
    """
    Perform 2D data filtration by selected exes.
    In:
        x : ndarray, X vector
        y : ndarray, Y vector
        axes : list, axes names which are used to choose filtered values. x, y or any combination
    Out:
        xnew : ndarray, filtered X
        ynew : ndarray, filtered Y
    """
    xnew = pl.array(x, dtype='float')
    ynew = pl.array(y, dtype='float')
    mask_x = pl.ones(len(x), dtype='bool')
    mask_y = pl.ones(len(y), dtype='bool')
    if 'y' in axes:
        mask_y = filter1d(y,algos=algos)        
    if 'x' in axes:
        mask_x = filter1d(x,algos=algos)
    mask = mask_x * mask_y
    xnew *= mask
    ynew *= mask
    
    xnew = pl.ma.masked_equal(xnew,0)
    xnew = pl.ma.compressed(xnew)
    ynew = pl.ma.masked_equal(ynew,0)
    ynew = pl.ma.compressed(ynew)

    assert pl.shape(xnew) == pl.shape(ynew)
    return xnew, ynew

def filter1d(x, mask_only=True, algos=['2sigma']):
    """
    Filter vector with selected algorithms.
    In:
        x : ndarray, input vector
        mask_only : bool, do not touch input vector, just find "bad" values
        algos : list of str, algos list to apply to input vector. The sequence is the same as in the list
    Out:        
        xnew : ndarray, filtered input vector (returned only if mask_only=False)
        mask : ndarray of bool, vector of the same length as x, where "0" represents masked values
    """
    xnew = pl.array(x, dtype='float')
    mask = pl.ones(len(x), dtype='bool')
    for algo in algos:
        if algo == 'diff02':
            for i in range(0, len(xnew)-1):
                if (abs(xnew[i+1] / xnew[i] - 1) > .2): # current rr differs more then 20% of previous one
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '2sigma':
            mean = pl.mean(xnew)
            std = pl.std(xnew)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < mean - 2*std, mean + 2*std < xnew[i]):
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '5per95':
            per5 = pl.percentile(xnew,5)
            per95 = pl.percentile(xnew,95)
            #print per5,per95
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < per5, per95 < xnew[i]):
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '3per97':
            per3 = pl.percentile(xnew,3)
            per97 = pl.percentile(xnew,97)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < per3, per97 < xnew[i]):
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '1per99':
            per1 = pl.percentile(xnew,1)
            per99 = pl.percentile(xnew,99)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] <= per1, per99 <= xnew[i]):
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == 'ho_moody':
            for i in range(2, len(xnew)-2):
                mean = pl.mean([xnew[i-2],xnew[i-1],xnew[i+1],xnew[i+2]])
                if pl.logical_or(xnew[i] < .8 * mean, 1.2 * mean < xnew[i]):
                    mask[i] = False
            if not mask_only:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        else:
           warn("Donno anything about such filtration algorithm")

    print "NB: Deleted %0.3f%% of array" % (float(len(x)-pl.sum(mask))/len(x)*100,)
    if mask_only:
        return mask
    else:
        return xnew, mask


def signal_to_csv(record,time,hrv,info):
    """
    Writes HRV signal to CSV.
    In:
        record : str, record name in WFDB format
        time : ndarray, R occurance time vector
        hrv : ndarray, HRV signal
        info : dict, record info
    """
    outfile = open("%s.csv" % (record,), "w")
    for i in range(len(hrv)):
        t_full = info['base_time'] + dt.timedelta(seconds=time[i])
        line = "%s %s %s\n" % (i+1,t_full.strftime("%H:%M:%S.%f"),hrv[i])
        outfile.write(line)
    outfile.close()


def elapsed_to_abs_time(seconds, base_time=None):
    """
    Convert relative (elapsed) time in seconds to astronomical datetime.
    In:
        seconds : float, elapsed seconds
        base_time : datetime, recording start time (if available)
    Out:        
        datetime
    """
    if not base_time:
        base_time = dt.datetime(1900, 01, 01, 0, 0, 0,tzinfo=None)
    return base_time + dt.timedelta(seconds=pl.floor(seconds),\
        microseconds=(seconds - pl.floor(seconds))*(10**6))


def approximate(x,y):
    """
    Linear approximation of y=f(x) using least square estimator.
    In:
        x : ndarray
        y : ndarray
    Out:
        a, b : float, as in a*x+b=y
    """
    assert pl.shape(x) == pl.shape(y)
    A = pl.vstack([x, pl.ones(len(x))]).T
    a, b = pl.lstsq(A, y)[0]
    return a, b
