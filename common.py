#!/usr/bin/env python
import csv
import pylab as pl
from pprint import pprint
from warnings import warn

def load_config(file_name="wfdb.conf"):
    config = {}
    execfile(file_name, config)
    return config


def read_signal(file_name=None):
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
                    result['gender'] = row[1]
                    result['age'] = row[2]
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


def filter2d(x, y, axes=['y'], algos=['2sigma'], valid_delta_ratio=.2):
    xnew = pl.array(x, dtype='float32')
    ynew = pl.array(y, dtype='float32')
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

def filter1d(x, return_mask=True, algos=['2sigma']):
    xnew = pl.array(x, dtype='float32')
    mask = pl.ones(len(x), dtype='bool')
    for algo in algos:
        if algo == 'diff02':
            for i in range(0, len(xnew)-1):
                if (abs(xnew[i+1] / xnew[i] - 1) > valid_delta_ratio): # current rr differs more then 20% of previous one
                    mask[i] = False
            if not return_mask:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '2sigma':
            mean = pl.mean(xnew)
            std = pl.std(xnew)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < mean - 2*std, mean + 2*std < xnew[i]):
                    mask[i] = False
            if not return_mask:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '5per95':
            per5 = pl.percentile(xnew,1)
            per95 = pl.percentile(xnew,99)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < per5, per95 < xnew[i]):
                    mask[i] = False
            if return_mask:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        elif algo == '3per97':
            per3 = pl.percentile(xnew,3)
            per97 = pl.percentile(xnew,97)
            for i in range(0, len(xnew)):
                if pl.logical_or(xnew[i] < per3, per97 < xnew[i]):
                    mask[i] = False
            if return_mask:
                xnew = xnew * mask
                xnew = pl.ma.masked_equal(xnew,0)
                xnew = pl.ma.compressed(xnew)

        else:
           warn("Donno anything about such filtration algorithm")

    print "NB: Deleted %0.2f%% of array" % (float(len(x)-len(xnew))/len(x)*100,)
    if return_mask:
        return mask
    else:
        return xnew, mask


def signal_to_csv(record,time,hrv,info):
    outfile = open("%s.csv" % (record,), "w")
    for i in range(len(hrv)):
        t_full = info['base_time'] + dt.timedelta(seconds=time[i])
        line = "%s %s %s\n" % (i+1,t_full.strftime("%H:%M:%S.%f"),hrv[i])
        outfile.write(line)
    outfile.close()

