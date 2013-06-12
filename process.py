#!/usr/bin/env python
import pylab as pl
from wfdbtools import rdhdr, rdann
from sys import exit
import datetime as dt
import time
from scipy import interpolate as ipl
from pprint import pprint
from warnings import warn
import subprocess
import common
import plot_results
global config

"""
+ Summary plot for all signals
  + beta_std by desease
  + beta_cv by desease
  + Boxed for every desease
+ Boxed plot for signal, all signals
+ Count stats for all signals
+ Different interpolations (esp. spline)
+ Sliding window
Compare window size
+ 2-sigma cutting signals
+ Grid
+ 3D interpolation
+ Cluster analysis
+ Optimize memory usage
+ Config files
Spectrogram
Lomb FFT
Parameters as command Parameters
Do not use pylab (or fix case if using Gtk)
"""


def variability(r_times):
    """
    Get HRV from RR times array.
    Parameters:
        r_times : ndarray
    Returns:
        time : ndarray
            Time vector
        rrs : ndarray
            HRV vector
    """

    time = pl.delete(r_times,-1)
    hrv = pl.zeros(len(r_times)-1, dtype='float')

    for i in range(0, len(r_times)-1):
        hrv[i] = (r_times[i+1]-r_times[i])* 1000
    
    assert pl.shape(time) == pl.shape(hrv)
    return time, hrv


def interpolate(x, y, order=2):
    """
    Interpolates x/y using cubic interpolation
    In:
        x, y : ndarray, Input data
    Out:
        xnew, ynew : ndarray, Interpolated data
    """
    assert pl.shape(x) == pl.shape(y)
    xnew = pl.arange(x.min(), x.max(), 1.0/config['INTERP_FREQ'], dtype='float')

    if order == 0:
        f = ipl.interp1d(x, y, kind='linear')
        return xnew, f(xnew)
    else:
        tck = ipl.splrep(x, y, k=order, s=0)
        ynew = ipl.splev(xnew, tck, der=0)
        return xnew, ynew


def fft_filter(time, signal, result, window_func=pl.hanning, \
               samp_freq=8, freq_limit = (0.0033, 0.04), preview=False):
    """
    Applies an Band Pass filter to signal in the frequency domain and plots
        signals and their spectrals.
    In:
        time : ndarray, relative RR times vector
        signal : ndarray, HRV vector
        result : dict, current fragment info
        freq_limit : frequencies for band pass filter, analyzed freq range
    Out:
        10*log10 of frequency vector
        10*log10 of specturm
    """

    window = window_func(len(signal)) # window
    signal_wed = signal * window

    n = len(signal)/2

    spec = (pl.absolute(pl.rfft(signal_wed))[:-1]/(n))**2 # entering to frequency domain        
    freq = pl.fftfreq(len(signal), 1.0/samp_freq)[:n] # get freqs axis values
    spec_filt = pl.zeros(spec.shape, dtype='float')
    freq_filt = pl.zeros(freq.shape, dtype='float') # same for filtered

    for i in range(n):  # filter by frequency
        if pl.logical_and(freq_limit[0] <= abs(freq[i]), abs(freq[i]) <= freq_limit[1]):
            freq_filt[i] = freq[i]
            spec_filt[i] = spec[i]
        else:
            freq_filt[i] = 0 # fill invalid frequencies with 0 value
            spec_filt[i] = 0

    spec_filt = pl.ma.masked_equal(spec_filt,0) # cutt off invalid values
    spec_filt = pl.ma.compressed(spec_filt)
    freq_filt = pl.ma.masked_equal(freq_filt,0)
    freq_filt = pl.ma.compressed(freq_filt)

    spec_filt_log = pl.log10(spec_filt) # for output
    freq_filt_log = pl.log10(freq_filt)

    return freq_filt_log, spec_filt_log


def process_signal(record, annotator, diagnosis=None, slide_rate=.5,\
        preview=False):
    global config
    info = rdhdr(record)
    #pprint(info)
    if diagnosis:
        info['diagnosis'] = diagnosis
    print "Loading record %s: %s %s %s" % (record, info['gender'], info['age'],\
        info['diagnosis'])

    rrs = read_rrs(record, annotator)

    time_filt, hrv_filt = common.filter2d(rrs[:,0], rrs[:,1],\
        algos=config['HRV_FILTER_ALGO'])
    time_interp, hrv_interp = interpolate(time_filt, hrv_filt, config['SPLINE_ORDER'])


    if preview:
        plot_results.plot_hrv([rrs[:,0], rrs[:,1]], [time_filt, hrv_filt], [time_interp, hrv_interp], record,\
            preview=preview)

    del time_filt, hrv_filt

    #signal_to_csv(record, time, hrv, info)

    window = config['INTERP_FREQ'] * config['WINDOW']

    if window % 2 != 0:
        warn("Window (%s) is not even. Adjust INTERP_FREQ and/or WINDOW parameters to be x**2" % (window, ))
    if not (window != 0 and ((window & (window - 1)) == 0)):
        warn("Window (%s) is not power of 2. Adjust INTERP_FREQ and/or WINDOW parameters" % (window, ))

    frag_count = int((len(time_interp) - window) / window / slide_rate)
    results = []
    for frag in range(0, frag_count): # get WINDOW secs exactly
        result = {'record': record, 'frag': frag+1}

        frag_beg = frag * window * slide_rate
        frag_end = frag * window * slide_rate + window

        result['time'] = common.elapsed_to_abs_time(pl.mean((time_interp[frag_beg],time_interp[frag_end])), info['base_time'])
        result['time_beg'] = common.elapsed_to_abs_time(time_interp[frag_beg], info['base_time'])
        result['time_end'] = common.elapsed_to_abs_time(time_interp[frag_end], info['base_time'])
        #result['time_from'] = elapsed_to_abs_time(time_interp[frag_beg], info['base_time']).strftime("%H:%M:%S.%f")
        #result['time_to'] = elapsed_to_abs_time(time_interp[frag_end], info['base_time']).strftime("%H:%M:%S.%f")
        
        #result['sec_from'] = time_interp[frag_beg]
        #result['sec_to'] = time_interp[frag_end]

        time_frag = pl.array(time_interp[frag_beg: frag_end])
        time_frag -= time_frag.min()
        hrv_frag = pl.array(hrv_interp[frag_beg: frag_end])

        result['std'] = hrv_frag.std()
        result['mean'] = hrv_frag.mean()
        result['cov'] = result['std']/result['mean']*100

        hrv_frag -= result['mean']

        freq_filt, fft_filt = fft_filter(time_frag, \
                hrv_frag,\
                result,\
                window_func=getattr(pl, config['FFT_WINDOW_FUNC']),\
                samp_freq=config['INTERP_FREQ'],\
                freq_limit = config['FREQ_LIMIT'],\
                preview=preview)
        del time_frag, hrv_frag

        # Dirty hack just for case when we have constant signal in window, 
        # and signal - mean = 0, so filter in fft_filter returns empty PSD
        if pl.shape(freq_filt) != pl.shape(fft_filt):
            print "shape(freq) != shape(spec). Maybe, got constant signal in window. Passing by..."
            continue
        
        _a,_b = common.approximate(freq_filt,fft_filt)

        result['beta'] = -_a
        if frag == 50 and preview:
            line_y = [_a*x+_b for x in freq_filt]
            plot_results.plot_beta(freq_filt, fft_filt, line_y, result, preview=preview)
            del line_y

        results.append(result)
        print "%(frag)03d/%(frag_count)03d: %(time_beg)s - %(time_end)s\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % dict(result.items()+{'frag_count':frag_count}.items())
        del frag_beg, frag_end, freq_filt, fft_filt

    if results:
        stats(results)
        results_to_csv(results, record, info)
        plot_results.plot_time(rrs, info, results, config['WINDOW'], preview=preview)
        if preview:
            plot_results.plot_beta_std(results, preview=preview)
            plot_results.plot_beta_cv(results, preview=preview)
    else:
        warn("No results")

    return rrs, info, results

def read_rrs(record, annotator):
    rrs = []
    #proc = subprocess.Popen(['rdann','-r',record,'-a',annotator,'-p','N','-c','0'],bufsize=-1,stdout=subprocess.PIPE)
    proc = subprocess.Popen(['/usr/bin/ann2rr','-r',record,'-a',annotator,'-p','N','-P','N','-V','s6','-i','s6','-c'],bufsize=-1,stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        rr = line.split()
        if rr[0][0].isdigit() and rr[1][0].isdigit():
            rrs.append([float(rr[0]),float(rr[1])*1000])
    assert rrs
    return pl.array(rrs)


def results_to_csv(results,record,info):
    outfile = open("results_%s.csv" % (record,), "w")
    outfile.write("%s,%s,%s,%s\n" % (record, info['gender'], info['age'], info['diagnosis']))
    for r in results:
        r_ = dict(r)
        r_['time_beg'] = r_['time_beg'].strftime("%H:%M:%S.%f")
        r_['time_end'] = r_['time_end'].strftime("%H:%M:%S.%f")
        outfile.write("%(time_beg)s,%(time_end)s,%(beta)s,%(std)s,%(cov)s,%(mean)s\n" % r_)
    outfile.close()


def stats(results):
    means = {'beta':0.0,'std':0.0, 'cov':0.0, 'mean':0.0}
    mins = means.copy()
    maxs = means.copy()
    medians = means.copy()
    stdevs = means.copy()
    for param in means.keys():
        mins[param] = min([v[param] for v in results])
        maxs[param] = max([v[param] for v in results])
        means[param] = pl.mean([v[param] for v in results])
        medians[param] = pl.median([v[param] for v in results])
        stdevs[param] = pl.std([v[param] for v in results])
    print "Min:\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % mins
    print "Max:\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % maxs
    print "Mean:\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % means
    print "Median:\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % medians
    print "Stdev:\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % stdevs


if __name__ == '__main__':
    global config
    config = common.load_config()
    batch = False
    if batch:
        for diag in config['SIGNALS']:
            for record in diag['records'].split():
                process_signal(record, diag['annotator'], diag['diagnosis'],
                    slide_rate=config['SLIDE_RATE'], preview=config['PREVIEW'])
    else:
        #records = '16483 16773 16795 17453'
        records = 'chf201'
        for record in records.split():
            process_signal(record, 'ecg', "HF2", config['SLIDE_RATE'], preview=True)
