#!/usr/bin/env python
import pylab as pl
from wfdbtools import rdhdr, rdann
from sys import exit
import datetime as dt
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
    
    #if len(time) > len(hrv):
    #    hrv[-1] = 0
    #print pl.shape(r_times), pl.shape(time),pl.shape(hrv)
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

    n = len(signal)

    window = window_func(n) # window
    signal_wed = signal * window

    if False:
        spec = pl.fftshift(pl.fft(signal_wed)) # entering to frequency domain
        spec_filt = pl.zeros(spec.shape, dtype=complex)
        # fftshift moves zero-frequency component 
        # to the center of the array
        freq = pl.fftshift(pl.fftfreq(n, 1.0/samp_freq)) # get freqs axis values
        freq_filt = pl.zeros(freq.shape, dtype=float) # same for filtered

        for i in range(n):  # filter by frequency
            if freq_limit[0] <= abs(freq[i]) and abs(freq[i]) <= freq_limit[1]:
                freq_filt[i] = freq[i]
                spec_filt[i] = spec[i]
            else:
                freq_filt[i] = 0 # fill invalid frequencies with small value
                spec_filt[i] = 0

        signal_filt = pl.absolute(pl.ifft(pl.ifftshift(spec_filt))) # get filtered signal
        signal_filt_uwed = signal_filt.copy()
        #signal_filt_uwed[1:-1] /= window[1:-1]
        #signal_filt_uwed = signal_filt_uwed[:len(time)]
        #signal = signal[:len(time)]
       
        spec_abs = pl.absolute(spec)[n/2:n]*2#/n # half of absolute spectra for orig signal
        spec_filt_abs = pl.absolute(spec_filt)[n/2:n]*2#/n # half of absolute spectra for filt signal
        freq_abs = freq[n/2:n] # half of absolute freqs axis for orig signal
        freq_filt_abs = freq_filt[n/2:n] # half of absolute freqs axis for filt signal

        spec_filt_abs = pl.ma.masked_equal(spec_filt_abs,0) # cutt off invalid values
        spec_filt_abs = pl.ma.compressed(spec_filt_abs)
        freq_filt_abs = pl.ma.masked_equal(freq_filt_abs,0)
        freq_filt_abs = pl.ma.compressed(freq_filt_abs)

        #spec_filt_abs *= 2/n # we cut off half of spectra - needs to be compensated
        #spec_abs *= 2/n

        spec_filt_log = 20*pl.log10(spec_filt_abs) # for output
        freq_filt_log = 10*pl.log10(freq_filt_abs)
    else:

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
    

    if preview and (result['frag'] == 0):
        fig = pl.figure("fft", figsize=(10, 10), facecolor='white') # plotting

        sp_sig = fig.add_subplot(221)
        sp_sig.set_title("Original signal")
        sp_sig.set_xlabel(r"Time, s")
        sp_sig.set_ylabel(r"HRV, ms")
        sp_sig.set_xlim(min(time),max(time))
        sp_sig.set_ylim(min(signal),max(signal))
        sp_sig.plot(time, signal)

        sp_spec = fig.add_subplot(222)
        sp_spec.set_title("Original signal - Fourier Specturm")
        sp_spec.set_xlabel(r"f, Hz")
        #sp_spec.set_xlim(min(freq_abs), max(freq_abs))
        sp_spec.set_ylabel(r"$\lg S(f)$")
        #sp_spec.set_ylim(min(spec_abs), max(spec_abs))
        sp_spec.loglog(freq_abs, spec_abs)
        
        sp_sig_filt = pl.subplot(223, sharex=sp_sig)
        sp_sig_filt.set_title("Filtered signal")
        sp_sig_filt.set_xlabel(r"Time, s")
        sp_sig_filt.set_ylabel(r"HRV, ms")
        #sp_sig_filt.set_ylim(min(signal_filt_uwed),max(signal_filt_uwed))
        sp_sig_filt.plot(time, signal_filt_uwed)
        
        sp_spec_filt = pl.subplot(224, sharex=sp_spec, sharey=sp_spec)
        sp_spec_filt.set_title("Filtered signal - Fourier Specturm")
        sp_spec_filt.set_xlabel(r"f, Hz")
        sp_spec_filt.set_ylabel(r"$\lg S(f)$")
        sp_spec_filt.loglog(freq_filt_abs, spec_filt_abs)
        
        pl.show()
        pl.savefig("bandpass_%(record)s_%(frag)s.png" % result, facecolor='w',
             dgecolor='k', transparent=True)
        pl.close()

    return freq_filt_log, spec_filt_log


def process_signal(record, annotator, diagnosis=None, slide_rate=.5,\
        preview=False):
    global config
    # rdhdr(record)
    info = rdhdr(record)
    if diagnosis:
        info['diagnosis'] = diagnosis
    print "Processing %s: %s %s %s" % (record, info['gender'], info['age'],\
        info['diagnosis'])

    # rdann(record, annotator, start=0, end=-1, types=[])
    #ann = rdann(record, annotator, start, end, [1])
    r_samples = read_r_samples(record, annotator)

    r_times = pl.empty(len(r_samples), dtype='float')
    r_times = r_samples / info['samp_freq']

    del r_samples

    #for i in range(10):
    #    print "{0:.10f}".format(r_times[i])
    
    time, hrv = variability(r_times)
    del r_times

    time_filt, hrv_filt = common.filter2d(time, hrv,\
        algos=config['HRV_FILTER_ALGO'])
    #print pl.shape(time_filt), pl.shape(hrv_filt)
    time_interp, hrv_interp = interpolate(time_filt, hrv_filt, config['SPLINE_ORDER'])

    if preview:
        plot_results.plot_hrv([time, hrv], [time_filt, hrv_filt], [time_interp, hrv_interp], record,\
            preview=preview)

    del time, hrv

    #signal_to_csv(record, time, hrv, info)

    window = config['INTERP_FREQ'] * config['WINDOW']

    if window % 2 != 0:
        warn("Window (%s) is not even. Adjust INTERP_FREQ and/or WINDOW parameters to be x**2" % (window, ))
    if not (window != 0 and ((window & (window - 1)) == 0)):
        warn("Window (%s) is not power of 2. Adjust INTERP_FREQ and/or WINDOW parameters" % (window, ))

    # pl.figure()
    # sp = pl.subplot(111)
    # Pxx, freq, bins, im = pl.specgram(hrv_interp, NFFT=window, Fs=config['INTERP_FREQ'], detrend=pl.detrend_linear, window=pl.window_none, pad_to=4*window, sides='onesided', noverlap=32)
    # sp.pcolor(bins, freq, Pxx)
    # sp.set_yscale('log')
    # #sp.set_xscale('log')
    # #pl.colorbar()
    # pl.show()
    # pl.close()

    frag_count = int((len(time_interp) - window) / window / slide_rate)
    results = []
    for frag in range(0, frag_count): # get WINDOW secs exactly
        result = {'record': record, 'frag': frag+1}
        #print "fragment %s of %s" % (frag+1, int(len(time_interp) / INTERP_FREQ / WINDOW))
        #chunk_first = ann[c][0] # get sample number of first RR interval
        #chunk_last = ann[c+w][0] # get sample number of last RR interval
        #for key in ['age','gender']:
        #  if key in info:
        #    r[key] = info[key]
        #  else:
        #    r[key] = ''

        frag_beg = frag * window * slide_rate
        frag_end = frag * window * slide_rate + window

        if "base_time" in info:
            r_first_full = info['base_time'] + \
                dt.timedelta(seconds=float(time_interp[frag_beg]))
            r_last_full = info['base_time'] + \
                dt.timedelta(seconds=float(time_interp[frag_end]))
        else:
            r_first_full = dt.datetime(1900, 01, 01, 10, 0, 0) + \
                dt.timedelta(seconds=float(time_interp[frag_beg]))
            r_last_full = dt.datetime(1900, 01, 01, 10, 0, 0) + \
                dt.timedelta(seconds=float(time_interp[frag_end]))
    
        result['time_from'] = r_first_full.strftime("%H:%M:%S")
        result['time_to'] = r_last_full.strftime("%H:%M:%S")

        time_frag = pl.empty(window, dtype='float')
        time_frag = time_interp[frag_beg: frag_end]
        time_frag -= time_frag.min()
        hrv_frag = pl.empty(window, dtype='float')
        hrv_frag = hrv_interp[frag_beg: frag_end]

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
        
        #result['freq'] = freq_filt
        #result['spec'] = fft_filt

        _a,_b = common.approximate(freq_filt,fft_filt)

        result['beta'] = -_a
        if result['frag'] % 100 == 0:
            line_y = [_a*x+_b for x in freq_filt]
            plot_results.plot_beta(freq_filt, fft_filt, line_y, result, preview=preview)
            del line_y

        #if beta_limit[0] <= result['beta'] and result['beta'] <= beta_limit[1]:
        results.append(result)
        print "%(frag)03d/%(frag_count)03d: %(time_from)s - %(time_to)s\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % dict(result.items()+{'frag_count':frag_count}.items())
        del frag_beg, frag_end, freq_filt, fft_filt

    if results:
        stats(results)
        results_to_csv(results, record, info)
        #pl.figure()
        #pl.yscale('symlog', linthreshy=0.01)
        #freqs = pl.fftfreq(window, 1.0/config['INTERP_FREQ'])[:window/2]
        #pl.pcolormesh(pl.arange(frag_count), freqs, [r['spec'] for r in results])
        #ax2.axis('tight')

        plot_results.plot_beta_std(results, preview=preview)
        plot_results.plot_beta_cv(results, preview=preview)
    else:
        warn("No results")

    return info, results

def read_r_samples(record, annotator):
    r_samples = []
    # Faster for 25% ann2rr -r 16265 -a atr -v s8 -i s8 -P N -p N
    proc = subprocess.Popen(['rdann','-r',record,'-a',annotator,'-p','N','-c','0'],bufsize=-1,stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        sample = line.split('  ')[1]
        if sample.isdigit():
            #time = (ann[0])[1:-1]
            r_samples.append(int(sample))

    return pl.array(r_samples, dtype=pl.uint32)


def results_to_csv(results,record,info):
    outfile = open("results_%s.csv" % (record,), "w")
    outfile.write("%s,%s,%s,%s\n" % (record, info['gender'], info['age'], info['diagnosis']))
    for r in results:
        outfile.write("%(time_from)s,%(time_to)s,%(beta)s,%(std)s,%(cov)s,%(mean)s\n" % r)
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
    print "Min:\t\t\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % mins
    print "Max:\t\t\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % maxs
    print "Mean:\t\t\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % means
    print "Median:\t\t\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % medians
    print "Stdev:\t\t\t\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % stdevs


if __name__ == '__main__':
    global config
    config = common.load_config()
    #pprint(config)
    batch = False
    if batch:
        for diag in config['SIGNALS']:
            for record in diag['records'].split():
                process_signal(record, diag['annotator'], diag['diagnosis'],
                    slide_rate=config['SLIDE_RATE'], preview=config['PREVIEW'])
    else:
        records = 'nsr001 nsr002 nsr003 nsr004 nsr005 nsr006 nsr007 nsr008 nsr009 nsr010 nsr011 nsr012 nsr013 nsr014 nsr015 nsr016 nsr017 nsr018 nsr019 nsr020 nsr021 nsr022 nsr023 nsr024 nsr025 nsr026 nsr027 nsr028 nsr029 nsr030 nsr031 nsr032 nsr033 nsr034 nsr035 nsr036 nsr037 nsr038 nsr039 nsr040 nsr041 nsr042 nsr043 nsr044 nsr045 nsr046 nsr047 nsr048 nsr049 nsr050 nsr051 nsr052 nsr053 nsr054'
        #records = '16483 16773 16795 17453'
        #signals = ['16273.atr', '16272.atr']
        #signals = ['chf04.ecg']
        #signals = ['chf05.ecg', 'chf06.ecg']
        #signals = ['nsr004.ecg']
        #for signal in signals:
        #process_signal('chf05','ecg', "CHF", preview=True)
        for record in records.split():
            process_signal(record, 'ecg', "NM2", config['SLIDE_RATE'], preview=False)
        #process_signal('16273','atr', "Normal", preview=True)
        #process_signal('nsr006','ecg', "Normal", preview=True)
        pass
    #save_params("config.txt")