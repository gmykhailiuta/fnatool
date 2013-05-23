#!/usr/bin/env python
import pylab as pl
from wfdbtools import rdhdr, rdann
from sys import exit
import datetime as dt
from scipy.interpolate import interp1d,splrep,splev
from pprint import pprint
from warnings import warn
import common
#from tools import draw

"""
+ Summary plot for all signals
  + beta_std by desease
  + beta_cv by desease
  + Boxed for every desease
+ Boxed plot for signal, all signals
+ Count stats for all signals
Different interpolations (esp. spline)
+ Sliding window
Compare window size
2-sigma cutting signals
Lomb FFT
Parameters as command Parameters
Do not use pylab (or fix case if using Gtk)
"""

WINDOW = 4096 # should be 2**x
INTERP_FREQ = 4 # should be 2**x
INTERP_KIND = 'linear' # linear/cubic
FREQ_LIMIT = (0.0033, 0.04) # (min_freq, max_freq)
BETA_LIMIT = (0.2, 1.8) # used to be 0.5 < beta < 1.5
VALID_RR_RATIO = 0.2 # 0 < x < 1
SLIDE_RATE = 1 # 0 < x < 1

def results_to_csv(results,record,info,diagnosis="?"):
  outfile = open("results_%s.csv" % (record,), "w")
  if 'diagnosis' in info:
    diagnosis = info['diagnosis']

  outfile.write("%s,%s,%s,%s\n" % (record, info['gender'], info['age'], diagnosis))
  for r in results:
    outfile.write("%(time_from)s,%(time_to)s,%(beta)s,%(std)s,%(cov)s,%(mean)s\n" % r)
  outfile.close()


def fft_filter(time, signal, result, samp_freq=INTERP_FREQ, freq_limit = FREQ_LIMIT, preview=False):
  """
   Applies an Band Pass filter to signal in the frequency domain and plots
      signals and their spectrals.
   In:  time : ndarray, relative RR times vector
        signal : ndarray, HRV vector
        result : dict, current fragment info
        freq_limit : frequencies for band pass filter, analyzed freq range
   Out: 20*log10 of frequency vector
        20*log10 of specturm
  """
  time -= min(time) # make relative time for this fragment
  n = len(signal)

  window = pl.hanning(n) # window
  signal_wed = signal * window

  spec = pl.fftshift(pl.fft(signal_wed)) # entering to frequency domain
  spec_filt = pl.zeros(spec.shape,dtype=complex)
  # fftshift moves zero-frequency component 
  # to the center of the array
  freq = pl.fftshift(pl.fftfreq(n, 1.0/samp_freq)) # get freqs axis values
  freq_filt = pl.zeros(freq.shape,dtype=float) # same for filtered

  for i in range(n):  # filter by frequency
    if freq_limit[0] < abs(freq[i]) and abs(freq[i]) < freq_limit[1]:
      freq_filt[i] = freq[i]
      spec_filt[i] = spec[i]
    else:
      freq_filt[i] = 0 # fill invalid frequencies with small value
      spec_filt[i] = 0

  signal_filt = pl.absolute(pl.ifft(pl.ifftshift(spec_filt))) # get filtered signal
  signal_filt_uwed = signal_filt * window * (1-window)
  #signal_filt_uwed = signal_filt_uwed[:len(time)]
  #signal = signal[:len(time)]
 
  spec_abs = pl.absolute(spec)[n/2:n]*2/n # half of absolute spectra for orig signal
  spec_filt_abs = pl.absolute(spec_filt)[n/2:n]*2/n # half of absolute spectra for filt signal
  freq_abs = freq[n/2:n] # half of absolute freqs axis for orig signal
  freq_filt_abs = freq_filt[n/2:n] # half of absolute freqs axis for filt signal

  spec_filt_abs = pl.ma.masked_equal(spec_filt_abs,0) # cutt off invalid values
  spec_filt_abs = pl.ma.compressed(spec_filt_abs)
  freq_filt_abs = pl.ma.masked_equal(freq_filt_abs,0)
  freq_filt_abs = pl.ma.compressed(freq_filt_abs)

  #spec_filt_abs *= 2/n # we cut off half of spectra - needs to be compensated
  #spec_abs *= 2/n

  spec_filt_log = 20*pl.log10(spec_filt_abs) # for output
  freq_filt_log = 20*pl.log10(freq_filt_abs)

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
  
  if preview and result['frag'] == 1:
     pl.show()
  pl.savefig("bandpass_%(record)s_%(frag)s.png" % result, facecolor='w', edgecolor='k', transparent=True)
  pl.close()

  return freq_filt_log, spec_filt_log

def approximate(x,y):
    """
    Linear approximation of y=f(x)
    In:
      x : ndarray
      y : ndarray
    Out:
      a, b : float, as in a*x+b=y
    """
    A = pl.vstack([x, pl.ones(len(x))]).T
    a, b = pl.lstsq(A, y)[0]
    return a, b


def plot_beta_cv(results, preview=False):
  x = [x['cov'] for x in results]
  y = [y['beta'] for y in results]
  pl.figure(figsize=(6, 6), facecolor='white')
  pl.title(r'%(record)s: $\sigma/\bar{x}$' % results[0])
  pl.plot(x,y,'.r')
  pl.axhline(y=1,color='b')
  pl.ylim(0,2)
  pl.xlabel(r"$\sigma/\bar{x}$")
  pl.ylabel(r"$\beta$")
  pl.savefig("cv_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
  if preview:  
    pl.show()
  pl.close()


def plot_beta_std(results, preview=False):
  x = [x['std'] for x in results]
  y = [y['beta'] for y in results]
  pl.figure(figsize=(6, 6), facecolor='white')
  pl.title(r'%(record)s: $\beta/\sigma$' % results[0])
  pl.plot(x,y,'.r')
  pl.ylim(0,2)
  pl.xlim(0,140)
  pl.axhline(y=1,color='b')
  pl.axvline(x=70,color='b')
  pl.xlabel(r"$\sigma$, ms")
  pl.ylabel(r"$\beta$")
  pl.savefig("std_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
  if preview:
    pl.show()
  pl.close()


def plot_beta(freq, fft, aprox_y, result, preview=False):
  pl.figure(figsize=(6, 6),facecolor='white')
  pl.plot(freq,fft,'b')
  pl.plot(freq,aprox_y,'r')
  pl.title('%(record)s: %(time_from)s - %(time_to)s' % result)
  pl.xlim(min(freq),max(freq))
  pl.xlabel(r"$\lg f$")
  pl.ylabel(r"$\lg S$")
  pl.savefig("beta_%(record)s_%(frag)s.png" % result,facecolor='w',edgecolor='k',transparent=True)
  if preview and result['frag'] == 1:
    pl.show()
  pl.close()


def variability(r_times, valid_rr_ratio=VALID_RR_RATIO):
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
  rrs = pl.zeros(len(r_times), dtype='float32')
  time = pl.array(r_times, dtype='float32')
  #last_r_time = 0
  last_rr = 0
  time[0] = 0
  for i in range(1, len(r_times)):
    rr = (r_times[i]-r_times[i-1])* 1000
    if last_rr: # if first rr computed, go with validation
      if (abs(last_rr / rr - 1) <= valid_rr_ratio) and (rr > 100) and (rr < 2000): # current rr differs less then 20% of previous one
        rrs[i] = rr
        last_rr = rr
      else:
        rrs[i] = 0
        time[i] = 0
    else: # no rrs computed yet. need first for validation to work
      rrs[i] = rr
      last_rr = rr
    #last_rr = rr


  rrs = pl.ma.masked_equal(rrs,0)
  rrs = pl.ma.compressed(rrs)
  time = pl.ma.masked_equal(time,0)
  time = pl.ma.compressed(time)
  #rrs *= 1000
  print "NB: Deleted %0.2f%% intervals" % (float(len(r_times)-len(rrs))/len(r_times)*100,)
  return time, rrs

def plot_hrv(hrv, hrv_interp, preview=False):
    """
    Plot HRV and it's approximation.
    In:
        hrv : list, [time, hrv], ndarray,ndarray
            Time & hrv vectors
        hrv_interp : list, [time_interp, hrv_interp], ndarray,ndarray
            Interpolated time & hrv vectors
    """
    fig = pl.figure("signals")
    pl.plot(hrv[0], hrv[1], 'x')
    pl.ylabel('HRV (s)') 
    pl.plot(hrv_interp[0], hrv_interp[1], 'k',ls='dotted')
    pl.xlabel('Time (s)')
    if preview:
        pl.show()
    pl.close()


def interpolate(x, y, kind='linear'):
    """
    Interpolates x/y using cubic interpolation
    In:
      x, y : ndarray, Input data
    Out:
      xnew, ynew : ndarray, Interpolated data
    """
    xnew = pl.arange(x.min(), x.max(), 1.0/INTERP_FREQ, dtype='float32')
    if kind == 'linear':
      f = interp1d(x, y, kind='linear')
      ynew = f(xnew)
    elif kind == 'cubic':
      tck = splrep(x, y, s=0)
      ynew = splev(xnew, tck, der=0)
    else:
      raise BaseException("Interpolation kind not supported")
      return
    return xnew, ynew


def process_signal(record, annotator, diagnosis=None, start=0, end=-1, slide_rate=SLIDE_RATE, beta_limit=BETA_LIMIT, preview=False):
  # rdsamp(record, start=0, end=-1, interval=-1)
  info = rdhdr(record)
  if diagnosis:
    info['diagnosis'] = diagnosis
  print "Processing %s: %s %s %s" % (record, info['gender'], info['age'], info['diagnosis'])
  
  # rdann(record, annotator, start=0, end=-1, types=[])
  ann = rdann(record, annotator, start, end, [1])
  
  time, hrv = variability(ann[:,1]-ann[0,1])
  time_interp, hrv_interp = interpolate(time, hrv, INTERP_KIND)

  if preview:
    plot_hrv([time, hrv], [time_interp, hrv_interp], preview=preview)

  #signal_to_csv(record, time, hrv, info)
  
  window = INTERP_FREQ * WINDOW

  if window % 2 != 0:
    warn("Window is not even. Adjust INTERP_FREQ and/or WINDOW parameters to be x**2")
  if window and not window & (window - 1):
    warn("Window is not power of 2. Adjust INTERP_FREQ and/or WINDOW parameters")
  frag_count = int(len(time_interp) / INTERP_FREQ / WINDOW)
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

    frag_beg = frag * slide_rate * window
    frag_end = frag * slide_rate * window + window
    #print frag_beg, frag_end, len(time_interp)

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

    result['std'] = (hrv_interp[frag_beg: frag_end]).std()
    result['mean'] = (hrv_interp[frag_beg: frag_end]).mean()
    result['cov'] = result['std']/result['mean']*100

    freq_filt, fft_filt = fft_filter(time_interp[frag_beg: frag_end], \
        hrv_interp[frag_beg: frag_end], result)

    _a,_b = approximate(freq_filt,fft_filt)
    result['beta'] = -_a
    line_y = [_a*x+_b for x in freq_filt]
    plot_beta(freq_filt, fft_filt, line_y, result, preview=preview)

    if beta_limit[0] <= result['beta'] and result['beta'] <= beta_limit[1]:
      results.append(result)
      print "%(frag)02d/%(frag_count)02d: %(time_from)s - %(time_to)s\t%(beta)0.2f\t%(std)d\t%(cov)0.2f%%\t%(mean)d" % dict(result.items()+{'frag_count':frag_count}.items())
  if results:
    stats(results)
    results_to_csv(results, record, info, preview=preview)
    plot_beta_std(results, preview=preview)
    plot_beta_cv(results, preview=preview)
  else:
    warn("No results")

  return info, results


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
  batch = False
  if batch:
    diagnosis = 'Hypertension'
    for record in common.SIGNALS[diagnosis].split('  '):
      process_signal(record, diagnosis)
  else:
    #signals = ['16483.atr', '16773.atr']
    #signals = ['16273.atr', '16272.atr']
    #signals = ['chf04.ecg']
    #signals = ['chf05.ecg', 'chf06.ecg']
    #signals = ['nsr004.ecg']
    #for signal in signals:
    process_signal('chf05','ecg', "CHF", preview=True)