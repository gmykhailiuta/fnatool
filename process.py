#!/usr/bin/env python
import pylab as pl
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
from sys import exit
#import math
from scipy.signal import fftconvolve
#import anfft
import datetime as dt
import numpy as np
#import time

WINDOW = 1024
INTERP_FREQ = 16

def write_csv(records,record,info):
  outfile = open(record+".csv", "w")
  if 'diagnosis' in info:
    d = info['diagnosis']
  outfile.write("%s,%s,%s,%s\n" % (record, info['gender'], info['age'], d))
  for r in records:
    outfile.write("%(time_from)s,%(time_to)s,%(beta)s,%(std)s,%(cov)s,%(mean)s\n" % r)
  outfile.close()

def sig2csv(record,t,v,info):
#  t = [row[1] for row in ann]
#  v = delta(t)*1000
#  t.pop()
  v = v*1000
  outfile = open("rr_"+record+".csv", "w")
  for i in range(len(v)):
    t_full = info['base_time'] + dt.timedelta(seconds=t[i])
    line = "%s %s %s\n" % (i+1,t_full.strftime("%H:%M:%S.%f"),v[i])
    outfile.write(line)
  outfile.close()

def autocor(data):
  print "Autocorrelation"
  data_length = len(data)
  in2 = pl.zeros(data_length * 2)

  in2[data_length/2:data_length/2+data_length] = data # This works for data_length being even

  # Do an array flipped convolution, which is a correlation.
  cor = fftconvolve(in2, data[::-1], mode='valid') 
  cor = cor[len(cor)/2:]
  return cor

def bandPass(T, V, fname):
  """
   Applies an Band Pass filter to V in the frequency domain.
   In:  T, list with time axis
        V, list with RR length signal
        #TD, list with time axis for ECG signal
        #D, list with data signal
        fname, file name to save figure to
   Out: Frequency in lg scale
        Specturm power in lg scale
  """
  n = len(V)
  #Ts = pl.mean(V)
  T -= min(T)
  #t = pl.arange(0,sum(V),float(Ts))
  fig = pl.figure(figsize=(10, 10), facecolor='white')
  def onclick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata)
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  def on_key(event):
    print('you pressed', event.key, event.xdata, event.ydata)

  cid = fig.canvas.mpl_connect('key_press_event', on_key)

  #pl.subplot(2,2,1)
  #print len(TD)
  #print len(D)
  #pl.plot(TD, D)


  pl.subplot(2,2,1)
  pl.title("Original signal")
  pl.xlabel(r"Time, s")
  pl.ylabel(r"RR, s")
  pl.xlim(min(T),max(T))
  pl.ylim(min(V),max(V))
  pl.plot(T, V)
  #Vzeroed = pl.zeros(2**(round(pl.log2(len(V)))+0))
  #Vzeroed[:len(V)] = V
  #while len(T) < len(Vzeroed):
  #  T = pl.append(T,T[len(T)-1]+Ts)
  #pl.plot(T, Vzeroed)
  I = pl.fftshift(pl.fft(V)) # entering to frequency domain
  ffs = pl.fftshift(pl.fftfreq(len(V), INTERP_FREQ))
  # fftshift moves zero-frequency component 
  # to the center of the array
  #P = pl.zeros(2**(round(pl.log2(len(I)))+8), dtype=complex)
  #n = len(P)
  P = pl.zeros(I.shape,dtype=complex)
  pl.subplot(2,2,2)
  pl.title("Original signal - Fourier Specturm")
  pl.xlabel(r"f, Hz")
  pl.ylabel(r"$\lg S(f)$")
  F = pl.log10(abs(I))
  pl.xlim(min(ffs),max(ffs))
  pl.ylim(min(F),max(F))
  pl.plot(ffs,F)

  for i in range(n):  # frequency cutting
    if 0.003 < abs(ffs[i]) and abs(ffs[i]) < 0.04:
    #if abs(ffs[i]) < 0.4:
      P[i] = I[i]

  Vfilt = np.real(pl.ifft(pl.ifftshift(P)))
  pl.subplot(2,2,3)
  pl.title("Filtered signal")
  pl.xlabel(r"Time, s")
  pl.ylabel(r"RR, s")
  pl.xlim(min(T),max(T))
  pl.ylim(min(Vfilt),max(Vfilt))
  pl.plot(T, Vfilt)
  pl.subplot(2,2,4)
  pl.title("Filtered signal - Fourier Specturm")
  pl.xlabel(r"f, Hz")
  pl.ylabel(r"$\lg S(f)$")
  Ffilt = pl.log10(abs(P))
  pl.xlim(min(ffs),max(ffs))
  pl.ylim(min(F),max(F))
  pl.plot(ffs,Ffilt)
  #pl.figure()
  #pl.plot(pl.log10(ffs[range(n/2,n)]),pl.log10(abs(P))[range(n/2,n)])
  
  #pl.show()
 # exit(0)
  
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
  pl.close()
  flog = np.real(pl.log10(ffs[range(n/2,n)]))
  Flog = np.real(pl.log10(abs(2*P)/len(P))[range(n/2,n)])
  flog_ = list()
  Flog_ = list()
  for i in range(len(flog)):
    #if np.isinf(flog[i]):
    #  flog[i] = 10**-10
    #if np.isinf(Flog[i]):
    #  Flog[i] = 10**-10

    if (not np.isinf(flog[i])) and (not np.isinf(Flog[i])):
      flog_.append(flog[i])
      Flog_.append(Flog[i])
  #fig(flog_, Flog_, show=1)
  del flog, Flog
  return Vfilt,flog_,Flog_

def aprox(x,y):
  A = pl.vstack([x, pl.ones(len(x))]).T
  a, b = pl.lstsq(A, y)[0]
  return a, b

def fig(in1,in2=None,label="",show=False):  
  pl.figure(facecolor='white')
  if in2 is None:
    pl.plot(in1,label=label)
  else:
    pl.plot(in1,in2,label=label)
  if show:
    pl.show()
    exit(0)

def draw_cov(rr, fname):
  x = [x['cov'] for x in rr]
  y = [y['beta'] for y in rr]
  pl.figure(figsize=(6, 6), facecolor='white')
  pl.plot(x,y,'.r')
  pl.axhline(y=1,color='b')
  pl.ylim(0,2)
  pl.xlabel(r"$\sigma/\bar{x}$")
  pl.ylabel(r"$\beta$")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
#  pl.show()
#  exit(0)
  pl.close()


def draw_std(rr, fname):
  x = [x['std'] for x in rr]
  y = [y['beta'] for y in rr]
  pl.figure(figsize=(6, 6), facecolor='white')
  pl.plot(x,y,'.r')
  pl.ylim(0,2)
  pl.xlim(0,140)
  pl.axhline(y=1,color='b')
  pl.axvline(x=70,color='b')
  pl.xlabel(r"$\sigma$, ms")
  pl.ylabel(r"$\beta$")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
  pl.close()

def draw_flicker(f,F,y, fname):
  pl.figure(figsize=(6, 6),facecolor='white')
  pl.plot(f,F,'b')
  pl.plot(f,y,'r')
  pl.xlim(min(f),max(f))
  pl.xlabel(r"$\lg f$")
  pl.ylabel(r"$\lg S$")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
  pl.close()

def draw_sig(fname, x, y):
  pl.figure(figsize=(100, 6), facecolor='white')
  pl.plot(x,y,'k')
  pl.xlabel(r"t, s")
  pl.ylabel(r"t, s")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
  pl.close()

def variability(r_times):
  rrs = pl.zeros(len(r_times), dtype='float32')
  time = pl.array(r_times, dtype='float32')
  last_r_time = 0
  last_rr = 0
  time[0] = 0
  for i in range(len(r_times)):
    if last_r_time: # do not count delta for 1st value
      rr = r_times[i]-last_r_time
      if abs(last_rr / rr - 1) <= 0.2: # current rr differs less then 20% of previous one
        rrs[i] = rr
      else:
        rrs[i] = 0
        time[i] = 0
      last_rr = rr       
    last_r_time = r_times[i]
  rrs = pl.ma.masked_equal(rrs,0)
  rrs = pl.ma.compressed(rrs)
  time = pl.ma.masked_equal(time,0)
  time = pl.ma.compressed(time)
  return time, rrs

"""Ideas:
interpolate rrs signal
fuck all: return to working version
get betas for sleeping|unsleeping periods - mark them "S"|"US"
compare betas depending on hour
"""

def plot_signals(data, info, hrv_time, hrv_data, hrv_interp_time, hrv_interp_data, ann=None):
    """
    Plot the signals HRV.

    Parameters
    ----------
    data : (N, 4) ndarray
         Output array from rdsamp.
    info : dict
         Header information as a dictionary.
         Output from rdsamp
    ann : (N, 2) ndarray, optional
         Output from rdann

    Returns
    -------
    None
    Matplotlib figure is plotted with the signals and annotations.
    """

    time = data[:, 1] #in seconds. use data[:, 0] to use sample no.
    sig1 = data[:, 2]
    sig2 = data[:, 3]
    
    pl.figure()
    pl.subplot(311)
    pl.plot(time, sig1, '-')
    pl.xticks([])
    pl.ylabel('%s (mV)' %(info['signal_names'][0]))
    
    pl.subplot(312)
    pl.plot(time, sig2, 'k')
    pl.xticks([])
    pl.ylabel('%s (mV)' %(info['signal_names'][1])) 
    

    if ann != None:
        # annotation time in samples from start
        ann_x = (ann[:, 0] - data[0, 0]).astype('int')
        pl.plot(ann[:, 1], data[ann_x, 3], 'xr')
        pl.subplot(311)
        pl.plot(ann[:, 1], data[ann_x, 2], 'xr')

    pl.subplot(313)
    pl.plot(hrv_time, hrv_data, '.')
    pl.xticks([])
    pl.ylabel('HRV (s)') 
    #pl.subplot(414)
    pl.plot(hrv_interp_time, hrv_interp_data, 'k-')

    pl.xlabel('Time (s)')

    pl.show()

def interpolate(x, y):
    """
    x : ndarray
    y : ndarray
        Input data
    returns
    xnew, ynew : ndarray, ndarray
        Interpolated data
    """
    from scipy import interpolate

    tck = interpolate.splrep(x, y, s=0)
    xnew = pl.arange(x.min(), x.max(), 1.0/INTERP_FREQ, dtype='float32')
    ynew = interpolate.splev(xnew,tck,der=0)
    print "old: %s %s", y.min(), y.max()
    print "new: %s %s", ynew.min(), ynew.max()
    return xnew, ynew

def process(record, annotator="atr", start=0, end=-1):
  print "Processing %s" % (record,)

  # Read in the data from 0 to 10 seconds
  # rdsamp(record, start=0, end=-1, interval=-1)
  data, info = rdsamp(record, start, end)
  pprint(info)
  #pprint(data[:100])

  print "total time ", int(info['samp_count'])/int(info['samp_freq'])

# rdann(record, annotator, start=0, end=-1, types=[])
  ann = rdann(record, annotator, start, end)
  #print ann[:100]
# annotation time in samples from start
  ann_x = (ann[:, 0] - data[0, 0]).astype('int')

  ann_time = [row[1] for row in ann]
  time, hrv = variability(ann_time)
  
  time_interp, hrv_interp = interpolate(time, hrv)

  #plot_signals(data, info, T, V, ann)
  #exit(0)

  #fig(T,V,show=1)
  #sig2csv(record,T,V,info)
  #draw_sig(record+".png",T,V)
  
  ## write data
#  write_csv_data(record,ann,info)
  ## / write data

#ann_y = v[ann_y,2]
#  plot_data(data, info, ann)
#  exit(0)

  print "RR count", len(ann)
  c = 0 # window first RR
  results = []
  for frag in range(0, int(len(time_interp) / INTERP_FREQ / WINDOW)): # get WINDOW secs exactly
    result = {}
    print "fragment %s of %s" % (frag+1, int(len(time_interp) / INTERP_FREQ / WINDOW))
    #chunk_first = ann[c][0] # get sample number of first RR interval
    #chunk_last = ann[c+w][0] # get sample number of last RR interval
    #for key in ['age','gender']:
    #  if key in info:
    #    r[key] = info[key]
    #  else:
    #    r[key] = ''

    frag_beg = frag * INTERP_FREQ * WINDOW
    frag_end = (frag+1) * INTERP_FREQ * WINDOW
    print "interval %s - %s = %s s" % (time_interp[frag_beg], time_interp[frag_end],
        time_interp[frag_end] - time_interp[frag_beg])
    if "base_time" in info:
      r_first_full = info['base_time'] + \
          dt.timedelta(seconds=float(time_interp[frag_beg]))
      r_last_full = info['base_time'] + \
          dt.timedelta(seconds=float(time_interp[frag_end]))
      result['time_from'] = r_first_full.strftime("%H:%M:%S.%f")
      result['time_to'] = r_last_full.strftime("%H:%M:%S.%f")
    else:
      result['time_from'] = str(dt.timedelta(0,time_interp[frag_beg]))
      result['time_to'] = str(dt.timedelta(0,time_interp[frag_end]))
    
    #print S[c]
    #print S[c+w]
    #print len(data)
    #d = data[S[c]:S[c+w],2]
    #td = data[S[c]:S[c+w],1]

    v_filt, fc, Fc = bandPass(time_interp[frag_beg: frag_end], \
        hrv_interp[frag_beg: frag_end], "%s_s%s.png" % (record, frag))
    
    result['std'] = (hrv[frag_beg: frag_end]).std()*1000
    result['mean'] = (hrv[frag_beg: frag_end]).mean()*1000
    result['cov'] = result['std']/result['mean']*100

#    fig(t,v)
#    import pywt
#    v = pl.append(v,[0.])
#    print len(v)
#    output = pywt.swt(v, 'db6', level=pywt.swt_max_level(len(v)))
#    output = pywt.thresholding.soft(outl r['std']*pylab.sqrt(2*pylab.log(len(v))))
#    output = pywt.thresholding.zero(output)
#    final = pywt.iswt(output, 'db6')
#    cA = pl.zeros(len(cA))
#    fig(cA)
#    fig(cD)
#    cD = pywt.thresholding.soft(cD,0.4)
#    vt = pywt.idwt(cA,cD,'db6')
#    res += v + final[:len(v)]
#    fig(res, show=True)

#    from denoise import denoise
#    J_dn=2
#    nm_dn='sym4'
#    mode_denoise='swt'
#    dn_method='visushrink'
#    dn_thresh='soft'
#    out=denoise(v,mode_denoise,nm_dn,J_dn,dn_method,dn_thresh)
#    fig(out, show=True)

#      t1=time.time()
#    cor = autocor(v)
#    cor = v
#    cor = pl.correlate(v,v, mode="full")
#      cor = pl.correlate(v,v)

#      print "Correlation took %s" % (time.time()-t1, )

#      exit(0)
#    fig(cor)

    
    #Ff, ff = fft2(cor, info['samp_freq'])
    #Ff, ff = fft2(v, (r_last - r_first)/1000)

    #fig(ff,Ff,show=True)
#    fig(v)
#    F = abs(pl.fft(v))
#    
#    f = pl.fftfreq(len(v))
#    fig(F)
#    Fc, fc = cut_freq0(Ff,ff)
#
#    fig(fc,Fc,show=True)
#    del cor, v

    #Fc, fc = cut_freq(Ff,ff)
    #fig(fc,Fc,show=False)

#    v_filt = abs(anfft.ifft(Fc))
#    fig(v_filt,show=True)
    #del Ff, ff
    
    #Fc = pl.log10(Fc)
    #fc = pl.log10(fc)
    #fig(fc,Fc,show=False)
    #print("aprox")
    _a,_b = aprox(fc,Fc)
    result['beta'] = -_a
    #print("done")
    y = [_a*x+_b for x in fc]
    draw_flicker(fc,Fc,y,"%s_flicker_i%s.png" % (record, frag))
    #exit(0)

#    pylab.show()
 #   exit(0)

    del fc, Fc, y, v_filt
    pprint(result)
    results.append(result)

  draw_std(results, "%s_std.png" % (record, ))
  draw_cov(results, "%s_cov.png" % (record, ))

  write_csv(results, record, info)
  del data, info, ann, results




if __name__ == '__main__':
  normal = True
  if normal:
    record = '16265'
    annotator = 'atr'
  else:
    record = 'chf01'
    annotator = 'ecg'
  process(record, annotator)
