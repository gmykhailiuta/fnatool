#!/usr/bin/env python
import pylab as pl
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
from sys import exit
#import math
from scipy.signal import fftconvolve
import anfft
import datetime as dt
#import time

def write_csv(records,record,info):
  outfile = open(record+".csv", "w")
  if 'diagnosis' in info:
    d = info['diagnosis']
  outfile.write("%s,%s\n" % (record, d))
  for r in records:
    line = "%s,%s,%s,%s,%s,%s,%s,%s\n" % (r['time_from'],r['time_to'],r['beta'],r['std'],r['sdn'],r['mean'],r['gender'],r['age'])
    outfile.write(line)
  outfile.close()

def sig2csv(record,t,v,info):
#  t = [row[1] for row in ann]
#  v = delta(t)*1000
#  t.pop()
  v = v*1000
  outfile = open("rr_"+record+".csv", "w")
  for i in range(len(v)):
    t_full = info['base_time'] + dt.timedelta(seconds=t[i])
    line = "%s,%s,%s\n" % (i+1,t_full.strftime("%H:%M:%S.%f"),v[i])
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

def fft(x,samp_freq=250):
  print "FFT"
  arr = pl.zeros(2**(round(pylab.log2(len(x)))+8))
  arr[:len(x)] = x
  F = abs(anfft.fft(arr))
#  F = abs(pl.fftpack.fft(arr))
  f = pl.fftfreq(len(arr), 1/samp_freq)
  F = F[:len(F)/2]
  f = f[:len(f)/2]
  return F, f

def fft2(x,sample_spacing):
  print "FFT"
  arr = pl.zeros(2**(round(pylab.log2(len(x)))+8))
  arr[:len(x)] = x
  F = abs(anfft.fft(arr))
#  F = abs(pl.fftpack.fft(arr))
  f = pl.fftfreq(len(arr), sample_spacing)
  F = F[:len(F)/2]
  f = f[:len(f)/2]
  return F, f

def cut_freq(s,f,fmin=0.004,fmax=0.4):
  print "Cutting freqs %s - %s" % (fmin, fmax)
  imin = f.searchsorted(fmin)
  imax = f.searchsorted(fmax)
  return s[imin:imax], f[imin:imax]

def cut_freq0(s,f,fmin=0.004,fmax=0.4):
  print "Cutting freqs %s - %s" % (fmin, fmax)
  imin = f.searchsorted(fmin)
  imax = f.searchsorted(fmax)
  s[:imin] = 0
  f[:imin] = 0
  s[imax:] = 0
  f[imax:] = 0
  return s, f

def aprox(x,l  A = pylab.vstack([x, pylab.ones(len(x))]).T
  a, b = pl.lstsq(A, y)[0]
  return a, b

#def beta(x,y):
#  a, b = aprox(x,y)
#  print "beta = %s" % (-_a, )
#  return -a

def fig(in1,in2=None,label="",show=False):  
  pl.figure(facecolor='white')
  if in2 is None:
    pl.plot(in1,label=label)
  else:
    pl.plot(in1,in2,label=label)
  if show:
    pl.show()
    exit(0)

def draw_sdn(rr, fname):
  x = [x['sdn'] for x in rr]
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
  x = [x['std']*1000 for x in rr]
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
  pl.xlabel(r"$\log{f}$")
  pl.ylabel(r"$\log{S}$")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
#  pl.close()

def draw_sig(fname, x, y):
  pl.figure(figsize=(100, 6), facecolor='white')
  pl.plot(x,y,'k')
  pl.xlabel(r"t, s")
  pl.ylabel(r"t, s")
  pl.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
  pl.close()

def delta(data):
  d = pl.zeros(len(data)-1)
  for i in range(1,len(data)):
    d[i-1] = data[i]-data[i-1]
  return d

def process(record, annotator="atr", end=-1):
  print "Processing %s" % (record,)

  # Read in the data from 0 to 10 seconds
  # rdsamp(record, start=0, end=-1, interval=-1)
  data, info = rdsamp(record, 0, 10)
  pprint(info)
#  pprint(data[:100])

  print "total time ", int(info['samp_count'])/int(info['samp_freq'])

# rdann(record, annotator, start=0, end=-1, types=[])
  ann = rdann(record, annotator, 0, end)
#  print ann[:100]
# annotation time in samples from start
  ann_x = (ann[:, 0] - data[0, 0]).astype('int')


  T = [row[1] for row in ann]
  V = delta(T)
  T.pop()
  sig2csv(record,T,V,info)
  draw_sig(record+".png",T,V)
  
  ## write data
#  write_csv_data(record,ann,info)
  ## / write data

#ann_y = v[ann_y,2]
#  plot_data(data, info, ann)
#  exit(0)

  w = 1000 # window last RR

  print "RR count", len(ann)

  c = 0 # window first RR
  rr = []
  while c+w < len(T):
    r = {}
  # get 10 RRs
    print "RR %s - %s / %s" % (c,c+w,len(ann))
    #chunk_first = ann[c][0] # get sample number of first RR interval
    #chunk_last = ann[c+w][0] # get sample number of last RR interval
    for key in ['age','gender']:
      if key in info:
        r[key] = info[key]
      else:
        r[key] = ''

    r_first = T[c]
    r_last = T[c+w]
    print "interval %s - %s = %s s" % (r_first, r_last, r_last - r_first)
    if "base_time" in info:
      r_first_full = info['base_time'] + dt.timedelta(seconds=r_first)
      r_last_full = info['base_time'] + dt.timedelta(seconds=r_last)
      r['time_from'] = r_first_full.strftime("%H:%M:%S.%f")
      r['time_to'] = r_last_full.strftime("%H:%M:%S.%f")
    else:
      r['time_from'] = str(dt.timedelta(0,r_first))
      r['time_to'] = str(dt.timedelta(0,r_last))
    
     t = T[c:c+w]
    v = V[c:c+w]

#    fig(t,v,show=True)
#    write_csv_data(record,t,v)
    r['std'] = pl.std(v)
    r['mean'] = pl.mean(v)
    r['sdn'] = r['std']/r['mean']

#    vx = v/max(v) - r['mean']
#    vx = v - min(v)
#    vx = vx / max(v)
#    v = vx

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
    cor = v
#    cor = pl.correlate(v,v, mode="full")
#      cor = pl.correlate(v,v)

#      print "Correlation took %s" % (time.time()-t1, )

#      exit(0)
#    fig(cor)

    
    #Ff, ff = fft2(cor, info['samp_freq'])
    Ff, ff = fft2(v, (r_last - r_first)/1000)

    fig(ff,Ff,show=True)
#    fig(v)
#    F = abs(pl.fft(v))
#    
#    f = pl.fftfreq(len(v))
#    fig(F)
#    Fc, fc = cut_freq0(Ff,ff)
#
#    fig(fc,Fc,show=True)
#    del cor, v

    Fc, fc = cut_freq(Ff,ff)
    fig(fc,Fc,show=False)

#    v_filt = abs(anfft.ifft(Fc))
#    fig(v_filt,show=True)
    del Ff, ff
    
    Fc = pl.log10(Fc)
    fc = pl.log10(fc)
    fig(fc,Fc,show=False)

    _l = aprox(fc,Fc)
    r['beta'] = -_a
    draw_flicker(fc,Fc,_a*fc+_b,"%s_flicker_i%s.png" % (record, int(c/w)))
    #exit(0)

    pylab.show()
    exit(0)


    pprint(r)
    rr.append(r)
    c +=w

  draw_std(rr, "%s_std.png" % (record, ))
  draw_sdn(rr, "%s_sdn.png" % (record, ))

  write_csv(rr, record, info)
  del rr, data, info, ann, ann_x  




if __name__ == '__main__':
  record = 's20011'
  annotator = 'atr'
  process(record, annotator, 1010)
