#!/usr/bin/env python
import pylab
import numpy
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
from sys import exit
from math import log
from scipy import signal
import anfft

def write_csv(records,record):
	outfile = open(record+".csv", "w")
	for r in records:
		line = "%s,%s,%s,%s,%s,%s,%s\n" % (record,r['time_from'],r['time_to'],r['betta'],r['std'],r['stn'],r['mean'])
		outfile.write(line)
	outfile.close()

def autocor(data):
	print "Autocorrelation"
	data_length = len(data)
	in2 = numpy.zeros(data_length * 2)

	in2[data_length/2:data_length/2+data_length] = data # This works for data_length being even

	# Do an array flipped convolution, which is a correlation.
	cor = signal.fftconvolve(in2, data[::-1], mode='valid') 
	cor = cor[len(cor)/2:]
	return cor

def fft(x,samp_freq=250):
	print "FFT"
	arr = numpy.zeros(2**(round(log(len(x),2))+3))
	arr[:len(x)] = x
	F = abs(anfft.fft(arr))
#	F = abs(pylab.fftpack.fft(arr))
	f = pylab.fftfreq(len(arr), 1/samp_freq)
	F = F[:len(F)/2]
	f = f[:len(f)/2]
	return F, f

def cut_freq(s,f,fmin=0.003,fmax=0.04):
	print "Cutting freqs %s - %s" % (fmin, fmax)
	imin = f.searchsorted(fmin)
	imax = f.searchsorted(fmax)
	return s[imin:imax], f[imin:imax]

def betta(x,y):
	A = pylab.vstack([x, pylab.ones(len(x))]).T
	a, b = pylab.lstsq(A, y)[0]
	_b = -a
	print "betta = %s" % (_b, )	
	y = a*x+b
#	fig(x,y)
	return _b

def fig(in1,in2=None,label="",show=False):	
	pylab.figure()
	if in2 is None:
		pylab.plot(in1,label=label)
	else:
		pylab.plot(in1,in2,label=label)
	if show:
		pylab.show()
		exit(0)

def process(record, write_file=True):
	rr = []

	# Read in the data from 0 to 10 seconds
	# rdsamp(record, start=0, end=-1, interval=-1)
	data, info = rdsamp(record, 0)
	pprint(info)

	print "total time ", int(info['samp_count'])/int(info['samp_freq'])

# rdann(record, annotator, start=0, end=-1, types=[])
	ann = rdann(record, 'atr', 0)
# annotation time in samples from start
	ann_x = (ann[:, 0] - data[0, 0]).astype('int')

#ann_y = v[ann_y,2]
#plot_data(data, info, ann)

	c = 0 # window first RR
	w = 1000 # window last RR

	print "RR count", len(ann)

	while c+w < len(ann):
		r = {}
	# get 10 RRs
		print "RR ",c," - ",c+w
		chunk_first = ann[c][0] # get sample number of first RR interval
		chunk_last = ann[c+w][0] # get sample number of last RR interval
	#time_interval = (chunk_last-chunk_first)/info['samp_freq']
		r['time_from'] = data[chunk_first,1]
		r['time_to'] = data[chunk_last,1]
		print "interval ", r['time_from'], " - ", r['time_to'], " = ", (chunk_last-chunk_first)/info['samp_freq'], " s"

		t = data[chunk_first:chunk_last,1] # time array of target window
		v = data[chunk_first:chunk_last,2] # data value of target window

#		fig(t,v,show=True)
		r['std'] = pylab.std(v)
		r['mean'] = pylab.mean(v)
		r['stn'] = r['std']/r['mean']
		
#		import pywt
#		cA, cD = pywt.dwt(v, 'coif3')
#		fig(cA)
#		fig(cD, show=True)

		cor = autocor(v)

		Ff, ff = fft(cor, info['samp_freq'])
		Fc, fc = cut_freq(Ff,ff)

		Fc = pylab.log(Fc)
		fc = pylab.log(fc)
#		fig(fc,Fc,show=False)

		r['betta'] = betta(fc,Fc)
#		fig(fc,Fc,show=True)
	
		pprint(r)
		rr.append(r)
		c +=w

	if write_file:
		write_csv(rr,record)

if __name__ == '__main__':
	process('s20171')
