#!/usr/bin/env python
#import scipy
#import scipy.fftpack
import pylab
import numpy
#import wfdbtools
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
from sys import exit
#from scipy.linalg import lstsq
from math import atan, log, pi
from scipy import signal
#from numpy import correlate, std, mean, ones

record = 's20011'
#outfile = open(record+".csv", "a+")

def kalman(z):
	"""Apply kalman filter to z"""
	# intial parameters
	n_iter = len(z)
	sz = (n_iter,) # size of array
	Q = 1e-5 # process variance

	# allocate space for arrays
	xhat=numpy.zeros(sz)      # a posteri estimate of x
	P=numpy.zeros(sz)         # a posteri error estimate
	xhatminus=numpy.zeros(sz) # a priori estimate of x
	Pminus=numpy.zeros(sz)    # a priori error estimate
	K=numpy.zeros(sz)         # gain or blending factor

	R = 0.1**16 # estimate of measurement variance, change to see effect

	# intial guesses
	xhat[0] = 0.0
	P[0] = 1.0

	for k in range(1,n_iter):
		# time update
		xhatminus[k] = xhat[k-1]
		Pminus[k] = P[k-1]+Q

		# measurement update
		K[k] = Pminus[k]/( Pminus[k]+R )
		xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
		P[k] = (1-K[k])*Pminus[k]
	return xhat

def autocor(data):
	print "autocorrelation"

	data_length = len(data)
	in2 = numpy.zeros(data_length * 2)

	in2[data_length/2:data_length/2+data_length] = data # This works for data_length being even

	# Do an array flipped convolution, which is a correlation.
	cor = signal.fftconvolve(in2, data[::-1], mode='valid') 
	cor = cor[len(cor)/2:]
#	pylab.figure()
#	pylab.plot(cor)
#pylab.show()
#exit(0)
	return cor

#corrfft = lambda x,y : pylab.irfft(pylab.rfft(x)*pylab.rfft(y[::-1]))

#	cor = pylab.correlate(v,v,mode='full')
#	cor = cor[len(cor)/2:]
#	cor = xhat
	#print "cor_len=",len(cor)
#	pylab.figure()
#	pylab.plot(cor,'r',label='correlation2')


def cut_freq(spower,freqs,freq_from=0.003,freq_to=0.04):
	print "Cutting freqs %s - %s" % (freq_from, freq_to)
	f_min, f_max = freqs.min(), freqs.min()
	i_min, i_max, i = 0,0,0
	for f in freqs:
		if f <= freq_from:
			f_min = f
			i_min = i
		if f <= freq_to:
			f_max = f
			i_max = i
		if f > freq_to:
			break
		i+=1

	#print "freqs " ,f_min, " - ", f_max
	#fmin = next(x for x in freqs if x >= 0.003)
	#fmax = next(x for x in freqs if x >= 0.04)

	spower = spower[i_min:i_max]
	#print F[:10]
	freqs = freqs[i_min:i_max] # right freqs only
	#print freqs[:10]

	#print len(v), len(F)
	#pprint(F[:10])
	#pprint(freqs[:10])
	return spower, freqs

# Read in the data from 0 to 10 seconds
# rdsamp(record, start=0, end=-1, interval=-1)
data, info = rdsamp(record, 0, 1000)
pprint(info)

print "total time ", int(info['samp_count'])/int(info['samp_freq'])

# rdann(record, annotator, start=0, end=-1, types=[])
ann = rdann(record, 'atr', 0, 1000)
# annotation time in samples from start
ann_x = (ann[:, 0] - data[0, 0]).astype('int')

#ann_y = v[ann_y,2]
#plot_data(data, info, ann)

c = 0 # window first RR
w = 1000 # window last RR

print "RR count", len(ann)

while c+w < len(ann):
	# get 10 RRs
	print "RR ",c," - ",c+w
	chunk_first = ann[c][0] # get sample number of first RR interval
	chunk_last = ann[c+w][0] # get sample number of last RR interval
	#time_interval = (chunk_last-chunk_first)/info['samp_freq']
	print "interval ", data[chunk_first,1], " - ", data[chunk_last,1], " = ", (chunk_last-chunk_first)/info['samp_freq'], " s"

	t = data[chunk_first:chunk_last,1] # time array of target window
	v = data[chunk_first:chunk_last,2] # data value of target window

#	pylab.figure()
#	pylab.plot(t,v)

	sigma = pylab.std(v)
	v_mean = pylab.mean(v)
	sdn = sigma/v_mean
	print "std=",sigma
	print "sdn=", sdn
	print "mean_RR=",v_mean


#	print "correlation"
#	cor = pylab.correlate(v,v,mode='full')
#	cor = cor[len(cor)/2:]
	#print "cor_len=",len(cor)
#	pylab.figure()
#	pylab.plot(cor,'r+',label='correlation')

#	print "fft"
#	F = abs(pylab.fftpack.rfft(cor,len(cor)*10))
#	freqs = pylab.fftfreq(len(cor)*10, 1/info['samp_freq'])
#	pylab.figure()
#	pylab.plot(F,'k+',label='fft2')

	print "len v = %s" % (len(v),)
	cor = autocor(v)

	print "fft"
	F = abs(pylab.fftpack.fft(cor, len(cor)*2**3))
#	F = F[len(F)/2+1:]
	freqs = pylab.fftfreq(len(cor)*2**3, 1/info['samp_freq'])
#	freqs = freqs[len(freqs)/2:]
	print "len_f=",len(freqs),"len_F=",len(F)
#	pylab.figure()
#	pylab.plot(freqs,F,'k',label='fft22')
#	pylab.show()
#	exit(0)


	F, freqs = cut_freq(F,freqs)
#	F2, f2 = cut_freq(F,freqs,0,1)
#	pylab.figure()
#	pylab.plot(f2,F2,'k')

	F = pylab.log(F)
	freqs = pylab.log(freqs)

#	pylab.figure()
#        pylab.plot(freqs, F)
#	pylab.show()
#	exit(0)

	A = pylab.vstack([freqs, pylab.ones(len(freqs))]).T
	#print "A=",A
	a, b = pylab.lstsq(A, F)[0]
	#print "a=",a," b=",b
	betta = -a
	print "betta=",betta
	#print "alpha=",atan(-a)/pi*180
	
	y = a*freqs+b
#	pylab.plot(freqs, y, 'r', label='Fitted line')

#	pylab.show()
#	exit(0)
	
	ostr=record+","+str(data[chunk_first,1])+","+str(data[chunk_last,1])+","+str(betta)+","+str(sigma)+","+str(sdn)+","+str(v_mean)+"\n"
	print ostr
	outfile.write(ostr)
	
	c +=w 
#	break


outfile.close()

