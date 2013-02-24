#!/usr/bin/env python
import scipy
import scipy.fftpack
import pylab
#import numpy
#import wfdbtools
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
import sys
from scipy.linalg import lstsq
from math import atan, log, pi
from numpy import correlate, std, mean, ones

record = 's20011'

outfile = open(record+".csv", "a+")


# Read in the data from 0 to 10 seconds
data, info = rdsamp(record, 0)
pprint(info)

print "total time ", int(info['samp_count'])/int(info['samp_freq'])

ann = rdann(record, 'atr', 0)
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
	pprint(v[400:700])
        sys.exit(0)

	sigma = std(v)
	v_mean = mean(v)
	sdn = sigma/v_mean
	print "std=",sigma
	print "sdn=", sdn
	print "mean_RR=",v_mean

	print "correlation"
	cor = correlate(v,v,mode='full')
	cor = cor[len(cor)/2:]
	#print "cor_len=",len(cor)
	#pylab.figure(1)
	#pylab.plot(cor)
	#pylab.show()

	print "fft"
	F = abs(scipy.fftpack.fft(cor,len(cor)*10))
	F = F[:len(F)/2]
	freqs = scipy.fftpack.fftfreq(len(cor)*10, 1/info['samp_freq'])
	freqs = freqs[:len(freqs)/2]
	#freqs *= info['samp_freq']

	print "freqs cutting"
	f_min = freqs.min()
	f_max = freqs.min()
	i_min = 0
	i_max = 0
	i = 0
	for f in freqs:
		if f <= 0.003:
			f_min = f
			i_min = i
		if f <= 0.04:
			f_max = f
			i_max = i
		if f > 0.04:
			break
		i+=1

	#print "freqs " ,f_min, " - ", f_max
	#fmin = next(x for x in freqs if x >= 0.003)
	#fmax = next(x for x in freqs if x >= 0.04)

	F = F[i_min:i_max]
	#print F[:10]
	freqs = freqs[i_min:i_max] # right freqs only
	#print freqs[:10]

	#print len(v), len(F)
	#pprint(F[:10])
	#pprint(freqs[:10])

	F = scipy.log(F)
	#F /= max(F)
	#F -= min(F)
	freqs = scipy.log(freqs)
	#freqs /= max(freqs)
	#freqs -= min(freqs)


	A = scipy.vstack([freqs, ones(len(freqs))]).T
	#print "A=",A
	a, b = lstsq(A, F)[0]
	#print "a=",a," b=",b
	betta = -a
	print "betta=",betta
	#print "alpha=",atan(-a)/pi*180
	
	ostr=record+","+str(data[chunk_first,1])+","+str(data[chunk_last,1])+","+str(betta)+","+str(sigma)+","+str(sdn)+","+str(v_mean)+"\n"
	print ostr
	outfile.write(ostr)
	
	c +=w 


outfile.close()

def draw():

	#pprint(freqs[:10])

	fig = pylab.figure(2)
	#fig.add_subplot(211)
	#ax1.bar(1280, 128, color='red', edgecolor='black', hatch="o")
	sp1 = pylab.subplot(211)
	sp1.set_ylabel('Signal, mV')
	sp1.set_xlabel('Time, s')
	sp1.plot(t, v)
	sp1.plot(ann[c:c+w, 1], data[ann_x[c:c+w], 2], 'xr')
	#pylab.subplot(212)
	sp2 = fig.add_subplot(212)
	sp2.set_ylabel('ln(F)')
	sp2.set_xlabel('ln(freq)')
	#ax2.bar(1280, 128, color='red', edgecolor='black', hatch="o")
	#ax.plot(freqs,20*scipy.log10(F),',')
	#pylab.plot(freqs,20*scipy.log10(F),',')
	sp2.plot(freqs,F,'.')
	y = a*freqs+b
	pylab.plot(freqs, y, 'r', label='Fitted line')

	#ax.set_xscale('log')
	#ax.set_yscale('log')
	pylab.show()
