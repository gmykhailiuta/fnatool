#!/usr/bin/env python
import pylab
from wfdbtools import rdsamp, rdann, plot_data
from pprint import pprint
from sys import exit
#import math
from scipy.signal import fftconvolve
import anfft

def write_csv(records,fname):
	outfile = open(fname, "w")
	for r in records:
		line = "%s,%s,%s,%s,%s,%s,%s\n" % (record,r['time_from'],r['time_to'],r['beta'],r['std'],r['sdn'],r['mean'])
		outfile.write(line)
	outfile.close()

def autocor(data):
	print "Autocorrelation"
	data_length = len(data)
	in2 = pylab.zeros(data_length * 2)

	in2[data_length/2:data_length/2+data_length] = data # This works for data_length being even

	# Do an array flipped convolution, which is a correlation.
	cor = fftconvolve(in2, data[::-1], mode='valid') 
	cor = cor[len(cor)/2:]
	return cor

def fft(x,samp_freq=250):
	print "FFT"
	arr = pylab.zeros(2**(round(pylab.log2(len(x)))+3))
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

def aprox(x,y):
	A = pylab.vstack([x, pylab.ones(len(x))]).T
	a, b = pylab.lstsq(A, y)[0]
	return a, b

#def beta(x,y):
#	a, b = aprox(x,y)
#	print "beta = %s" % (-_a, )
#	return -a

def fig(in1,in2=None,label="",show=False):	
	pylab.figure(figsize=(6, 6), facecolor='white')
	if in2 is None:
		pylab.plot(in1,label=label)
	else:
		pylab.plot(in1,in2,label=label)
	if show:
		pylab.show()
		exit(0)

def draw_sdn(rr, fname):
	x = [x['sdn'] for x in rr]
	y = [y['beta'] for y in rr]
	pylab.figure(figsize=(6, 6), facecolor='white')
	pylab.plot(x,y,'.r')
	pylab.axhline(y=1,color='b')
	pylab.ylim(0,2)
	pylab.xlabel(r"$\sigma/\bar{x}$")
	pylab.ylabel(r"$\beta$")
	pylab.savefig(fname,facecolor='w',edgecolor='k',transparent=True)
#	pylab.show()
#	exit(0)

def draw_std(rr, fname):
	x = [x['std'] for x in rr]
	y = [y['beta'] for y in rr]
	pylab.figure(figsize=(6, 6), facecolor='white')
	pylab.plot(x,y,'.r')
	pylab.axhline(y=1,color='b')
	pylab.ylim(0,2)
	pylab.xlabel(r"$\sigma$")
	pylab.ylabel(r"$\beta$")
	pylab.savefig(fname,facecolor='w',edgecolor='k',transparent=True)

def draw_flicker(f,F,y, fname):
	pylab.figure(figsize=(6, 6),facecolor='white')
	pylab.plot(f,F,'b')
	pylab.plot(f,y,'r')
	pylab.xlabel(r"$\log{f}$")
	pylab.ylabel(r"$\log{S}$")
	pylab.savefig(fname,facecolor='w',edgecolor='k',transparent=True)

def process(record, end=-1):
	print "Processing %s" % (record,)
	

	# Read in the data from 0 to 10 seconds
	# rdsamp(record, start=0, end=-1, interval=-1)
	data, info = rdsamp(record, 0, end)
	pprint(info)
#	pprint(data[:100])

	print "total time ", int(info['samp_count'])/int(info['samp_freq'])

# rdann(record, annotator, start=0, end=-1, types=[])
	ann = rdann(record, 'atr', 0, end)
# annotation time in samples from start
	ann_x = (ann[:, 0] - data[0, 0]).astype('int')


#ann_y = v[ann_y,2]
	#plot_data(data, info, ann)

	w = 1000 # window last RR

	print "RR count", len(ann)

	for channel in range(2):
		print "Channel %s" % (channel,)
		c = 0 # window first RR
		rr = []
		while c+w < len(ann):
			r = {}
		# get 10 RRs
			print "RR %s - %s / %s" % (c,c+w,len(ann))
			chunk_first = ann[c][0] # get sample number of first RR interval
			chunk_last = ann[c+w][0] # get sample number of last RR interval
		#time_interval = (chunk_last-chunk_first)/info['samp_freq']
			r['time_from'] = data[chunk_first,1]
			r['time_to'] = data[chunk_last,1]
			#print "interval ", r['time_from'], " - ", r['time_to'], " = ", (chunk_last-chunk_first)/info['samp_freq'], " s"

			t = data[chunk_first:chunk_last,1] # time array of target window
			v = data[chunk_first:chunk_last,channel+2] # data value of target window

	#		fig(t,v,show=True)
			r['std'] = pylab.std(v)
			r['mean'] = pylab.mean(v)
			r['sdn'] = r['std']/r['mean']

	#		vx = v/max(v) - r['mean']
	#		vx = v - min(v)
	#		vx = vx / max(v)
	#		v = vx

	#		fig(v)
	#		import pywt
	#		v = pylab.append(v,[0.])
	#		print len(v)
	#		output = pywt.swt(v, 'db6', level=pywt.swt_max_level(len(v)))
	#		output = pywt.thresholding.soft(output, r['std']*pylab.sqrt(2*pylab.log(len(v))))
	#		output = pywt.thresholding.zero(output)
	#		final = pywt.iswt(output, 'db6')
	#		cA = pylab.zeros(len(cA))
	#		fig(cA)
	#		fig(cD)
	#		cD = pywt.thresholding.soft(cD,0.4)
	#		vt = pywt.idwt(cA,cD,'db6')
	#		res += v + final[:len(v)]
	#		fig(res, show=True)

	#		from denoise import denoise
	#		J_dn=2
	#		nm_dn='sym4'
	#		mode_denoise='swt'
	#		dn_method='visushrink'
	#		dn_thresh='soft'
	#		out=denoise(v,mode_denoise,nm_dn,J_dn,dn_method,dn_thresh)
	#		fig(out, show=True)

			cor = autocor(v)

			Ff, ff = fft(cor, info['samp_freq'])
			Fc, fc = cut_freq(Ff,ff)

			Fc = pylab.log10(Fc)
			fc = pylab.log10(fc)
	#		fig(fc,Fc,show=False)

			_a,_b = aprox(fc,Fc)
			r['beta'] = -_a
			#draw_flicker(fc,Fc,_a*fc+_b,"%s_c%s_flicker_i%s.png" % (record, channel, int(c/w)))
			#exit(0)

			#pylab.show()
	#exit(0)


			pprint(r)
			rr.append(r)
			c +=w

		draw_std(rr, "%s_c%s_std.png" % (record, channel))
		draw_sdn(rr, "%s_c%s_sdn.png" % (record, channel))

		#write_csv(rr,record+"_c"+channel+".csv")    
		del rr
	del data, info, ann, ann_x	




if __name__ == '__main__':
	record = 's20161'
	process(record)
