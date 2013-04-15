import pylab as pl
import numpy as np
import sys

def bandPass(V):
 """
  Applies an Band Pass filter to V in the frequency domain.
  In:  V, list with signal
  Out: Frequency in lg scale
       Specturm power in lg scale
 """
 n = len(V)
 Ts = pl.mean(V)
 t = pl.arange(0,sum(V),float(Ts))
 pl.figure()
 pl.subplot(2,2,1)
 pl.title("Original signal")
 pl.plot(t, V)
 I = pl.fftshift(pl.fft(V)) # entering to frequency domain
 ffs = pl.fftshift(pl.fftfreq(len(V), Ts))
 # fftshift moves zero-frequency component 
 # to the center of the array
 P = pl.zeros(I.shape,dtype=complex)
 pl.subplot(2,2,2)
 pl.title("Original signal - Fourier Specturm")
 pl.plot(ffs,pl.log10(abs(I)))

 for i in range(len(V)):  # frequency cutting
  if 0.004 < abs(ffs[i]) and abs(ffs[i]) < 0.4:
   P[i] = I[i]

 pl.subplot(2,2,3)
 pl.title("Filtered signal")
 pl.plot(t, np.real(pl.ifft(pl.ifftshift(P))))
 pl.subplot(2,2,4)
 pl.title("Filtered signal - Fourier Specturm")
 pl.plot(ffs,pl.log10(abs(P)))
 pl.figure()
 pl.plot(pl.log10(ffs[range(n/2,n)]),pl.log10(abs(P))[range(n/2,n)])
 pl.show()

 return pl.log10(ffs[range(n/2,n)]),pl.log10(abs(P))[range(n/2,n)]

f = open("rr_16265.csv", "r")
V = []
for line in f:
 V.append(float(line.rstrip().split(",")[2])/1000)
f.close()


bandPass(V)