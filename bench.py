#!/usr/bin/env python

#from process import *
#import os
import timeit



t = timeit.timeit(stmt='process("chf01","ecg",1,1100)', setup='from process import *', number=10)
print t

