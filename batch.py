#!/usr/bin/env python

from process import *
import os
import sys

#files = ['16265', '16272', '16273', '16483', '16539', '16773', '16786', '17052', '17453', '18177', '18184', '19090', '19093', '19140', '19830']
#files = ['16265', '16272', '16273']
#files = ['s20011','s20021','s20031']
#files = ['s20011','s20021','s20031','s20041','s20051','s20061','s20071','s20081','s20091','s20101','s20111','s20121','s20131','s20141','s20151','s20161','s20171','s20181','s20191','s20201','s20211','s20221','s20231','s20241','s20251','s20261','s20271','s20281','s20291','s20301','s20311','s20321','s20331','s20341','s20351','s20361']

for f in os.listdir('.'):
    path = os.path.splitext(f)
    if path[1] == '.hea' and not os.path.exists("%s_std.png" % path[0]):
	print "*******************  Processing: ", path[0], "********************"
        try:
            if path[0][0] == 'c':
                process(path[0], 'ecg')
            else:
                process(path[0], 'atr')
        except:
            print "FAILED: ", sys.exc_info()

#for f in files:
#    process(f)
