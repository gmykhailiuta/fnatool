#!/usr/bin/env python
import urllib2
import os

fpath = 'http://www.physionet.org/physiobank/database/chfdb/'

def get_file(fpath,fname):
    try:
        u = urllib2.urlopen(fpath+fname)
	if u.code == 200:
            localFile = open(fname, 'w')
            localFile.write(u.read())
            localFile.close()
            print 'Got %s' % (fpath+fname, )
        else:
            print "File %s not found (return code %s)" % (fpath+fname, u.code)
    except urllib2.HTTPError:
        print "Can't open %s" % (fpath+fname,)                                                                                                              
	pass

def get_file2(fpath,fname):
    os.system("wget -c "+fpath+fname)

u = urllib2.urlopen(fpath+"RECORDS")
records = u.read().split('\n')
u.close()

print records

for r in records:
    get_file2(fpath,r+".hea")
    get_file2(fpath,r+".dat")
    get_file2(fpath,r+".ecg")
