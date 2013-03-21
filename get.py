#!/usr/bin/env python
import urllib2

fpath = 'http://www.physionet.org/physiobank/database/ltstdb/'

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
        #print "Can't open %s" % (fpath+fname,)                                                                                                              
	pass

for a in range(2,4):
    for b in range(41,81):
        for c in range(1,5):
            for ext in ('dat','atr','hea'):
                fname = 's'+str(a)+'0'+str(b).zfill(2)+str(c)+'.'+ext
                get_file(fpath,fname)
