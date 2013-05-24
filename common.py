#!/usr/bin/env python
import csv
import pylab as pl

SIGNALS = [
    {'diagnosis': 'Normal',
      'annotator': 'atr',
      'plot_color': 'g',
      'records': '16265 16273 16483 16773 16795 17453 18184 19090 19140 16272 16420 16539 16786 17052 18177 19088 19093 19830'},
    {'diagnosis': 'CHF',
      'annotator': 'ecg',
      'plot_color': 'r',
      'records': 'chf01 chf03 chf05 chf07 chf09 chf11 chf13 chf15 chf02 chf04 chf06 chf08 chf10 chf12 chf14'},
    {'diagnosis': 'Hypertension',
      'annotator': 'atr',
      'plot_color': 'b',
      'records': 's20031 s20121 s20221 s20471 s20551 s20651 s30691 s30751 s30791 s20051 s20131 s20341 s20481 s20561 s30661 s30741 s30752 s30801 s20101 s20171 s20411 s20501 s20581 s30681 s30742 s30761'
    }]

def read_signal(file_name=None):
    beta, std, cov, mean = [], [], [], []
    result = {'record':False}
    try:
        with open(file_name, 'rb') as f:
            reader = csv.reader(f)
            for row in reader:
                if result['record']:
                    beta.append(float(row[2]))
                    std.append(float(row[3]))
                    cov.append(float(row[4]))
                    mean.append(float(row[5]))
                else:
                    result['record'] = row[0]
                    result['gender'] = row[1]
                    result['age'] = row[2]
                    result['diagnosis'] = row[3]
    except IOError:
        print 'Could not read %s' % (file_name,)
        return None
    result['beta'] = beta
    result['std'] = std
    result['cov'] = cov
    result['mean'] = mean
    #pprint(result)
    return result


def filter2d(x, y, filtration_algo='2sigma', valid_delta_ratio=.2):
    xnew = pl.array(x, dtype='float32')
    ynew = pl.array(y, dtype='float32')
    if filtration_algo == 'rate20':
        for i in range(0, len(ynew)-1):
            if (abs(ynew[i+1] / ynew[i] - 1) > valid_delta_ratio): # current rr differs more then 20% of previous one
                ynew[i] = 0
                xnew[i] = 0
        ynew = pl.ma.masked_equal(ynew,0)
        ynew = pl.ma.compressed(ynew)
        xnew = pl.ma.masked_equal(xnew,0)
        xnew = pl.ma.compressed(xnew)
    elif filtration_algo == '2sigma':
        mean = pl.mean(ynew)
        std = pl.std(ynew)
      
        for i in range(0, len(ynew)):
            if pl.logical_or(ynew[i] < mean - 2*std, mean + 2*std < ynew[i]):
                ynew[i] = 0
                xnew[i] = 0

        ynew = pl.ma.masked_equal(ynew,0)
        ynew = pl.ma.compressed(ynew)
        xnew = pl.ma.masked_equal(xnew,0)
        xnew = pl.ma.compressed(xnew)

    elif filtration_algo == 'rate20_2sigma':
        for i in range(0, len(ynew)-1):
            if (abs(ynew[i+1] / ynew[i] - 1) > valid_delta_ratio): # current rr differs more then 20% of previous one
                ynew[i] = 0
                xnew[i] = 0
        ynew = pl.ma.masked_equal(ynew,0)
        ynew = pl.ma.compressed(ynew)
        xnew = pl.ma.masked_equal(xnew,0)
        xnew = pl.ma.compressed(xnew)

        mean = pl.mean(ynew)
        std = pl.std(ynew)
      
        for rr in range(0, len(ynew)):
            if pl.logical_or(ynew[i] < mean - 2*std, mean + 2*std < ynew[i]):
                ynew[i] = 0
                xnew[i] = 0
    else:
       warn("Donno anything about such filtration algorithm")

    print "NB: Deleted %0.2f%% intervals" % (float(len(y)-len(ynew))/len(y)*100,)
    return xnew, ynew



def signal_to_csv(record,time,hrv,info):
    outfile = open("%s.csv" % (record,), "w")
    for i in range(len(hrv)):
        t_full = info['base_time'] + dt.timedelta(seconds=time[i])
        line = "%s %s %s\n" % (i+1,t_full.strftime("%H:%M:%S.%f"),hrv[i])
        outfile.write(line)
    outfile.close()

