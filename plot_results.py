#!/usr/bin/env python
import matplotlib as mp
#mp.use('GtkAgg')
import pylab as pl
from pprint import pprint
from scipy.cluster.vq import kmeans2
from matplotlib.dates import HourLocator, DateFormatter
import datetime as dt
import scipy.interpolate as ip
from warnings import warn
from convex_hull import convex_hull
from theil_sen import theil_sen
import itertools
import common
global config
A4_WIDTH = 11.7
A4_HEIGHT = 8.3

def draw(in1,in2=None,label="",show=False):  
    pl.figure(facecolor='white')
    if in2 is None:
        pl.plot(in1,label=label)
    else:
        pl.plot(in1,in2,label=label)
    if show:
        pl.show()
    exit(0)

def plot_hrv(hrv, hrv_filt, hrv_interp, record, preview=False):
    """
    Plot HRV and it's approximation.
    In:
        hrv : list, [time, hrv], ndarray,ndarray
            Time & hrv vectors
        hrv_interp : list, [time_interp, hrv_interp], ndarray,ndarray
            Interpolated time & hrv vectors
    """
    fig = pl.figure("signals", figsize=(2*A4_WIDTH, 2*A4_HEIGHT), facecolor='white')
    pl.suptitle(record)
    sp_orig = pl.subplot(311)
    pl.title("Original signal")
    pl.ylabel(r"$RR_i$ (ms)") 
    pl.plot(hrv[0], hrv[1], color='b', marker='.', linestyle='--', alpha=.5)
    sp_filt = pl.subplot(312,sharex=sp_orig,sharey=sp_orig)
    pl.title("Filtered signal")
    pl.ylabel("$RR_i$ (ms)")
    pl.plot(hrv_filt[0], hrv_filt[1], color='b', marker='.', linestyle='--', alpha=.5)
    sp_filt = pl.subplot(313,sharex=sp_orig,sharey=sp_orig)
    pl.title("Interpolated signal")
    pl.ylabel(r"$RR_j$ (ms)") 
    pl.xlabel('Time (s)')
    pl.plot(hrv_interp[0], hrv_interp[1], color='k')
    pl.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    if preview:
        pl.show()
    pl.close()

def plot_time(rrs, info, results, window, preview=False):
    rr_times = pl.array([common.elapsed_to_abs_time(t[0],info['base_time']) for t in rrs])
    beg_times = pl.array([x['time_beg'] for x in results])
    beta = pl.array([y['beta'] for y in results],dtype='float32')
    std = pl.array([y['std'] for y in results],dtype='float32')

    fig, (sp1,sp2,sp3) = pl.subplots(3,1,sharey=False,sharex=True,figsize=(12, 12),facecolor='white',num="beta_time")
    fig.suptitle('Parameters in time', fontsize=20)    
    sp1.plot(rr_times, rrs[:,1], color='k')
    sp1.vlines(beg_times,rrs[:,1].min(),rrs[:,1].min(),linestyle='solid',color=[0.7,0.7,0.7])
    sp1.set_ylabel(r"$RR_i$ (ms)", fontsize=16)
    sp1.grid(True)
    sp2.bar(beg_times, beta, width=(1./24./60/60*window),fill=False)
    sp2.axhline(y=1,color='r',linestyle='--')
    sp2.xaxis.set_major_locator(HourLocator())
    sp2.xaxis.set_major_formatter(DateFormatter("%H"))
    sp2.set_ylabel(r"$\beta$", fontsize=16)
    #sp2.xaxis.set_minor_locator(MinuteLocator(interval=15))
    #sp2.autoscale_view()
    #sp2.xaxis.grid(True, 'major')
    #sp2.xaxis.grid(True, 'minor')
    sp2.grid(True)
    #fig.autofmt_xdate()
    sp3.bar(beg_times, std, width=(1./24./60/60*window), fill=False)
    sp3.axhline(y=70,color='r',linestyle='--')
    sp3.set_ylabel(r"$\sigma$ (ms)", fontsize=16)
    sp3.set_xlabel("Time (hours)", fontsize=16)
    sp3.grid(True)
    fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    fig.savefig("time_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


# def plot_beta_cv(results, preview=False):
#     x = [x['cov'] for x in results]
#     y = [y['beta'] for y in results]
#     pl.figure(figsize=(6, 6), facecolor='white')
#     pl.title(r'%(record)s: $\sigma/\bar{x}$' % results[0])
#     pl.ylim(0,2)
#     pl.xlim(0,50)
#     pl.xlabel(r"$\sigma/\bar{x}$")
#     pl.ylabel(r"$\beta$")
#     pl.grid()
#     pl.axhline(y=1,color='k',linestyle='--')
#     pl.plot(x,y,'.r')
#     pl.savefig("cv_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
#     if preview:
#         pl.show()
#     pl.close()


# def plot_beta_std(results, preview=False):
#     x = pl.array([x['std'] for x in results])
#     y = pl.array([y['beta'] for y in results])
#     pl.figure(figsize=(6, 6), facecolor='white')
#     pl.plot(x,y,'.r')
#     #rx = pl.linspace(0,140,140)
#     #ra, rb = theil_sen(x,y)
#     #ra, rb = common.approximate(x,y)
#     #pl.plot(rx,ra*rx+rb,alpha=.9,lw=2)
#     pl.title(r'%(record)s: $\beta/\sigma$' % results[0])
#     pl.ylim(0,2)
#     pl.xlim(0,140)
#     pl.xlabel(r"$\sigma$, ms")
#     pl.ylabel(r"$\beta$")
#     pl.grid()
#     pl.axhline(y=1,color='k',linestyle='--')
#     pl.axvline(x=70,color='k',linestyle='--')
#     pl.savefig("std_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
#     if preview:
#         pl.show()
#     pl.close()

def plot_homeostasis_record(results, preview=False):
    std = [x['std'] for x in results]
    cov = [x['cov'] for x in results]
    beta = [y['beta'] for y in results]
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num="stochastic_homeostasis_record")
    fig.suptitle('Stochastic homeostasis for record %(record)s' % results[0], fontsize=20) 
    sp1.scatter(std,beta,color='b')
    sp1.scatter(pl.mean(std),pl.mean(beta),s=200,color='b',marker='o',alpha=.5)
    sp1.axhline(y=1,color='k',linestyle='--')
    sp1.axvline(x=70,color='k',linestyle='--')
    sp1.set_ylim(0,2)
    sp1.set_xlim(0,140)
    sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
    sp1.set_ylabel(r"$\beta$", fontsize=18)
    sp1.grid()
    sp2.scatter(cov,beta,color='b')
    sp2.scatter(pl.mean(cov),pl.mean(beta),s=200,color='b',marker='o',alpha=.5)
    sp2.axhline(y=1,color='k',linestyle='--')
    sp2.set_xlim(0,20)
    sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
    sp2.grid()
    fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    fig.savefig("homeostasis_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_beta(freq, fft, aprox_y, result, preview=False):
    pl.figure(figsize=(6, 6),facecolor='white')
    pl.plot(freq,fft,'b')
    pl.plot(freq,aprox_y,'r')
    pl.title('%(record)s: %(time_beg)s - %(time_end)s' % result)
    pl.xlim(min(freq),max(freq))
    pl.xlabel(r"$\lg f$")
    pl.ylabel(r"$\lg S$")
    pl.grid()
    pl.savefig("beta_%(record)s_%(frag)s.png" % result,facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_homeostasis(preview=False):
    global config
    fig_sum, (sp1_sum, sp2_sum) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num="stochastic_homeostasis")
    fig_sum.suptitle('Stochastic homeostasis by diagnosis', fontsize=20) 
    sp1_sum.axhline(y=1,color='k',linestyle='--')
    sp1_sum.axvline(x=70,color='k',linestyle='--')
    sp1_sum.set_ylim(0,2)
    sp1_sum.set_xlim(0,140)
    sp1_sum.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
    sp1_sum.set_ylabel(r"$\beta$", fontsize=18)
    sp1_sum.grid()
    sp2_sum.axhline(y=1,color='k',linestyle='--')
    sp2_sum.set_xlim(0,20)
    sp2_sum.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
    sp2_sum.grid()

    #legend = []
    diag_data = [[],[],[]]
    for s in range(len(config['SIGNALS'])):
        db = config['SIGNALS'][s]
        record_data = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_data[0].extend(signal['beta'])
                record_data[1].extend(signal['std'])
                record_data[2].extend(signal['cov'])

        #diag_beta.append(record_data[0])
        for i in range(3):
            diag_data[i].append(pl.array(record_data[i]))

        fig_diag, (sp1_diag, sp2_diag) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num=("stochastic_homeostasis_%s" % db['diagnosis']))
        fig_diag.suptitle('Stochastic homeostasis of diagnosis %s' % db['diagnosis'], fontsize=20)
        
        sp1_diag.scatter(record_data[1],record_data[0],label=db['diagnosis'],color=db['plot_color'])
        sp1_diag.scatter(pl.mean(record_data[1]),pl.mean(record_data[0]),s=200,color=db['plot_color'],marker='o',alpha=.5)
        sp1_diag.axhline(y=1,color='k',linestyle='--')
        sp1_diag.axvline(x=70,color='k',linestyle='--')
        sp1_diag.set_ylim(0,2)
        sp1_diag.set_xlim(0,140)
        sp1_diag.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
        sp1_diag.set_ylabel(r"$\beta$", fontsize=18)
        sp1_diag.grid()
        
        sp2_diag.scatter(record_data[2],record_data[0],label=db['diagnosis'],color=db['plot_color'])
        sp2_diag.scatter(pl.mean(record_data[1]),pl.mean(record_data[0]),s=200,color=db['plot_color'],marker='o',alpha=.5)
        sp2_diag.axhline(y=1,color='k',linestyle='--')
        sp2_diag.set_xlim(0,20)
        sp2_diag.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        sp2_diag.grid()

        fig_diag.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
        fig_diag.savefig("stochastic_homeostasis_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

        x, y = diag_data[1][s],diag_data[0][s]
        sp1_sum.scatter(x,y,label=db['diagnosis'],color=db['plot_color'],marker='.')
        x, y = diag_data[2][s],diag_data[0][s]
        sp2_sum.scatter(x,y,label=db['diagnosis'],color=db['plot_color'],marker='.')

    fig_sum.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    sp2_sum.legend(loc='best')
    fig_sum.savefig("stochastic_homeostasis_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    
    if preview:
        pl.show()
    pl.close()


def boxplot_diagnosis(preview=False):
    global config
    legend = []
    diagnosis_stat = [[],[],[],[]]
    #betas, stds, covs, means = [], [], [], []
    for db in config['SIGNALS']:
        record_stat = [[],[],[],[]]
        #record_stat_filt = [[],[],[],[]]
        records = db['records'].split()
        for record in records:
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_stat[0].append(signal['beta'] * common.filter1d(signal['beta'], algos=['5per95']))
                record_stat[1].append(signal['std'] * common.filter1d(signal['std'], algos=['5per95']))
                record_stat[2].append(signal['cov'] * common.filter1d(signal['cov'], algos=['5per95']))
                record_stat[3].append(signal['mean'] * common.filter1d(signal['mean'], algos=['5per95']))
        fig = pl.figure("boxplot_%s" % db['diagnosis'], figsize=(23, 16), facecolor='white')
        fig.suptitle('Statistics for diagnisis %s' % db['diagnosis'], fontsize=20)
        xticks = range(1,len(records)+1)
        sp_beta = pl.subplot(221)
        pl.title(r"$\beta$", fontsize=18)
        #pprint(record_stat[0])
        pl.axhline(y=1,color='g',linestyle='--')
        sp_beta.boxplot(record_stat[0])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[0]], marker='x',color='k')
        pl.scatter(xticks, [pl.std(x) for x in record_stat[0]], marker='o',color='k')
        #pl.ylim(0,pl.percentile())
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(-.2,2)
        sp_beta.yaxis.grid()
        sp_std = pl.subplot(222,sharex=sp_beta)
        pl.title(r"$\sigma$", fontsize=18)
        pl.axhline(y=70,color='g',linestyle='--')
        sp_std.boxplot(record_stat[1])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[1]], marker='x',color='k')
        pl.scatter(xticks, [pl.std(x) for x in record_stat[1]], marker='o',color='k')
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(0,400)
        sp_std.yaxis.grid()
        sp_cov = pl.subplot(223,sharex=sp_beta)
        pl.title(r"$\sigma/\bar{x}$", fontsize=18)
        sp_cov.boxplot(record_stat[2])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[2]], marker='x',color='k')
        pl.scatter(xticks, [pl.std(x) for x in record_stat[2]], marker='o',color='k')
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(0,50)
        sp_cov.yaxis.grid()
        sp_mean = pl.subplot(224,sharex=sp_beta)
        pl.title(r"$\bar{x}$", fontsize=18)
        sp_mean.boxplot(record_stat[3])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[3]], marker='x',color='k')
        pl.scatter(xticks, [pl.std(x) for x in record_stat[3]], marker='o',color='k')
        pl.ylim(400,1500)
        sp_mean.yaxis.grid()
        pl.xticks(xticks,records, rotation='vertical')
        pl.subplots_adjust(left=0.06, right=0.95, top=0.91, bottom=0.09, wspace=0.15, hspace=0.3)
        pl.savefig("boxplot_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()

        for i in range(4):
            diagnosis_stat[i].append(list(itertools.chain.from_iterable(record_stat[i])))
        legend.append(db['diagnosis'][:3])

    for s in range(len(diagnosis_stat)):
        for r in range(len(diagnosis_stat[s])):
            diagnosis_stat[s][r] *= common.filter1d(diagnosis_stat[s][r], algos=['5per95'])
    fig = pl.figure("boxplot_diagnosis_summary", figsize=(23, 16), facecolor='white')
    fig.suptitle('Statistics by diagnoses', fontsize=20)   
    xticks = range(1,len(config['SIGNALS'])+1)
    sp_beta = pl.subplot(221)
    pl.title(r"$\beta$", fontsize=18)
    pl.axhline(y=1,color='g',linestyle='--')
    sp_beta.boxplot(diagnosis_stat[0])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[0]], marker='x',color='k')
    pl.scatter(xticks, [pl.std(x) for x in diagnosis_stat[0]], marker='o',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(-.2,2)
    sp_beta.yaxis.grid()
    pl.ylabel(r"$\beta$", fontsize=18)
    sp_std = pl.subplot(222,sharex=sp_beta)
    pl.title(r"$\sigma$", fontsize=18)
    pl.axhline(y=70,color='g',linestyle='--')
    sp_std.boxplot(diagnosis_stat[1])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[1]], marker='x',color='k')
    pl.scatter(xticks, [pl.std(x) for x in diagnosis_stat[1]], marker='o',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(0,400)
    sp_std.yaxis.grid()
    sp_cov = pl.subplot(223,sharex=sp_beta)
    pl.title(r"$\sigma/\bar{x}$", fontsize=18)
    sp_cov.boxplot(diagnosis_stat[2])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[2]], marker='x',color='k')
    pl.scatter(xticks, [pl.std(x) for x in diagnosis_stat[2]], marker='o',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(0,50)
    sp_cov.yaxis.grid()
    sp_mean = pl.subplot(224,sharex=sp_beta)
    pl.title(r"$\bar{x}$", fontsize=18)
    sp_mean.boxplot(diagnosis_stat[3])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[3]], marker='x',color='k')
    pl.scatter(xticks, [pl.std(x) for x in diagnosis_stat[3]], marker='o',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(400,1500)
    sp_mean.yaxis.grid()
    pl.subplots_adjust(left=0.06, right=0.95, top=0.91, bottom=0.09, wspace=0.15, hspace=0.3)
    pl.savefig("boxplot_diagnosis_summary.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_clusters(preview=False):
    global config
    pl.figure("clusterization", figsize=(22, 11), facecolor='white')
    pl.title('Clusterization by desease')
    pl.suptitle('Clusterization', fontsize=20) 
    
    sp_beta_std = pl.subplot(121)
    pl.axhline(y=1,color='k')
    pl.axvline(x=70,color='k')
    pl.xlabel(r"$\sigma\ (ms)$", fontsize=18)
    pl.ylabel(r"$\beta$", fontsize=18)
    pl.grid()

    sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
    pl.axhline(y=1,color='k')
    pl.xlabel(r"$\sigma/\bar{x}$", fontsize=18)
    pl.grid()

    stat = [[],[],[]]
    #summary = []
    for db in config['SIGNALS']:
        record_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                stat[0].extend(signal['beta'])
                stat[1].extend(signal['std'])
                stat[2].extend(signal['cov'])


    std, beta = common.filter2d(stat[1], stat[0], axes=['x','y'], algos=['3per97'])
    res, idx = kmeans2(pl.array(zip(std, beta)),len(config['SIGNALS']))
    colors = ([config['SIGNALS'][i]['plot_color'] for i in idx])
    sp_beta_std.scatter(std, beta, c=colors)
    sp_beta_std.scatter(res[:,0],res[:,1], marker='o', s = 1500, linewidths=2, c='w', alpha=0.5)
    sp_beta_std.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2, c='k')

    cov, beta = common.filter2d(stat[2], stat[0], axes=['x','y'], algos=['3per97'])
    res, idx = kmeans2(pl.array(zip(cov, beta)),len(config['SIGNALS']))
    colors = ([config['SIGNALS'][i]['plot_color'] for i in idx])
    sp_beta_cov.scatter(cov, beta, c=colors)    
    sp_beta_cov.scatter(res[:,0],res[:,1], marker='o', s = 1500, linewidths=2, c='w', alpha=0.5)
    sp_beta_cov.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2, c='k')

    pl.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1)
    pl.savefig("clusterization.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_homeostasis_interp(preview=False):
    global config
    for db in config['SIGNALS']:
        record_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_stat[0].extend(signal['beta'])
                record_stat[1].extend(signal['std'])
                record_stat[2].extend(signal['cov'])
        
        fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(14, 6),facecolor='white',num=("stochastic_homeostasis_summary_interp"))
        fig.suptitle('Stochastic homeostasis: interpolation for diagnosis %s' % (db['diagnosis'],), fontsize=20)
        sp1.axhline(y=1,color='k',linestyle='--')
        sp1.axvline(x=70,color='k',linestyle='--')
        sp1.set_ylim(0,2)
        sp1.set_xlim(0,140)
        sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
        sp1.set_ylabel(r"$\beta$", fontsize=18)
        sp1.grid()
        sp2.axhline(y=1,color='k',linestyle='--')
        sp2.set_xlim(0,50)  
        sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        sp2.grid()

        #x, y = common.filter2d(record_stat[1], record_stat[0], axes=['x','y'], algos=['3per97'])
        x, y = record_stat[1], record_stat[0]
        s1 = sp1.hexbin(x, y,extent=[0,140,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
        s1.set_clim()
        fig.colorbar(s1,ax=sp1)
        #x, y = common.filter2d(record_stat[2], record_stat[0], axes=['x','y'], algos=['3per97'])
        x, y = record_stat[2], record_stat[0]
        s1 = sp2.hexbin(x, y,extent=[0,50,0,2],gridsize=20,mincnt=2,label=db['diagnosis'])
        s1.set_clim()
        fig.colorbar(s1,ax=sp2)
        #cb.set_label('log10(N)')

        fig.tight_layout()
        fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.1, wspace=0.05, hspace=.05)
        fig.savefig("stochastic_homeostasis_diagnosis_interp_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()
        pl.close()


def plot_homeostasis_median(preview=False):
    global config
    
    fig_sum, (sp1_sum, sp2_sum) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num=("stochastic_homeostasis_medians_summary"))
    fig_sum.suptitle('Stochastic homeostasis: medians by diagnosis', fontsize=20) 
    sp1_sum.axhline(y=1,color='k',linestyle='--')
    sp1_sum.axvline(x=70,color='k',linestyle='--')
    sp1_sum.set_ylim(0,2)
    sp1_sum.set_xlim(0,140)
    sp1_sum.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
    sp1_sum.set_ylabel(r"$\beta$", fontsize=18)
    sp1_sum.grid()
    sp2_sum.axhline(y=1,color='k',linestyle='--')
    sp2_sum.set_xlim(0,10)
    sp2_sum.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
    sp2_sum.grid()

    #legend = []
    diag_data = [[],[],[]]
    diag_median = [[],[],[]]
    for db in config['SIGNALS']:
        record_data = [[],[],[]]
        record_median = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_data[0].append(signal['beta'])
                record_data[1].append(signal['std'])
                record_data[2].append(signal['cov'])
                record_median[0].append(pl.median(signal['beta']))
                record_median[1].append(pl.median(signal['std']))
                record_median[2].append(pl.median(signal['cov']))

        fig_diag, (sp1_diag, sp2_diag) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),num=("stochastic_homeostasis_median_%s" % db['diagnosis']))
        fig_diag.suptitle('Stochastic homeostasis: medians by record of diagnosis %s' % db['diagnosis'], fontsize=20)
        
        sp1_diag.scatter(record_median[1],record_median[0],label=db['diagnosis'],color=db['plot_color'])
        sp1_diag.scatter(pl.mean(record_median[1]),pl.mean(record_median[0]),s=200,color=db['plot_color'],marker='o',alpha=.5)
        sp1_diag.axhline(y=1,color='k',linestyle='--')
        sp1_diag.axvline(x=70,color='k',linestyle='--')
        sp1_diag.set_ylim(0,2)
        sp1_diag.set_xlim(0,140)
        sp1_diag.set_xlabel(r"$\sigma\ (ms)$", fontsize=18)
        sp1_diag.set_ylabel(r"$\beta$", fontsize=18)
        sp1_diag.grid()
        
        sp2_diag.scatter(record_median[2],record_median[0],label=db['diagnosis'],color=db['plot_color'])
        sp2_diag.scatter(pl.mean(record_median[2]),pl.mean(record_median[0]),s=200,color=db['plot_color'],marker='o',alpha=.5)
        sp2_diag.axhline(y=1,color='k',linestyle='--')
        sp2_diag.set_xlim(0,20)
        sp2_diag.set_xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        sp2_diag.grid()

        fig_diag.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
        fig_diag.savefig("stochastic_homeostasis_median_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

        for i in range(3):
            diag_data_1 = pl.array(list(itertools.chain.from_iterable(record_data[i])))
            diag_data[i].append(diag_data_1)
            diag_median[i] = pl.median(diag_data_1)
        
        sp1_sum.scatter(diag_median[1],diag_median[0],s=200,label=db['diagnosis'],color=db['plot_color'],marker='o')
        sp2_sum.scatter(diag_median[2],diag_median[0],s=200,label=db['diagnosis'],color=db['plot_color'],marker='o')

    fig_sum.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    sp2_sum.legend(loc='best')
    fig_sum.savefig("stochastic_homeostasis_medians_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    
    pl.close()


def plot_homeostasis_contour(preview=False):
    global config
    #legend = []
    # diagnosis_stat = [[],[],[]]
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(24, 12),facecolor='white',num="stochastic_homeostasis_contours_summary")
    fig.suptitle('Stochastic homeostasis: contour by diagnosis', fontsize=36)
    for i in range(len(config['SIGNALS'])):
        db = config['SIGNALS'][i]
        records_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                records_stat[0].extend(signal['beta'])
                records_stat[1].extend(signal['std'])
                records_stat[2].extend(signal['cov'])

        std, beta = common.filter2d(records_stat[1], records_stat[0], axes=['x','y'], algos=['1per99'])
        ch1 = convex_hull(pl.transpose(zip(std,beta)),False)
        cov, beta = common.filter2d(records_stat[2], records_stat[0], axes=['x','y'], algos=['1per99'])
        ch2 = convex_hull(pl.transpose(zip(cov,beta)),False)
        sp1.fill(ch1[:,0],ch1[:,1],label=db['diagnosis'],color=db['plot_color'],alpha=.2,lw=2)
        sp2.fill(ch2[:,0],ch2[:,1],label=db['diagnosis'],color=db['plot_color'],alpha=.2,lw=2)


    sp1.axhline(y=1,color='k',linestyle='--')
    sp1.axvline(x=70,color='k',linestyle='--')
    sp1.set_ylim(0,2.5)
    sp1.set_xlim(0,140)
    sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=32)
    sp1.set_ylabel(r"$\beta$", fontsize=32)
    sp1.grid()
        
    sp2.axhline(y=1,color='k',linestyle='--')
    sp2.set_xlim(0,100)
    sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=32)
    sp2.grid()
    sp2.legend(loc='best',prop={'size':24})

    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.2)
    fig.savefig("stochastic_homeostasis_contours_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        fig.show()
    pl.close()


def plot_homeostasis_contour_median(preview=False):
    global config
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(24, 12),facecolor='white',num="stochastic_homeostasis_contours_median_summary")
    fig.suptitle('Stochastic homeostasis by diagnosis (contours_medians)', fontsize=36)
    for db in config['SIGNALS']:
        records_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                records_stat[0].append(pl.median(signal['beta']))
                records_stat[1].append(pl.median(signal['std']))
                records_stat[2].append(pl.median(signal['cov']))
        for i in range(3):
            records_stat[i] = pl.array(records_stat[i])

        ch1 = convex_hull(pl.transpose(zip(records_stat[1], records_stat[0])),False)
        ch2 = convex_hull(pl.transpose(zip(records_stat[2], records_stat[0])),False)
        sp1.fill(ch1[:,0],ch1[:,1],label=db['diagnosis'],color=db['plot_color'],alpha=.2,lw=1)
        sp2.fill(ch2[:,0],ch2[:,1],label=db['diagnosis'],color=db['plot_color'],alpha=.2,lw=1)
        sp1.fill(ch1[:,0],ch1[:,1],color=db['plot_color'],alpha=.5,lw=1,fill=False,linestyle='dashed')
        sp2.fill(ch2[:,0],ch2[:,1],color=db['plot_color'],alpha=.5,lw=1,fill=False,linestyle='dashed')

        x = pl.linspace(0,140,140)
        a,b = theil_sen(records_stat[1],records_stat[0])
        #a,b = common.approximate(records_stat[1],records_stat[0])
        sp1.plot(x,a*x+b,color=db['plot_color'],alpha=.9,lw=2)
        #a,b = common.approximate(records_stat[2],records_stat[0])
        a,b = theil_sen(records_stat[2],records_stat[0])
        sp2.plot(x,a*x+b,color=db['plot_color'],alpha=.9,lw=2)


    sp1.axhline(y=1,color='k',linestyle='--')
    sp1.axvline(x=70,color='k',linestyle='--')
    sp1.set_ylim(0,2.2)
    sp1.set_xlim(0,140)
    sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=32)
    sp1.set_ylabel(r"$\beta$", fontsize=32)
    sp1.grid()
        
    sp2.axhline(y=1,color='k',linestyle='--')
    sp2.set_xlim(0,60)
    sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=32)
    sp2.grid()
    sp2.legend(loc='best',prop={'size':24})

    pl.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.2)
    fig.savefig("stochastic_homeostasis_contours_median_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()

def plot_metric_median_age(preview=False):
    global config
    
    fig_sum, ((sp1_sum, sp2_sum), (sp3_sum, sp4_sum)) = pl.subplots(2,2,sharey=False,sharex=False,figsize=(9, 9),facecolor='white',num=("metric_median_age"))
    fig_sum.suptitle('Metric by age : medians by diagnosis', fontsize=20) 

    #legend = []
    diag_data = [[],[],[]]
    for db in config['SIGNALS']:
        record_data = []
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal and signal['age'] and signal['age']:
                record_data.append([pl.int8(signal['age']),pl.median(signal['beta']),pl.median(signal['std']),pl.median(signal['cov']),pl.median(signal['mean'])])

        fig, ((sp1, sp2), (sp3, sp4)) = pl.subplots(2,2,sharey=False,sharex=False,figsize=(9, 9),facecolor='white',num=("metric_median_age_%s" % db['diagnosis']))
        fig.suptitle('Metric medians by age of diagnosis %s' % db['diagnosis'], fontsize=16)
        record_data = pl.array(record_data)

        def _draw_sp(sp,sp_sum,x,y,ylabel):
            a,b = theil_sen(x, y)
            _x = pl.arange(90)
            _y = a * _x + b
            #pprint(pl.square(pl.corrcoef(x,y)))
            sp.scatter(x,y,color=db['plot_color'],alpha=.5)
            sp.plot(_x,_y,color=db['plot_color'],label=(r"$y=%0.3f%+0.3f*x$"%(b,a)  + '\n' + (r"$R^2=%0.3f$"%(pl.square(pl.corrcoef(x,y)[0,1]),))),lw=2)
            sp.set_xlim(0,90)
            sp.set_xlabel("Age (years)", fontsize=14)
            sp.set_ylabel(ylabel, fontsize=14)
            sp.grid()
            sp.legend(loc='best',prop={'size':12})
            sp_sum.plot(_x,_y,color=db['plot_color'],label=db['diagnosis'])
            sp_sum.set_xlabel("Age (years)", fontsize=14)
            sp_sum.set_ylabel(ylabel, fontsize=14)
            sp_sum.grid()

        _draw_sp(sp1,sp1_sum,record_data[:,0],record_data[:,1],r"$\beta$")
        _draw_sp(sp2,sp2_sum,record_data[:,0],record_data[:,2],r"$\sigma$")
        _draw_sp(sp3,sp3_sum,record_data[:,0],record_data[:,3],r"$\sigma/\bar{x}$")
        _draw_sp(sp4,sp4_sum,record_data[:,0],record_data[:,4],r"$\bar{x}$")
        fig.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.06, wspace=0.25, hspace=.2)
        fig.savefig("metric_median_age_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

    fig_sum.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.06, wspace=0.25, hspace=.2)
    sp1_sum.legend(loc='best',prop={'size':12})
    sp2_sum.legend(loc='best',prop={'size':12})
    sp3_sum.legend(loc='best',prop={'size':12})
    sp4_sum.legend(loc='best',prop={'size':12})
    fig_sum.savefig("metric_median_age.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    
    pl.close()


def plot_metric_age(preview=False):
    global config
    
    fig_sum, ((sp1_sum, sp2_sum), (sp3_sum, sp4_sum)) = pl.subplots(2,2,sharey=False,sharex=False,figsize=(9, 9),facecolor='white',num=("metric_age"))
    fig_sum.suptitle('Metric by age : medians by diagnosis', fontsize=20) 

    #legend = []
    diag_data = [[],[],[]]
    for db in config['SIGNALS']:
        record_data = []
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal and signal['age'] and signal['age']:
                for i in range(len(signal['beta'])):
                    record_data.append([pl.int8(signal['age']),signal['beta'][i],signal['std'][i],signal['cov'][i],signal['mean'][i]])

        fig, ((sp1, sp2), (sp3, sp4)) = pl.subplots(2,2,sharey=False,sharex=False,figsize=(9, 9),facecolor='white',num=("metric_age_%s" % db['diagnosis']))
        fig.suptitle('Metric by age of diagnosis %s' % db['diagnosis'], fontsize=16)
        record_data = pl.array(record_data)

        def _draw_sp(sp,sp_sum,x,y,ylabel):
            a,b = theil_sen(x, y)
            _x = pl.arange(90)
            _y = a * _x + b
            sp.scatter(x,y,color=db['plot_color'],alpha=.5)
            sp.plot(_x,_y,color=db['plot_color'],label=(r"$y=%0.3f%+0.3f*x$"%(b,a)  + '\n' + (r"$R^2=%0.3f$"%(pl.square(pl.corrcoef(x,y)[0,1]),))),lw=2)
            sp.set_xlim(0,90)
            sp.set_xlabel("Age (years)", fontsize=14)
            sp.set_ylabel(ylabel, fontsize=14)
            sp.grid()
            sp.legend(loc='best',prop={'size':12})

            sp_sum.plot(_x,_y,color=db['plot_color'],label=db['diagnosis'])
            sp_sum.set_xlabel("Age (years)", fontsize=14)
            sp_sum.set_ylabel(ylabel, fontsize=14)
            sp_sum.grid()
        _draw_sp(sp1,sp1_sum,record_data[:,0],record_data[:,1],r"$\beta$")
        _draw_sp(sp2,sp2_sum,record_data[:,0],record_data[:,2],r"$\sigma$")
        _draw_sp(sp3,sp3_sum,record_data[:,0],record_data[:,3],r"$\sigma/\bar{x}$")
        _draw_sp(sp4,sp4_sum,record_data[:,0],record_data[:,4],r"$\bar{x}$")

        fig.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.06, wspace=0.25, hspace=.2)
        fig.savefig("metric_age_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

    fig_sum.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.06, wspace=0.25, hspace=.2)
    sp1_sum.legend(loc='best',prop={'size':12})
    sp2_sum.legend(loc='best',prop={'size':12})
    sp3_sum.legend(loc='best',prop={'size':12})
    sp4_sum.legend(loc='best',prop={'size':12})
    fig_sum.savefig("metric_age.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    
    pl.close()


if __name__ == '__main__':
    global config
    config = common.load_config()
    preview = True
    # plot_homeostasis(preview)
    # boxplot_diagnosis(preview)
    # #plot_clusters(preview)
    # plot_homeostasis_interp(preview)
    plot_homeostasis_median(preview)
    # plot_homeostasis_contour(preview)
    #plot_homeostasis_contour_median(preview)
    #plot_metric_median_age(preview)
    #plot_metric_age(preview)