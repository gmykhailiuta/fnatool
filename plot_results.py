#!/usr/bin/env python
import matplotlib as mp
#mp.use('GtkAgg')
import pylab as pl
import matplotlib.mlab as mlab
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
    Plot pure, filtered and interpolated HRV.
    In:
        hrv : list, [time, hrv], ndarray,ndarray, Time & hrv vectors
        hrv_filt : list, [time_filt, hrv_filt], ndarray,ndarray, Filtered time & hrv vectors
        hrv_interp : list, [time_interp, hrv_interp], ndarray,ndarray, Interpolated time & hrv vectors
        record : str, record name in WFDB format
        preview : show plots while running
    """
    fig = pl.figure("signals", figsize=(22, 18), facecolor='white')
    pl.suptitle(record)
    sp_orig = pl.subplot(311)
    pl.title("Original signal")
    pl.ylabel(r"$RR_i$ (ms)") 
    pl.plot(hrv[0], hrv[1], color='k', marker='.', linestyle='--', alpha=.5)
    sp_filt = pl.subplot(312,sharex=sp_orig,sharey=sp_orig)
    pl.title("Filtered signal")
    pl.ylabel("$RR_i$ (ms)")
    pl.plot(hrv_filt[0], hrv_filt[1], color='k', marker='.', linestyle='--', alpha=.5)
    sp_filt = pl.subplot(313,sharex=sp_orig,sharey=sp_orig)
    pl.title("Interpolated signal")
    pl.ylabel(r"$RR_j$ (ms)") 
    pl.xlabel('Time (s)')
    pl.plot(hrv_interp[0], hrv_interp[1], color='k', marker='.', alpha=.5)
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


def plot_homeostasis_record(results, preview=False):
    std = [x['std'] for x in results]
    cov = [x['cov'] for x in results]
    beta = [y['beta'] for y in results]
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num="stochastic_homeostasis_record")
    fig.suptitle('Stochastic homeostasis for record %(record)s' % results[0], fontsize=20) 
    sp1.scatter(std,beta,color='b')
    sp1.scatter(pl.mean(std),pl.mean(beta),s=400,color='r',marker='o',alpha=.7)
    sp1.scatter(pl.mean(std),pl.mean(beta), marker='x', s=200, c='k',lw=2)
    sp1.axhline(y=1,color='k',linestyle='--')
    sp1.axvline(x=70,color='k',linestyle='--')
    sp1.set_ylim(0,2)
    sp1.set_xlim(0,140)
    sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
    sp1.set_ylabel(r"$\beta$", fontsize=20)
    sp1.grid()
    sp2.scatter(cov,beta,color='b')
    sp2.scatter(pl.mean(cov),pl.mean(beta),s=400,color='r',marker='o',alpha=.7)
    sp2.scatter(pl.mean(cov),pl.mean(beta), marker='x', s=200, c='k',lw=2)
    sp2.axhline(y=1,color='k',linestyle='--')
    sp2.set_xlim(0,20)
    sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
    sp2.grid()
    fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    fig.savefig("homeostasis_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_beta_hist(results, preview=False):
    beta = [y['beta'] for y in results]
    fig, sp = pl.subplots(1,1,sharey=False,sharex=False,figsize=(6, 6),facecolor='white',num="beta_hist_record")
    fig.suptitle('Beta histogram for record %(record)s' % results[0], fontsize=20) 
    sp.hist(beta, bins=10, normed=1, facecolor='green', alpha=0.75)
    sp.set_xlabel(r"$\beta$", fontsize=20)
    sp.set_ylabel(r"$Density$", fontsize=20)
    sp.set_xlim(0,2)
    sp.grid()
    fig.savefig("hist_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.show()
    pl.close()


def plot_beta(freq, fft, a, b, result, preview=False):
    y = a * freq + b
    pl.figure(figsize=(6, 6),facecolor='white')
    pl.plot(freq,fft,color='k',label="PSD")
    pl.plot(freq,y,color='r',ls='--',lw=2,label=(r"$y=%0.3f%+0.3f*x$"%(b,a)  + '\n' + (r"$R^2=%0.3f$"%(pl.square(pl.corrcoef(freq,fft)[0,1]),))))
    pl.title('%(record)s: %(time_beg)s' % result)
    pl.xlim(min(freq),max(freq))
    pl.xlabel(r"$\lg f$", fontsize=18)
    pl.ylabel(r"$\lg S$", fontsize=18)
    pl.grid()
    pl.legend(loc='best')
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
    sp1_sum.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
    sp1_sum.set_ylabel(r"$\beta$", fontsize=20)
    sp1_sum.grid()
    sp2_sum.axhline(y=1,color='k',linestyle='--')
    sp2_sum.set_xlim(0,20)
    sp2_sum.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
    sp2_sum.grid()

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

        for i in range(3):
            diag_data[i].append(pl.array(record_data[i]))

        fig_diag, (sp1_diag, sp2_diag) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num=("stochastic_homeostasis_%s" % db['diagnosis']))
        fig_diag.suptitle('Stochastic homeostasis of diagnosis %s' % db['diagnosis'], fontsize=20)
        
        sp1_diag.scatter(record_data[1],record_data[0],label=db['diagnosis'],color=db['plot_color'],alpha=.5,marker='x')
        sp1_diag.scatter(pl.mean(record_data[1]),pl.mean(record_data[0]), marker='o', s=100, c='w', alpha=0.5, lw=2)
        sp1_diag.scatter(pl.mean(record_data[1]),pl.mean(record_data[0]), marker='x', s=500, c='k', lw=2)
        sp1_diag.axhline(y=1,color='k',linestyle='--')
        sp1_diag.axvline(x=70,color='k',linestyle='--')
        sp1_diag.set_ylim(0,2)
        sp1_diag.set_xlim(0,140)
        sp1_diag.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
        sp1_diag.set_ylabel(r"$\beta$", fontsize=20)
        sp1_diag.grid()
        
        sp2_diag.scatter(record_data[2],record_data[0],label=db['diagnosis'],color=db['plot_color'],alpha=.5,marker='x')
        sp2_diag.scatter(pl.mean(record_data[2]),pl.mean(record_data[0]), marker='o', s=100, c='w', alpha=0.5, lw=2)
        sp2_diag.scatter(pl.mean(record_data[2]),pl.mean(record_data[0]), marker='x', s=500, c='k', lw=2)
        sp2_diag.axhline(y=1,color='k',linestyle='--')
        sp2_diag.set_xlim(0,20)
        sp2_diag.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
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
    for db in config['SIGNALS']:
        record_stat = [[],[],[],[]]
        records = db['records'].split()
        for record in records:
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_stat[0].append(signal['beta'] * common.filter1d(signal['beta'], algos=['5per95']))
                record_stat[1].append(signal['std'] * common.filter1d(signal['std'], algos=['5per95']))
                record_stat[2].append(signal['cov'] * common.filter1d(signal['cov'], algos=['5per95']))
                record_stat[3].append(signal['mean'] * common.filter1d(signal['mean'], algos=['5per95']))
        fig, (sp1,sp2,sp3,sp4) = pl.subplots(4,1,sharey=False,sharex=True,figsize=(18, 18),facecolor='white',num="boxplot_diagnosis_%s" % db['diagnosis'])
        fig.suptitle('Statistics for diagnisis %s' % db['diagnosis'], fontsize=20)
        xticks = range(1,len(records)+1)
        def _draw_sp(sp,y,ylabel):
            sp.set_ylabel(ylabel, fontsize=22)
            sp.boxplot(y)
            sp.scatter(xticks, [pl.mean(x) for x in y], marker='x',color='k')
            sp.scatter(xticks, [pl.std(x) for x in y], marker='o',color='k')
            sp.set_xticks(xticks)
            sp.set_xticklabels(records, rotation='vertical')
            sp.yaxis.grid()

        _draw_sp(sp1,record_stat[0],r"$\beta$")
        sp1.axhline(y=1,color='g',linestyle='--')
        sp1.set_ylim(-.2,2)
        _draw_sp(sp2,record_stat[1],r"$\sigma$")
        sp2.axhline(y=70,color='g',linestyle='--')
        sp2.set_ylim(10,140)
        _draw_sp(sp3,record_stat[2],r"$\sigma/\bar{x}$")
        sp3.set_ylim(0,20)
        _draw_sp(sp4,record_stat[3],r"$\bar{x}$")
        sp4.set_ylim(400,1200)
        sp4.set_xlabel(r"$Record$", fontsize=22)

        fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.12, wspace=0.15, hspace=0.05)
        fig.savefig("boxplot_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

        for i in range(4):
            diagnosis_stat[i].append(list(itertools.chain.from_iterable(record_stat[i])))
        legend.append(db['diagnosis'][:3])

    for s in range(len(diagnosis_stat)):
        for r in range(len(diagnosis_stat[s])):
            diagnosis_stat[s][r] *= common.filter1d(diagnosis_stat[s][r], algos=['5per95'])
    
    fig_sum, ((sp1_sum,sp2_sum),(sp3_sum,sp4_sum)) = pl.subplots(2,2,sharey=False,sharex=True,figsize=(12, 12),facecolor='white',num="boxplot_diagnosis_summary")
    fig_sum.suptitle('Statistics by diagnisis', fontsize=20)
    xticks = range(1,len(config['SIGNALS'])+1)
    def _draw_sp(sp,y,ylabel):
        sp.set_ylabel(ylabel, fontsize=22)
        sp.boxplot(y)
        sp.scatter(xticks, [pl.mean(x) for x in y], marker='x',color='k')
        sp.scatter(xticks, [pl.std(x) for x in y], marker='o',color='k')
        sp.set_xticks(xticks)
        sp.set_xticklabels(legend, rotation='vertical')
        sp.yaxis.grid()

    _draw_sp(sp1_sum,diagnosis_stat[0],r"$\beta$")
    sp1_sum.axhline(y=1,color='g',linestyle='--')
    sp1_sum.set_ylim(-.2,2)
    _draw_sp(sp2_sum,diagnosis_stat[1],r"$\sigma$")
    sp2_sum.axhline(y=70,color='g',linestyle='--')
    sp2_sum.set_ylim(0,140)
    _draw_sp(sp3_sum,diagnosis_stat[2],r"$\sigma/\bar{x}$")
    sp3_sum.set_ylim(0,20)
    sp3_sum.set_xlabel(r"$Diagnosis$", fontsize=20)
    _draw_sp(sp4_sum,diagnosis_stat[3],r"$\bar{x}$")
    sp4_sum.set_ylim(400,1200)
    sp4_sum.set_xlabel(r"$Diagnosis$", fontsize=20)

    fig_sum.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.1, wspace=0.2, hspace=0.05)
    fig_sum.savefig("boxplot_diagnosis_summary.png",facecolor='w',edgecolor='k',transparent=True)
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
        sp1.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
        sp1.set_ylabel(r"$\beta$", fontsize=20)
        sp1.grid()
        sp2.axhline(y=1,color='k',linestyle='--')
        sp2.set_xlim(0,20)  
        sp2.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
        sp2.grid()

        x, y = record_stat[1], record_stat[0]
        s1 = sp1.hexbin(x, y,extent=[0,140,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
        s1.set_clim()
        fig.colorbar(s1,ax=sp1)
        x, y = record_stat[2], record_stat[0]
        s1 = sp2.hexbin(x, y,extent=[0,20,0,2],gridsize=20,mincnt=2,label=db['diagnosis'])
        s1.set_clim()
        fig.colorbar(s1,ax=sp2)

        fig.tight_layout()
        fig.subplots_adjust(left=0.06, right=0.995, top=0.9, bottom=0.1, wspace=0.01, hspace=.01)
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
    sp1_sum.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
    sp1_sum.set_ylabel(r"$\beta$", fontsize=20)
    sp1_sum.grid()
    sp2_sum.axhline(y=1,color='k',linestyle='--')
    sp2_sum.set_xlim(0,10)
    sp2_sum.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
    sp2_sum.grid()

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

        fig_diag, (sp1_diag, sp2_diag) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(12, 6),facecolor='white',num=("stochastic_homeostasis_median_%s" % db['diagnosis']))
        fig_diag.suptitle('Stochastic homeostasis: medians by record of diagnosis %s' % db['diagnosis'], fontsize=20)
        
        sp1_diag.scatter(record_median[1],record_median[0],color=db['plot_color'],label='Data')
        sp1_diag.scatter(pl.mean(record_median[1]),pl.mean(record_median[0]),s=500,marker='x', c='k',lw=2,label='Centre of mass')
        sp1_diag.axhline(y=1,color='k',linestyle='--')
        sp1_diag.axvline(x=70,color='k',linestyle='--')
        sp1_diag.set_ylim(0,2)
        sp1_diag.set_xlim(0,140)
        sp1_diag.set_xlabel(r"$\sigma\ (ms)$", fontsize=20)
        sp1_diag.set_ylabel(r"$\beta$", fontsize=20)
        sp1_diag.grid()
        
        sp2_diag.scatter(record_median[2],record_median[0],color=db['plot_color'],label='Data')
        sp2_diag.scatter(pl.mean(record_median[2]),pl.mean(record_median[0]),s=500,marker='x', c='k',lw=2,label='Centre of mass')
        sp2_diag.axhline(y=1,color='k',linestyle='--')
        sp2_diag.set_xlim(0,20)
        sp2_diag.set_xlabel(r"$\sigma/\bar{x}$", fontsize=20)
        sp2_diag.grid()

        fig_diag.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.1)
        fig_diag.savefig("stochastic_homeostasis_median_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)

        for i in range(3):
            diag_data_1 = pl.array(list(itertools.chain.from_iterable(record_data[i])))
            diag_data[i].append(diag_data_1)
            diag_median[i] = pl.median(diag_data_1)
        
        sp1_sum.scatter(diag_median[1],diag_median[0],s=200,label=db['diagnosis'],color=db['plot_color'],marker='o')
        sp2_sum.scatter(diag_median[2],diag_median[0],s=200,label=db['diagnosis'],color=db['plot_color'],marker='o')

    fig_sum.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.1)
    sp2_sum.legend(loc='best')
    fig_sum.savefig("stochastic_homeostasis_medians_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    
    pl.close()


def plot_homeostasis_contour(diags_filter=('HYP','CHF'),preview=False):
    global config
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(24, 12),facecolor='white',num="stochastic_homeostasis_contours_%s" % ('_'.join(diags_filter),))
    fig.suptitle('Stochastic homeostasis: contour of diagnoses %s' % (','.join(diags_filter),), fontsize=36)
    for i in range(len(config['SIGNALS'])):
        db = config['SIGNALS'][i]
        if db['diagnosis'] not in diags_filter:
            continue
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

    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.2)
    fig.savefig("stochastic_homeostasis_contours_diagnoses_%s.png" % ('_'.join(diags_filter),),facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_homeostasis_contour_median(diags_filter=('NM1','HYP'),preview=False):
    global config
    fig, (sp1, sp2) = pl.subplots(1,2,sharey=True,sharex=False,figsize=(24, 12),facecolor='white',num='stochastic_homeostasis_contours_median_diags')
    fig.suptitle('Stochastic homeostasis: median contour of diagnoses %s' % (','.join(diags_filter),), fontsize=32)
    for db in config['SIGNALS']:
        if db['diagnosis'] not in diags_filter:
            continue
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

    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=.2)
    fig.savefig("stochastic_homeostasis_contours_median_diagnoses_%s.png" % ('_'.join(diags_filter),),facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_metric_median_age(preview=False):
    global config
    
    fig_sum, ((sp1_sum, sp2_sum), (sp3_sum, sp4_sum)) = pl.subplots(2,2,sharey=False,sharex=False,figsize=(9, 9),facecolor='white',num=("metric_median_age"))
    fig_sum.suptitle('Metric by age : medians by diagnosis', fontsize=20) 

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
    preview = config['PREVIEW']
    if config['DEBUG']:
        #plot_homeostasis(preview)
        boxplot_diagnosis(preview)
        # #plot_clusters(preview)
        #plot_homeostasis_interp(preview)
        #plot_homeostasis_median(preview)
        #plot_homeostasis_contour(preview=preview)
        #plot_homeostasis_contour_median(preview=preview)
        #plot_metric_median_age(preview)
        #plot_metric_age(preview)
    else:
        plot_homeostasis(preview)
        boxplot_diagnosis(preview)
        #plot_clusters(preview)
        plot_homeostasis_interp(preview)
        plot_homeostasis_median(preview)
        plot_homeostasis_contour(preview)
        plot_homeostasis_contour_median(preview)
        plot_metric_median_age(preview)
        plot_metric_age(preview)