#!/usr/bin/env python
import pylab as pl
from pprint import pprint
from scipy.cluster.vq import kmeans2
import itertools
import common

def draw(in1,in2=None,label="",show=False):  
    pl.figure(facecolor='white')
    if in2 is None:
        pl.plot(in1,label=label)
    else:
        pl.plot(in1,in2,label=label)
    if show:
        pl.show()
    exit(0)

def plot_hrv(hrv, hrv_interp, preview=False):
    """
    Plot HRV and it's approximation.
    In:
        hrv : list, [time, hrv], ndarray,ndarray
            Time & hrv vectors
        hrv_interp : list, [time_interp, hrv_interp], ndarray,ndarray
            Interpolated time & hrv vectors
    """
    fig = pl.figure("signals")
    pl.plot(hrv[0], hrv[1], 'or')
    pl.ylabel('HRV (s)') 
    pl.plot(hrv_interp[0], hrv_interp[1], 'k')
    pl.xlabel('Time (s)')
    if preview:
        pl.show()
    pl.close()


def plot_beta_cv(results, preview=False):
    x = [x['cov'] for x in results]
    y = [y['beta'] for y in results]
    pl.figure(figsize=(6, 6), facecolor='white')
    pl.title(r'%(record)s: $\sigma/\bar{x}$' % results[0])
    pl.plot(x,y,'.r')
    pl.axhline(y=1,color='b')
    pl.ylim(0,2)
    pl.xlabel(r"$\sigma/\bar{x}$")
    pl.ylabel(r"$\beta$")
    pl.grid()
    pl.savefig("cv_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:  
        pl.show()
    pl.close()


def plot_beta_std(results, preview=False):
    x = [x['std'] for x in results]
    y = [y['beta'] for y in results]
    pl.figure(figsize=(6, 6), facecolor='white')
    pl.title(r'%(record)s: $\beta/\sigma$' % results[0])
    pl.plot(x,y,'.r')
    pl.ylim(0,2)
    pl.xlim(0,140)
    pl.axhline(y=1,color='b')
    pl.axvline(x=70,color='b')
    pl.xlabel(r"$\sigma$, ms")
    pl.ylabel(r"$\beta$")
    pl.grid()
    pl.savefig("std_%(record)s.png" % results[0],facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_beta(freq, fft, aprox_y, result, preview=False):
    pl.figure(figsize=(6, 6),facecolor='white')
    pl.plot(freq,fft,'b')
    pl.plot(freq,aprox_y,'r')
    pl.title('%(record)s: %(time_from)s - %(time_to)s' % result)
    pl.xlim(min(freq),max(freq))
    pl.xlabel(r"$\lg f$")
    pl.ylabel(r"$\lg S$")
    pl.grid()
    if result['frag'] == 1:
        pl.savefig("beta_%(record)s_%(frag)s.png" % result,facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()
    pl.close()


def plot_homeostasis(preview=False):
    pl.figure("stochastic_homeostasis_summary", figsize=(12, 6), facecolor='white')
    pl.title(r'$\beta$ by desease')
    pl.suptitle('Stochastic homeostasis by diagnosis', fontsize=20) 
    
    sp_beta_std = pl.subplot(121)
    pl.ylim(0,2)
    pl.xlim(0,140)
    pl.axhline(y=1,color='k')
    pl.axvline(x=70,color='k')
    pl.xlabel(r"$\sigma\ (ms)$", fontsize=18)
    pl.ylabel(r"$\beta$", fontsize=18)
    pl.grid()

    sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
    pl.xlim(0,50)
    pl.axhline(y=1,color='k')
    pl.xlabel(r"$\sigma/\bar{x}$", fontsize=18)
    pl.grid()

    for db in common.SIGNALS:
        record_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_stat[0].extend(signal['beta'])
                record_stat[1].extend(signal['std'])
                record_stat[2].extend(signal['cov'])
        sp_beta_std.scatter(record_stat[1],record_stat[0],label=db['diagnosis'],color=db['plot_color'])
        sp_beta_cov.scatter(record_stat[2],record_stat[0],label=db['diagnosis'],color=db['plot_color'])

    pl.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.1)
    sp_beta_cov.legend(loc='best')
    pl.savefig("stochastic_homeostasis_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def boxplot_diagnosis(preview=False):
    legend = []
    by_diagnosis = [[],[],[],[]]
    #betas, stds, covs, means = [], [], [], []
    for db in common.SIGNALS:
        by_record = [[],[],[],[]]
        records = db['records'].split()
        for record in records:
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                by_record[0].append(signal['beta'])
                by_record[1].append(signal['std'])
                by_record[2].append(signal['cov'])
                by_record[3].append(signal['mean'])
        fig = pl.figure("boxplot_%s" % db['diagnosis'], figsize=(12, 12), facecolor='white')
        fig.suptitle('Statistics for %s database' % db['diagnosis'], fontsize=20)        
        sp_beta = pl.subplot(221)
        sp_beta.boxplot(by_record[0])
        pl.xticks(range(1,len(records)+1),records, rotation='vertical')
        pl.ylabel(r"$\beta$", fontsize=18)
        sp_std = pl.subplot(222,sharex=sp_beta)
        sp_std.boxplot(by_record[1])
        pl.xticks(range(1,len(records)+1),records, rotation='vertical')
        pl.ylabel(r"$\sigma$", fontsize=18)
        sp_cov = pl.subplot(223,sharex=sp_beta)
        sp_cov.boxplot(by_record[2])
        pl.xticks(range(1,len(records)+1),records, rotation='vertical')
        pl.ylabel(r"$\sigma/\bar{x}$", fontsize=18)
        sp_mean = pl.subplot(224,sharex=sp_beta)
        sp_mean.boxplot(by_record[3])
        pl.xticks(range(1,len(records)+1),records, rotation='vertical')
        pl.ylabel(r"$\bar{x}$", fontsize=18)
        pl.subplots_adjust(left=0.06, right=0.95, top=0.94, bottom=0.1)
        pl.savefig("boxplot_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()

        for i in range(4):
            by_diagnosis[i].append(list(itertools.chain.from_iterable(by_record[i])))
        legend.append(db['diagnosis'][:2])

    fig = pl.figure("boxplot_diagnosis_summary", figsize=(12, 12), facecolor='white')
    fig.suptitle('Summary statistics for diagnoses', fontsize=20)        
    sp_beta = pl.subplot(221)
    sp_beta.boxplot(by_diagnosis[0])
    pl.xticks(range(1,len(common.SIGNALS)+1),legend, rotation='vertical')
    pl.ylabel(r"$\beta$", fontsize=18)
    sp_std = pl.subplot(222,sharex=sp_beta)
    sp_std.boxplot(by_diagnosis[1])
    pl.xticks(range(1,len(common.SIGNALS)+1),legend, rotation='vertical')
    pl.ylabel(r"$\sigma$", fontsize=18)
    sp_cov = pl.subplot(223,sharex=sp_beta)
    sp_cov.boxplot(by_diagnosis[2])
    pl.xticks(range(1,len(common.SIGNALS)+1),legend, rotation='vertical')
    pl.ylabel(r"$\sigma/\bar{x}$", fontsize=18)
    sp_mean = pl.subplot(224,sharex=sp_beta)
    sp_mean.boxplot(by_diagnosis[3])
    pl.xticks(range(1,len(common.SIGNALS)+1),legend, rotation='vertical')
    pl.ylabel(r"$\bar{x}$", fontsize=18)
    pl.subplots_adjust(left=0.06, right=0.95, top=0.94, bottom=0.1)
    pl.savefig("boxplot_diagnosis_summary.png",facecolor='w',edgecolor='k',transparent=True)
    if preview:
        pl.show()
    pl.close()


def plot_clusters(preview=False):
    pl.figure("clusterisation", figsize=(12, 6), facecolor='white')
    pl.title('Clusterisation by desease')
    pl.suptitle('Clasterisation', fontsize=20) 
    
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
    for db in common.SIGNALS:
        record_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                stat[0].extend(signal['beta'])
                stat[1].extend(signal['std'])
                stat[2].extend(signal['cov'])
    
    pl.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.1)
    pl.savefig("clusterisation.png",facecolor='w',edgecolor='k',transparent=True)

    std, beta = common.filter2d(stat[1], stat[0], 1)
    res, idx = kmeans2(pl.array(zip(std, beta)),len(common.SIGNALS))
    colors = ([([0.4,1,0.4],[1,0.4,0.4],[0.1,0.8,1])[i] for i in idx])
    sp_beta_std.scatter(std, beta, c=colors)
    sp_beta_std.scatter(res[:,0],res[:,1], marker='o', s = 500, linewidths=2, c='none')
    sp_beta_std.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2, c='k')

    cov, beta = common.filter2d(stat[2], stat[0], 1)
    res, idx = kmeans2(pl.array(zip(cov, beta)),len(common.SIGNALS))
    colors = ([([0.4,1,0.4],[1,0.4,0.4],[0.1,0.8,1])[i] for i in idx])
    sp_beta_cov.scatter(cov, beta, c=colors)    
    sp_beta_cov.scatter(res[:,0],res[:,1], marker='o', s = 500, linewidths=2, c='none')
    sp_beta_cov.scatter(res[:,0],res[:,1], marker='x', s = 500, linewidths=2, c='k')

    if preview:
        pl.show()
    pl.close()


def plot_homeostasis_interp(preview=False):
    for db in common.SIGNALS:
        record_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                record_stat[0].extend(signal['beta'])
                record_stat[1].extend(signal['std'])
                record_stat[2].extend(signal['cov'])
        fig = pl.figure("stochastic_homeostasis_summary_interp", figsize=(12, 6), facecolor='white')
        fig.suptitle('Stochastic homeostasis of %s' % (db['diagnosis'],), fontsize=20) 

        sp_beta_std = pl.subplot(121)
        #pl.title('Stochastic homeostasis for %s' % (db['diagnosis'],), fontsize=18) 
        pl.ylim(0,2)
        pl.xlim(0,140)
        pl.axhline(y=1,color='k')
        pl.axvline(x=70,color='k')
        pl.xlabel(r"$\sigma\ (ms)$", fontsize=18)
        pl.ylabel(r"$\beta$", fontsize=18)
        pl.grid()
        #sp_beta_std.hist2d(record_stat[1],record_stat[0],range=[[0,140],[0,2]],bins=[140/2,200/2],cmin=1,label=db['diagnosis'])
        s = sp_beta_std.hexbin(record_stat[1],record_stat[0],extent=[0,140,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
        s.set_clim()
        cb = pl.colorbar(s)

        sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
        pl.xlim(0,50)
        pl.axhline(y=1,color='k')
        pl.xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        pl.grid()
        #sp_beta_std.hist2d(record_stat[1],record_stat[0],range=[[0,140],[0,2]],bins=[140/2,200/2],cmin=1,label=db['diagnosis'])
        s = sp_beta_cov.hexbin(record_stat[2],record_stat[0],extent=[0,50,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
        s.set_clim()
        cb = pl.colorbar(s)
        #cb.set_label('log10(N)')
        
        #sp_beta_std.scatter(record_stat[1],record_stat[0],label=db['diagnosis'],color=db['plot_color'])
        #sp_beta_cov.scatter(record_stat[2],record_stat[0],label=db['diagnosis'],color=db['plot_color'])

        pl.tight_layout()
        pl.subplots_adjust(left=0.06, right=0.97, top=0.9, bottom=0.1)
        pl.savefig("stochastic_homeostasis_diagnosis_interp_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()

    pl.close()


if __name__ == '__main__':
    #plot_homeostasis(True)
    #boxplot_diagnosis()
    #plot_clusters(True)
    plot_homeostasis_interp(True)