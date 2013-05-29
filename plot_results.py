#!/usr/bin/env python
import pylab as pl
from pprint import pprint
from scipy.cluster.vq import kmeans2
from warnings import warn
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

def plot_hrv(hrv, hrv_interp, record, preview=False):
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
    pl.title(record)
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
    pl.xlim(0,50)
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
    global config
    #legend = []
    diagnosis_stat = [[],[],[]]
    for db in config['SIGNALS']:
        records_stat = [[],[],[]]
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                records_stat[0].extend(signal['beta'])
                records_stat[1].extend(signal['std'])
                records_stat[2].extend(signal['cov'])

        fig = pl.figure("stochastic_homeostasis_%s" % db['diagnosis'], figsize=(12, 6), facecolor='white')
        fig.suptitle('Stochastic homeostasis of %s database' % db['diagnosis'], fontsize=20)        
        sp_beta_std = pl.subplot(121)
        pl.ylim(0,2)
        pl.xlim(0,140)
        pl.axhline(y=1,color='k')
        pl.axvline(x=70,color='k')
        pl.xlabel(r"$\sigma\ (ms)$", fontsize=18)
        pl.ylabel(r"$\beta$", fontsize=18)
        pl.grid()
        #print len(records_stat[0]),len(records_stat[1]),len(records_stat[2])
        std, beta = common.filter2d(records_stat[1], records_stat[0], axes=['x','y'], algos=['5per95'])
        #std, beta = records_stat[1], records_stat[0]
        sp_beta_std.scatter(std,beta,label=db['diagnosis'],color=db['plot_color'])

        sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
        pl.xlim(0,50)
        pl.axhline(y=1,color='k')
        pl.xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        pl.grid()
        cov, beta = common.filter2d(records_stat[2], records_stat[0], axes=['x','y'], algos=['5per95'])
        #cov, beta = records_stat[2], records_stat[0]
        sp_beta_cov.scatter(cov,beta,label=db['diagnosis'],color=db['plot_color'])

        pl.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
        pl.savefig("stochastic_homeostasis_diagnosis_%s.png" % (db['diagnosis'],),facecolor='w',edgecolor='k',transparent=True)
        if preview:
            pl.show()
        #pprint(records_stat[0][:10])
        for i in range(3):
            diagnosis_stat[i].extend(records_stat[i])
        #legend.append(db['diagnosis'][:3])
        #print db['diagnosis'][:2]
        

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
    for db in config['SIGNALS']:
        std, beta = common.filter2d(diagnosis_stat[1], diagnosis_stat[0], axes=['x','y'], algos=['5per95'])
        sp_beta_std.scatter(diagnosis_stat[1],diagnosis_stat[0],label=db['diagnosis'],color=db['plot_color'],marker='x')
        cov, beta = common.filter2d(diagnosis_stat[2], diagnosis_stat[0], axes=['x','y'], algos=['5per95'])
        sp_beta_cov.scatter(cov,beta,label=db['diagnosis'],color=db['plot_color'],marker='x')

    pl.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.2, hspace=.3)
    sp_beta_cov.legend(loc='best')
    pl.savefig("stochastic_homeostasis_diagnosis.png",facecolor='w',edgecolor='k',transparent=True)
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
        fig = pl.figure("boxplot_%s" % db['diagnosis'], figsize=(12, 12), facecolor='white')
        fig.suptitle('Statistics for %s database' % db['diagnosis'], fontsize=20)
        xticks = range(1,len(records)+1)
        sp_beta = pl.subplot(221)
        pl.title(r"$\beta$", fontsize=18)
        #pprint(record_stat[0])
        sp_beta.boxplot(record_stat[0])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[0]], marker='x',color='k')
        #pl.ylim(0,pl.percentile())
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(0,2)
        sp_beta.yaxis.grid()
        sp_std = pl.subplot(222,sharex=sp_beta)
        pl.title(r"$\sigma$", fontsize=18)
        sp_std.boxplot(record_stat[1])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[1]], marker='x',color='k')
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(0,400)
        sp_std.yaxis.grid()
        sp_cov = pl.subplot(223,sharex=sp_beta)
        pl.title(r"$\sigma/\bar{x}$", fontsize=18)
        sp_cov.boxplot(record_stat[2])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[2]], marker='x',color='k')
        pl.xticks(xticks,records, rotation='vertical')
        pl.ylim(0,50)
        sp_cov.yaxis.grid()
        sp_mean = pl.subplot(224,sharex=sp_beta)
        pl.title(r"$\bar{x}$", fontsize=18)
        sp_mean.boxplot(record_stat[3])
        pl.scatter(xticks, [pl.mean(x) for x in record_stat[3]], marker='x',color='k')
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

    #for i in range(len(diagnosis_stat)):
     #   diagnosis_stat_filt[i] = common.filter2d(diagnosis_stat[i], diagnosis_stat[i], axes=['x'], algos=['5per95'])
    for s in range(len(diagnosis_stat)):
        for r in range(len(diagnosis_stat[s])):
            diagnosis_stat[s][r] *= common.filter1d(diagnosis_stat[s][r], algos=['5per95'])
    fig = pl.figure("boxplot_diagnosis_summary", figsize=(12, 12), facecolor='white')
    fig.suptitle('Summary statistics for diagnoses', fontsize=20)   
    xticks = range(1,len(config['SIGNALS'])+1)
    sp_beta = pl.subplot(221)
    pl.title(r"$\beta$", fontsize=18)
    sp_beta.boxplot(diagnosis_stat[0])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[0]], marker='x',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(0,2)
    sp_beta.yaxis.grid()
    pl.ylabel(r"$\beta$", fontsize=18)
    sp_std = pl.subplot(222,sharex=sp_beta)
    pl.title(r"$\sigma$", fontsize=18)
    sp_std.boxplot(diagnosis_stat[1])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[1]], marker='x',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(0,400)
    sp_std.yaxis.grid()
    sp_cov = pl.subplot(223,sharex=sp_beta)
    pl.title(r"$\sigma/\bar{x}$", fontsize=18)
    sp_cov.boxplot(diagnosis_stat[2])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[2]], marker='x',color='k')
    pl.xticks(xticks,legend, rotation='vertical')
    pl.ylim(0,50)
    sp_cov.yaxis.grid()
    sp_mean = pl.subplot(224,sharex=sp_beta)
    pl.title(r"$\bar{x}$", fontsize=18)
    sp_mean.boxplot(diagnosis_stat[3])
    pl.scatter(xticks, [pl.mean(x) for x in diagnosis_stat[3]], marker='x',color='k')
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
    pl.figure("clusterization", figsize=(12, 6), facecolor='white')
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
        std, beta = common.filter2d(record_stat[1], record_stat[0], axes=['x','y'], algos=['3per97'])
        s = sp_beta_std.hexbin(std, beta,extent=[0,140,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
        s.set_clim()
        cb = pl.colorbar(s)

        sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
        pl.xlim(0,50)
        pl.axhline(y=1,color='k')
        pl.xlabel(r"$\sigma/\bar{x}$", fontsize=18)
        pl.grid()
        #sp_beta_std.hist2d(record_stat[1],record_stat[0],range=[[0,140],[0,2]],bins=[140/2,200/2],cmin=1,label=db['diagnosis'])
        cov, beta = common.filter2d(record_stat[2], record_stat[0], axes=['x','y'], algos=['3per97'])
        s = sp_beta_cov.hexbin(cov, beta,extent=[0,50,0,2],gridsize=20,mincnt=1,label=db['diagnosis'])
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
    global config
    config = common.load_config()
    plot_homeostasis(1)
    #boxplot_diagnosis(1)
    #plot_clusters(1)
    #plot_homeostasis_interp(1)