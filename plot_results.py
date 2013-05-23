#!/usr/bin/env python
import pylab as pl
import common
from pprint import pprint
import itertools

def plot_homeostasis(preview=False):
    pl.figure("stochastic_homeostasis_summary", figsize=(12, 6), facecolor='white')
    pl.title(r'$\beta$ by desease')
    
    sp_beta_std = pl.subplot(121)
    pl.ylim(0,2)
    pl.xlim(0,140)
    pl.axhline(y=1,color='k')
    pl.axvline(x=70,color='k')
    pl.xlabel(r"$\sigma\ (ms)$")
    pl.ylabel(r"$\beta$")

    sp_beta_cov = pl.subplot(122, sharey=sp_beta_std)
    pl.axhline(y=1,color='k')
    pl.xlabel(r"$\sigma/\bar{x}$")

    legend = []
    for db in common.SIGNALS:
        beta, std, cov = [], [], []
        for record in db['records'].split():
            signal = common.read_signal("results_%s.csv" % record)
            if signal:
                beta.extend(signal['beta'])
                std.extend(signal['std'])
                cov.extend(signal['cov'])
        sp_beta_std.plot(std,beta,db['plot_params'])
        sp_beta_cov.plot(cov,beta,db['plot_params'])
        legend.append(db['diagnosis'])
    
    pl.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.1)
    pl.legend(legend, loc='best')
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
    fig.suptitle('Statistics for all databases', fontsize=20)        
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


if __name__ == '__main__':
    #plot_homeostasis()
    boxplot_diagnosis()