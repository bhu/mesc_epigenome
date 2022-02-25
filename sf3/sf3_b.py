import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import glob
import re
import seaborn as sns
import functools
import math
plt.rc('font', family='Arial')
nfacs = pd.read_csv('../data/libsz.csv')
mfacs = pd.read_csv('../../massspec.facs.csv')

clrs = {'WT': '#0C8140', 'DNMT-TKO': '#3953A4' , 'sgNSD1': '#EE2E2E'}
samps = ['WT (A1)', 'WT (J1)', 'DNMT-TKO', 'sgNSD1', ]
chips = ['mESC_WT_A1_H3K27me3', 'mESC_WT_J1_H3K27me3', 'mESC_DNMT_TKO_H3K27me3_rep2',  'mESC_NSD1_KO_H3K27me3', ]
inps = ['mESC_WT_A1_input', 'mESC_WT_J1_input', 'mESC_DNMT_TKO_input_rep2', 'mESC_NSD1_KO_input', ]
conds = ['mESC_A1_par', 'mESC_J1_par', 'mESC_J1-DNMT-TKO', 'mESC_A1_sgNSD1']
marks = ['H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3']
kind = 'enr'

k = 80

for reg in ['Met', 'Unm']:
    ext = '.cgi' + reg + '.point.mat.gz'
    nrm = 'ms'
    kind = 'lfc'
    ms = []
    for c, i, cnd, mrk in zip(chips, inps, conds, marks):
        mc = pd.read_csv('../aggr/' + c + ext, sep="\t", skiprows=1,
                         header=None, usecols=np.arange(k) + 6).values
        mi = pd.read_csv('../aggr/' + i + ext, sep="\t", skiprows=1,
                         header=None, usecols=np.arange(k) + 6).values
        n = nfacs[nfacs['samp'] == c]
        rx = ((n['endo.chip'] / n['exo.chip']) / (n['endo.input'] / n['exo.input'])).values
        nc = n['endo.chip'].values
        ni = n['endo.input'].values
        fc, fi =  min(nc, ni) / [nc, ni]
        if nrm == 'rx':
            fac = rx
        elif nrm == 'cpm':
            fac = 1
        elif nrm == 'ms':
            fac = mfacs.loc[(mfacs['mark'] == mrk) & (mfacs['type'] == cnd),'f'].values[0]
        fc = (mc * fac * fc + 1) / (mi * fi + 1)
        lfc = np.log2(fc)
        cov = mc * 1e6 * fac / nc
        if kind == 'cov':
            ms.append(pd.DataFrame(cov))
        elif kind == 'lfc':
            ms.append(pd.DataFrame(lfc))
        elif kind == 'fc':
            ms.append(pd.DataFrame(fc))
    ms = pd.concat(ms,axis=1)
    ms = ms.assign(m=ms.mean(axis=1)).sort_values('m', ascending=False).drop('m', axis=1)
    ymi = []
    yma = []
    for i in np.arange(len(samps)):
        samp = samps[i]
        mm = ms.iloc[:,(i*k):(i*k+k)].dropna(axis=0,thresh=int(k/5)).mean(axis=0)
        ymi.append(min(mm))
        yma.append(max(mm))
    ymi = min(ymi)
    yma = max(yma)
    r = .1 * (yma - ymi)
    ymi = ymi - r
    yma = yma + r
    zma = np.nanquantile(ms.values, .99)
    if kind == 'cov' or kind == 'fc':
        zmi = np.nanquantile(ms.values, 0)
    elif kind == 'lfc':
        zmi = np.nanquantile(ms.values, 0.01)
    fig = plt.figure(constrained_layout=True,figsize=[8.5,4])
    widths = np.append(np.repeat(1, len(samps)), .1).tolist()
    heights = [1, 2.4]
    gs = fig.add_gridspec(ncols=len(samps)+1, nrows=2,
                          width_ratios=widths, height_ratios=heights)
    for i in np.arange(len(samps)):
        samp = samps[i]
        mat = ms.iloc[:,(i*k):(i*k+k)].dropna(axis=0,thresh=int(k/5))
        aggr = pd.DataFrame({'x': np.arange(k) + 1, 'y': mat.mean(axis=0)})
        ax = fig.add_subplot(gs[0, i])
        clr = clrs[samp.split(' ')[0]]
        ax.plot(aggr['x'], aggr['y'], color=clr)
        ax.set_title(samp, color=clr)
        ax.set_ylim([ymi,yma])
        ax.set_xticks([])
        if i == 0:
            if kind == 'cov':
                ax.set_ylabel('H3K27me3\ncoverage')
            elif kind == 'lfc':
                ax.set_ylabel('H3K27me3\nMS norm enrichment')
            elif kind == 'fc':
                ax.set_ylabel('H3K27me3\nChIP/input')
        else:
            ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticklabels([])
        ax.yaxis.grid(False)
        ax = fig.add_subplot(gs[1, i])
        ax.set_yticks([])
        ax.set_xticks([0,40,79])
        if reg == 'Met':
            ax.set_xticklabels(['-20kb','Methyl\'d\nCGI','+20kb'])
        elif reg == 'Unm':
            ax.set_xticklabels(['-20kb','Unmethyl\'d\nCGI','+20kb'])

        im = ax.imshow(mat.values, aspect='auto', cmap='magma',
                       vmin = zmi, vmax = zma)
    cax = fig.add_subplot(gs[:,-1])
    fig.colorbar(im, cax=cax, orientation="vertical")
    plt.savefig('sf3_c.%s.%s.%s.pdf' % (reg, nrm, kind), transparent=True,dpi=600)
    plt.close()
