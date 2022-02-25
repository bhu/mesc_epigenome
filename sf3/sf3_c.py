from deeptools import heatmapper
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
k = 4000
samps = ['Serum mESC', '2i mESC']
for reg in ['Met', 'Unm']:
    hm = heatmapper.heatmapper()
    hm.read_matrix_file('../aggr/mESC_K27me3_CPM_serum2i.%s.mat.gz' % reg)
    ms = pd.DataFrame(hm.matrix.matrix)
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
    zma = np.nanquantile(ms.values, .98)
    zmi = np.nanquantile(ms.values, 0)
    fig = plt.figure(constrained_layout=True,figsize=[5,4])
    widths = np.append(np.repeat(1, len(samps)), .1).tolist()
    heights = [1, 2.4]
    gs = fig.add_gridspec(ncols=len(samps)+1, nrows=2,
                          width_ratios=widths, height_ratios=heights)
    for i in np.arange(len(samps)):
        samp = samps[i]
        mat = ms.iloc[:,(i*k):(i*k+k)].dropna(axis=0,thresh=int(k/5))
        aggr = pd.DataFrame({'x': np.arange(k) + 1, 'y': mat.mean(axis=0)})
        ax = fig.add_subplot(gs[0, i])
        #clr = clrs[samp.split(' ')[0]]
        ax.plot(aggr['x'], aggr['y'])
        ax.set_title(samp)
        ax.set_ylim([ymi,yma])
        ax.set_xticks([])
        if i == 0:
            ax.set_ylabel('H3K27me3 (CPM)')
        else:
            ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticklabels([])
        ax.yaxis.grid(False)
        ax = fig.add_subplot(gs[1, i])
        ax.set_yticks([])
        ax.set_xticks([0,2000,3999])
        if reg == 'Met':
            ax.set_xticklabels(['-20kb','Methyl\'d\nCGI','+20kb'])
        elif reg == 'Unm':
            ax.set_xticklabels(['-20kb','Unmethyl\'d\nCGI','+20kb'])

        im = ax.imshow(mat.values, aspect='auto', cmap='magma',
                       vmin = zmi, vmax = zma)
    cax = fig.add_subplot(gs[:,-1])
    fig.colorbar(im, cax=cax, orientation="vertical")
    plt.savefig('sf3_b.%s.pdf' % reg, transparent=True,dpi=600)
    plt.close()
