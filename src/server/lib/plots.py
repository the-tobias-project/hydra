# libs
import os

# third party lib
import h5py
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.cbook
import warnings
import pandas as pd
from scipy.stats import chi2

chi2sf = chi2.sf
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def qc_plots(plinkFile, outfile):
    af = []
    hz = []
    # missing_p_s = []
    missing_p_l = []
    with h5py.File(plinkFile, 'r') as store:
        for chrom in store.keys():
            if chrom == "meta" or chrom == "pca":
                continue
            group = store[chrom]
            missing_p_l.append(group["missing_rates"].value)
            hz.append(group["counts"].value[:, :3])
            af.append(group["allele_freq"].value)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 20))
    af = np.concatenate(af)
    missing_p_l = np.concatenate(missing_p_l)
    sb.set()
    sb.set_style("white")

    # missing per loci
    sb.distplot(missing_p_l, kde=False, norm_hist=False, ax=axs[0])
    axs[0].set_title("Missing per snp")
    axs[0].set_ylabel('counts')
    axs[0].set_xlabel("Fraction missing")
    axs[0].set_yscale('log')

    # hwe
    hz = np.concatenate(hz)
    sums = np.sum(hz, axis=1)
    ax = sb.regplot(y=hz[:, 0]/sums, x=af, ax=axs[1], fit_reg=False, color='blue',
                    label="Homozygotes ref", scatter_kws={'s': 5, 'alpha': .5})
    ax = sb.regplot(y=hz[:, 2]/sums, x=af, ax=ax, fit_reg=False, color='green',
                    label="Homozygotes alt", scatter_kws={'s': 5, 'alpha': .5})
    ax = sb.regplot(y=hz[:, 1]/sums, x=af, ax=ax, fit_reg=False, color='red',
                    label="Heterozygotes", scatter_kws={'s': 5, 'alpha': .5})
    af_ = np.array(list(range(100)))/100
    ax.plot(af_, af_*af_, label="Expected hom. ref.", linewidth=4, color='black', linestyle='--')
    ax.plot(af_, 2*af_*(1-af_), label="Expected het", color='black', linewidth=4, linestyle='--')
    ax.plot(af_, (1-af_)*(1-af_), label="Expected hom. alt.", color='black', linewidth=4, linestyle='--')
    ax.legend(bbox_to_anchor=(0, 1), loc='best', ncol=1, prop={'size': 15}, markerscale=8)
    plt.tight_layout()

    fig.savefig(outfile, bbox_inches='tight')


def manhattan_plot(plinkFile, outfile):
    pos = []
    pval = []
    chroms = []
    betas = []
    chromosome = []
    with h5py.File(plinkFile, 'r') as store:
        for chrom in store.keys():
            if chrom == "meta"or chrom == "pca":
                continue
            coef = store[f"meta/{chrom}/newton_coef"].value[:, 1]
            betas.append(coef)
            pos.append(store[f"{chrom}/positions"].value[coef[:, 0] != 0])
            pval.append(chi2sf(-2*store[f"meta/{chrom}/newton_ell"].value, 1))
            chroms.append([int(chrom)] * coef.shape[0])
            chromosome.append(chrom)
    fbase = os.path.splitext(outfile)[0]
    df = pd.DataFrame(data={"chrom": np.concatenate(chroms),  "positions": np.concatenate(pos),
                      "pval": np.concatenate(pval)[:, 0], "coef": np.concatenate(betas)[:, 0]})
    df.to_csv(fbase + ".txt", sep='\t', index=False)
    ax, fig = manhattan(pos, pval, chromosome, alpha=.7)
    plt.savefig(outfile)

def manhattan(positions, y, chroms, ax=None, alpha=.7, log_trans=True, **kwargs):
    sb.set(font_scale=1)
    sb.set_style("white")
    sb.despine(left=True)
    colors = sb.color_palette("hls", 8)[5:7]
    _eps = np.min(y[y != 0])
    tot = 0
    for i, position in enumerate(positions):
        positions[i] = position + tot
        tot += position[-1]
        tmp = y[i]
        tmp[tmp == 0] = _eps/10
        if log_trans:
            y[i] = -np.log10(tmp)
    edge_clearance = .01
    fig, ax = plt.subplots(figsize=(12, 6))
    tick_pos = []
    start = 0
    for i, chrom in enumerate(chroms):
        tick_pos.append((positions[i][-1]+start)/2/tot)
        start = positions[i][-1]
        ax.scatter(positions[i]/tot, y[i], c=np.array([colors[i % 2]]), alpha=alpha, edgecolors='none', **kwargs)
        ax.set_xlim(-edge_clearance, 1+edge_clearance)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(labels=chroms)
    ax.set_ylabel(r'$-\log_{10}($P$)$', fontsize=16)
    sb.despine(fig=None, ax=ax, top=True, right=True, left=True, bottom=True, offset=None, trim=False)
    return ax, fig
