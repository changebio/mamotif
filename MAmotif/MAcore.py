#
# Copyright @ 2014, 2015 Jiawei Wang <jerryeah@gmail.com>, Zhen Shao <shao@enders.tch.harvard.edu>
#
# Licensed under the GPL License, Version 3.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#   http://www.gnu.org/copyleft/gpl.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 

"""The core functions of MAmotif.
"""

__authors__ = ['"Jiawei Wang" <jerryeah@gmail.com>']


# -----------------------------------------
#  modules
# -----------------------------------------
import numpy as np
import pandas as pd
from scipy import stats
import sys

import matplotlib
import rpy2.robjects as robjects
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import copy
import re
import math


#################################################################################
#  plot functions
#################################################################################
def fc_plot(peak_table, rand_table, plot_out_dir):
    """
    """
    motif_name_list = extract_motif_name_from_peak_result_table(peak_table)
    npeak = len(peak_table)
    nrand = len(rand_table)
    nsamp = round(nrand / npeak)
    bin = 1500
    half_bin = bin/2
    peak_table.sort('value', ascending=False, inplace=True)
    peak_table.reset_index(inplace=True)
    n_motif = len(motif_name_list)
    for idx, motif_name in enumerate(motif_name_list):
        tarnum_fc_smooth = np.zeros(npeak)
        tarnum_smooth = np.zeros(npeak)
        for pi in np.arange(npeak):
            peak_start_idx = max(0, pi - half_bin)
            peak_end_idx = min(npeak, pi + half_bin)
            peak_tarnum = peak_table['%s.tarnum' %
                                     motif_name].iloc[peak_start_idx:peak_end_idx]

            peak_tarnum_smooth = float(
                len(peak_tarnum[peak_tarnum > 0])) / len(peak_tarnum)
            rnd_idx = peak_table.iloc[peak_start_idx:peak_end_idx].index*int(nsamp)
            rand_tarnum = rand_table['%s.tarnum' %
                                     motif_name].ix[rnd_idx]
            rand_tarnum_smooth = float(
                len(rand_tarnum[rand_tarnum > 0])) / len(rand_tarnum)

            tarnum_smooth[pi] = peak_tarnum_smooth
            if rand_tarnum_smooth == 0:
                tarnum_fc_smooth[pi] = peak_tarnum_smooth
            else:
                tarnum_fc_smooth[pi] = peak_tarnum_smooth / rand_tarnum_smooth

            #raw_input()
        tarnum_fc_smooth = smooth(tarnum_fc_smooth, 50)
        tarnum_fc_smooth = pd.Series(tarnum_fc_smooth)
        tarnum_smooth = pd.Series(tarnum_smooth)

        value = peak_table['value'].tolist()
        value = smooth(value, 50)
        value = pd.Series(value)
        #plt.cla()
        fig = plt.figure(figsize=(12, 6))
        plt.plot(tarnum_fc_smooth,
                 color='#4169E1',
                 lw=2,
                 label='fc of %s' % motif_name)
        plt.hlines(1, 1, npeak, linestyles='dashed', colors='#4169E1')
        ax1 = plt.gca()
        xtick_step = _get_highest_digit(npeak/6)*10**int(math.log10(npeak/6))
        ax1.set_xticks(np.arange(0, npeak, xtick_step))
        ax1.set_xlim(xmin=1, xmax=npeak)
        ax1.set_ylim(ymin=0, ymax=1.3*max(tarnum_fc_smooth))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(axis='x', which='both', top='off')
        ax1.tick_params(axis='y', which='both', right='off')
        ax1.set_ylabel('Fold Change', weight='bold', color='#4169E1')
        ax1.set_xlabel('Peak Rank (Sorted by Value in Descending Order)', weight='bold')
        for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontweight('bold')

        plt.title(motif_name, weight='bold')
        ax2 = ax1.twinx()
        ax2.plot(value, color='#E25668', lw=2)
        ax2.set_ylabel('Mvalue', weight='bold', color='#E25668')
        ax2.set_ylim(ymin=-10, ymax=10)
        ax2.set_xlim(xmin=1, xmax=npeak)
        plt.hlines(-1, 1, npeak, linestyles='dashed', colors='#E25668')
        plt.hlines(1, 1, npeak, linestyles='dashed', colors='#E25668')

        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontweight('bold')
        plt.vlines(_find_nearest_number(1, peak_table['value'])+1, -10, 10, linestyles='dashed', lw=2)
        plt.vlines(_find_nearest_number(-1, peak_table['value'])+1, -10, 10, linestyles='dashed', lw=2)
        fig.savefig('%s/%s_fc_dist.png' % (plot_out_dir,  re.sub(r'[\:\-\.]{1, 2}', '_', motif_name)), dpi=600)
        plt.close()

        sys.stdout.write("Fold change distribution on all peaks for %s...%s%s           \r" % (motif_name, (idx+1)*100/n_motif, "%"))
        sys.stdout.flush()
    print 'Fold change distribution on all peaks...Done!         \r'


def _find_nearest_number(n, sl):
    d = abs(n-sl[0])
    m = 0
    for idx, i in sl.iteritems():
        if abs(n-i) < d:
            d = abs(n-i)
            m = idx
    return m


def h_cluster(peak_result, peak_test_result, output,
              fc_cutoff=1.3, pvalue_cutoff=0.01, pcc_cutoff=0.75, negative=False):
    mvalues = np.array([peak_result['value']])
    if negative:
        mvalues = -mvalues
    #otif_name_list = extract_motif_name_from_peak_result_table(peak_result)
    if len(peak_test_result) > 25:
        motif_name_list = peak_test_result.iloc[:25]["Motif"].tolist()
    else:
        motif_name_list = peak_test_result["Motif"].tolist()
    tarnum_list = []
    for idx, name in enumerate(motif_name_list):
        tarnum_list.append(peak_result['%s.tarnum' % name])
    tarnum_mat = np.array(tarnum_list)
    tarnum_mat[tarnum_mat > 0] = 1  # convert to 0/1matrix

    cluster_mat = tarnum_mat
    if np.shape(cluster_mat)[0] < 5:
        print 'Two few motif candidate for clustering.'
        return 0
    cluster_mat = np.concatenate((mvalues, cluster_mat), axis=0)
    # prepare label info
    motif_name_list = np.array(motif_name_list)
    motif_label = motif_name_list
    label = np.insert(motif_label, 0, 'Mvalue').tolist()
    # PW distance and linkage
    D = pdist(cluster_mat, 'correlation')  # correlation
    L = linkage(D, method='average')  # average
    # plot dendrogram
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_axes([0.1, 0.71, 0.75, 0.28])
    dd = dendrogram(L, orientation='top', labels=label)
    ax1.set_xticks([])
    ax1.set_yticks([])
    label_dendro = dd['ivl']

    cluster_matrix_dendro = copy.deepcopy(cluster_mat)
    for i in np.arange(len(label_dendro)):
        idx = label.index(label_dendro[i])
        cluster_matrix_dendro[i] = cluster_mat[idx]

    D_corr = pdist(cluster_matrix_dendro, 'correlation')
    # plot correlation matrix

    #register the colormap
    startcolor = '#3D3888'
    midcolor = '#ffffff'
    endcolor = '#B32E27'
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list('own', [startcolor, midcolor, endcolor])

    axmatrix = fig.add_axes([0.1, 0.01, 0.75, 0.75])
    mat = 1 - squareform(D_corr)
    im = axmatrix.matshow(mat, aspect='auto', cmap=cmap1)
    im.set_clim(vmin=-1, vmax=1)
    axmatrix.set_yticks(range(0, len(label_dendro), 1))
    axmatrix.set_yticklabels(label_dendro)
    axmatrix.set_xticklabels([])

    cbax = fig.add_axes([0.85, 0.1, 0.1, 0.1])
    cbax.set_visible(False)
    fig.colorbar(im, ax=cbax, aspect='auto', ticks=[-1, 0, 1])
    fig.savefig('%s.png' % output, dpi=600)


def peak_set_enrichment_analysis(peak_table, plot_output_dir, permutation=1000):
    motif_names = extract_motif_name_from_peak_result_table(peak_table)
    n_motif = len(motif_names)
    for idx, motif_name in enumerate(motif_names):
    # walk through the sorted list
        target_yes, target_no = split_manorm_peaks_motif(
            peak_table, motif_name)
        peak_set = target_yes.index
        L = peak_table['value'].copy().squeeze()
        L.sort(ascending=False)

        ES, hit_cum = calculate_ES(L, peak_set)
        if ES and hit_cum:
            ES_plot(hit_cum, L, peak_set, motif_name, '%s/%s_ES.png' % (plot_output_dir, motif_name))
        sys.stdout.write("Peak set enrichment analysis... %s...%s%s           \r" % (motif_name, (idx+1)*100/n_motif, "%"))
        sys.stdout.flush()
    print 'Peak set enrichment analysis...Done!         \r'
    # ES_null = []
    # for i in np.arange(permutation):
    #     L.index = np.random.permutation(L.index) # permutation method 2
    #     ES_null.append(calculate_ES(L, peak_set))
    # ES_null = pd.Series(ES_null)
    # NES = ES/np.mean(ES_null)
    # pvalue = ES_null.loc[ES_null>=ES].count()/float(permutation)
    # return ES


def calculate_ES(L, N):
    """
    @para L: ranked gene list, datatype: pandas Series
    @para N: gene set, DataType: pandas index
    """
    # assert unique
    if len(N.unique()) != len(N):
        print "Elements in N should be unique! skipping..."
        return None, None
    elif len(L[N].dropna()) != len(N) or len(L) <= len(N) or len(N) == 0:
        print "N should be a subset of L! skipping..."
        return None, None

    total_hit_sum = pd.Series(map(abs, L[N])).sum()
    total_miss = len(L) - len(N)
    hit_cum = []
    iloc = 0
    for idx, value in L.iteritems():
        iloc += 1
        hit_idx = L[:iloc].index.intersection(N)
        hit_fraction = float(pd.Series(map(abs, L[hit_idx])).sum())/total_hit_sum  # which is weighted
        miss_fraction = float(iloc-len(hit_idx))/total_miss
        hit_cum.append(hit_fraction-miss_fraction)
    return np.max(map(abs, hit_cum)), hit_cum


def ES_plot(hit_cum, L, peak_set, motif_name, plot_output):
    max_site = np.argmax(map(abs, hit_cum))
    max_value = hit_cum[max_site]
    fig = plt.figure(figsize=(12, 6))
    idx = map(pd.Index.get_loc,  [L.index]*len(peak_set), peak_set.tolist())
    L = smooth(L.tolist(), 20)

    plt.bar(np.arange(len(L)), L, 1,
            color='#C6C6C6', lw=0)
    ax1 = plt.gca()

    ax1.bar(idx, [2]*len(peak_set), 1, [-10]*len(peak_set), color='#000000', lw=0)
    ax1.set_ylabel("mvalue", color='#BFBFBF')
    ax1.vlines(max_site, -10, 10, linestyles='dashed', lw=2)
    ax1.set_xlabel("Peak ranked by Mvalue (Descending)")
    ax1.set_xlim([0, len(L)])
    ax1.set_ylim(ymin=-10, ymax=10)

    ax2 = ax1.twinx()
    hit_cum = smooth(hit_cum, 20)
    ax2.plot(hit_cum, color="#E25668", lw=2)
    ax2.set_xlim([0, len(L)])
    ax2.set_ylim(ymin=-abs(max_value)*1.2, ymax=abs(max_value)*1.2)
    ax2.set_ylabel("Running Enrichment Score", color="#E25668")
    ax2.tick_params(axis='x', which='both', top='off')
    ax2.tick_params(axis='y', which='both', right='off')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.annotate('ES=%.2f' % max_value, xy=(max_site, max_value), xytext=(max_site+5, max_value*1.05))
    plt.title(motif_name)
    fig.savefig(plot_output, dpi=600)
    plt.close()


def peak_result_test(peak_table, rnd_table=None, negative=False):
    '''
    doing test for each motif
    '''
    npeak = len(peak_table)

    motif_names = extract_motif_name_from_peak_result_table(peak_table)
    n_motif = len(motif_names)
    yes_mean, no_mean, yes_std, no_std = [], [], [], []
    t_statistic, t_pvalue_right, t_pvalue_corrected = [], [], []
    r_statistic, r_pvalue_right, r_pvalue_corrected = [], [], []
    comment = []
    tarnum, no_tarnum = [], []

    cnt = 0
    # if not os.path.exists(plot_output_dir):
    #     os.mkdir(plot_output_dir)

    for motif_name in motif_names:
        sys.stdout.write("Statistical test for %s...%s%s           \r" % (motif_name, (cnt)*100/n_motif, "%"))
        sys.stdout.flush()
        target_yes, target_no = split_manorm_peaks_motif(
            peak_table, motif_name)
        yes_mean.append(np.mean(target_yes))
        no_mean.append(np.mean(target_no))
        yes_std.append(np.std(target_yes))
        no_std.append(np.std(target_no))
        tarnum.append(len(target_yes))
        no_tarnum.append(len(target_no))
        # mean and std
        if negative:
            target_yes = -target_yes
            target_no = -target_no
        target_yes = target_yes.tolist()
        target_no = target_no.tolist()
        # motif validation I: motifs that have too few or too much targets will be excluded.
        if len(target_yes) > npeak*0.5:
            comment.append('Too much motif discoveries on peaks;')
        elif len(target_yes) < npeak*0.0005 and len(target_yes) > 0:
            comment.append('Too few motif discoveries on peaks;')
        else:
            comment.append('')
        # motif validation II: motifs that not enriched on total peaks will be excluded.
        if rnd_table:
            nrand = len(rnd_table)
            bin = 1000
            half_bin = bin/2
            nsamp = round(nrand / npeak)
            tarnum_fc_smooth = np.zeros(npeak)
            for pi in np.arange(npeak):
                peak_start_idx = max(0, pi - half_bin)
                peak_end_idx = min(npeak, pi + half_bin)
                peak_tarnum = peak_table['%s.tarnum' %
                                         motif_name].iloc[peak_start_idx:peak_end_idx]

                peak_tarnum_smooth = len(peak_tarnum[peak_tarnum > 0])
                rnd_idx = peak_table.iloc[peak_start_idx:peak_end_idx].index*int(nsamp)
                rand_tarnum = rnd_table['%s.tarnum' % motif_name].ix[rnd_idx]
                rand_tarnum_smooth = len(rand_tarnum[rand_tarnum > 0])

                if rand_tarnum_smooth == 0:
                    tarnum_fc_smooth[pi] = peak_tarnum_smooth
                else:
                    tarnum_fc_smooth[pi] = float(peak_tarnum_smooth) / rand_tarnum_smooth
            tarnum_fc_smooth = pd.Series(tarnum_fc_smooth)
            if tarnum_fc_smooth.quantile(0.95) < 1.3 and len(target_yes) > 0:
                comment[-1] += 'Not enriched across peaks; '

        # statistical test
        if len(target_yes) == 0 or len(target_no) == 0:
            t = l = r = 'NaN'
        else:
            t, l, r = t_test_two_sample(target_yes, target_no)
        t_pvalue_right.append(r)
        t_statistic.append(t)
        if len(target_yes) == 0 or len(target_no) == 0:
            t = l = r = 'NaN'
        else:
            t, l, r = ranksum_test_two_sample(target_yes, target_no)
        r_pvalue_right.append(r)
        r_statistic.append(t)
        cnt += 1
    t_pvalue_corrected = _benjamini_r(t_pvalue_right)
    r_pvalue_corrected = _benjamini_r(r_pvalue_right)
    maximal_pvalue = map(max, t_pvalue_corrected, r_pvalue_corrected)
    test_result_tb = pd.DataFrame({
        'Motif': motif_names,
        'Target Number': tarnum,
        'Mean of Target Mvalue': yes_mean,
        'SD of Target Mvalue': yes_std,
        'Non-Target Number': no_tarnum,
        'Mean of Non-Target Mvalue': no_mean,
        'SD of Non-Target Mvalue': no_std,
        'T Statistic': t_statistic,
        'T Test Pvalue (Right Tailed)': t_pvalue_right,
        'Corrected T Test Pvalue (Benjamini)': t_pvalue_corrected,
        'Z Statistic': r_statistic,
        'Ranksum Test Pvalue (Right Tailed)': r_pvalue_right,
        'Corrected Ranksum Test Pvalue (Benjamini)': r_pvalue_corrected,
        'Maximal Pvalue': maximal_pvalue,
        'Comment': comment
        })
    test_result_tb.sort("Maximal Pvalue", inplace=True)
    return test_result_tb

def mark_both_enriched_motifs(test_result_A, test_result_B, pvalue_cutoff=.05):
    t_test_enrich_motifs_A = test_result_A.loc[test_result_A['T Test Pvalue (Right Tailed)']
                                               < pvalue_cutoff]['Motif'].tolist()
    t_test_enrich_motifs_B = test_result_B.loc[test_result_B['T Test Pvalue (Right Tailed)']
                                               < pvalue_cutoff]['Motif'].tolist()
    t_test_both_enriched = set(t_test_enrich_motifs_A+t_test_enrich_motifs_B)

    ranksum_test_enrich_motifs_A = test_result_A.loc[test_result_A['Ranksum Test Pvalue (Right Tailed)']
                                                     < pvalue_cutoff]['Motif'].tolist()
    ranksum_test_enrich_motifs_B = test_result_B.loc[test_result_B['Ranksum Test Pvalue (Right Tailed)']
                                                     < pvalue_cutoff]['Motif'].tolist()
    ranksum_test_both_enriched = set(ranksum_test_enrich_motifs_A+ranksum_test_enrich_motifs_B)

    both_enriched = t_test_both_enriched.union(ranksum_test_both_enriched)
    for motif in both_enriched:
        test_result_A.loc[test_result_A['Motif'] ==
                          motif]['Comment'] += 'Enriched on both samples; '
    return test_result_A, test_result_B


def split_manorm_peaks_motif(peak_result_table, motif_name):
    '''
    split manorm peaks into two sets, have or not have motif target
    '''
    target_yes = peak_result_table.ix[
        peak_result_table['%s.tarnum' % motif_name] > 0, 'value']
    target_no = peak_result_table.ix[
        peak_result_table['%s.tarnum' % motif_name] == 0, 'value']
    return target_yes, target_no


def t_test_two_sample(sample_a, sample_b):
    '''
    doing student's t test on two samples, sample_a and sample_b should be array_like
    '''
    t_statistic, two_tailed_pvalue = stats.ttest_ind(sample_a, sample_b)
    t_statistic = t_statistic.tolist()
    if t_statistic < 0:
        pvalue_left = two_tailed_pvalue / 2
        pvalue_right = 1 - two_tailed_pvalue / 2
    else:
        pvalue_left = 1 - two_tailed_pvalue / 2
        pvalue_right = two_tailed_pvalue / 2
    return t_statistic, pvalue_left, pvalue_right


def _benjamini_r(pvalues):

    """
    compute benjenmini correction value by calling R function
    @param pvalue: a list of pvalues
    @return: a list of corrected p values
    """
    r = robjects.r
    return list(r['p.adjust'](robjects.FloatVector(pvalues), 'BY'))


def ranksum_test_two_sample(sample_a, sample_b):
    '''
    doing ranksum test on two samples, sample_a and sample_b should be array_like
    '''
    z_statistic, two_tailed_pvalue = stats.ranksums(sample_a, sample_b)
    z_statistic = z_statistic.tolist()
    if z_statistic < 0:
        pvalue_left = two_tailed_pvalue / 2
        pvalue_right = 1 - two_tailed_pvalue / 2
    else:
        pvalue_left = 1 - two_tailed_pvalue / 2
        pvalue_right = two_tailed_pvalue / 2
    return z_statistic, pvalue_left, pvalue_right


def extract_motif_name_from_peak_result_table(peak_result_table):
    '''
    extract the motif name from the peak_result_table
    '''
    motif_name = []
    for col_name in peak_result_table.columns:
        if re.search(r'.*\.tarnum', col_name):
            motif_name.append(col_name[:-7])
    return motif_name


def smooth(x, window_len=5, window='hanning'):
    if type(x) == type([]):
        x = np.array(x)
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='same')
    # return the smoothed signal,  chopping off the ends so that it has the previous size.
    return y[window_len-1:-window_len+1]


def _get_highest_digit(n):
    if n / 100000000:
        n /= 100000000
    if n / 10000:
        n /= 10000
    if n / 100:
        n /= 100
    if n / 10:
        n /= 10
    return n
