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

"""The core functions of MotifScan.
    
    motif_scan(): the interface of MotifScan core
    target_enrichment(): perform the enrichment analysis between sample and random
    target_enrichment_peak2peak(): perform enrichment analysis between two samples
    fc_tarnum_distribution(): fold change plot graph of each motif ranked by peak mvalue
    target_site_distribution(): motif target distribution graph around peak summit

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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import math
import os
import os.path
import ctypes
import MotifScan


class MAT(ctypes.Structure):
    _fields_ = [("n", ctypes.c_int),
                ("a_arr", ctypes.POINTER(ctypes.c_double)),
                ("c_arr", ctypes.POINTER(ctypes.c_double)),
                ("g_arr", ctypes.POINTER(ctypes.c_double)),
                ("t_arr", ctypes.POINTER(ctypes.c_double))]


class MOTIF_RES(ctypes.Structure):
    _fields_ = [("tarnum", ctypes.c_int),
                ("ratio", ctypes.c_double),
                ("tarsite", ctypes.POINTER(ctypes.c_int)),
                ("tarratio", ctypes.POINTER(ctypes.c_double))]


def arr2MAT(arr):
    n = np.shape(arr)[1]
    a = np.ctypeslib.as_ctypes(arr[0])
    c = np.ctypeslib.as_ctypes(arr[1])
    g = np.ctypeslib.as_ctypes(arr[2])
    t = np.ctypeslib.as_ctypes(arr[3])
    return MAT(n, a, c, g, t)


def motif_scan(peak_table, motif_table, background, tmp_output):
    """ The motif scan interface

        Perform motif scanning by motifs

    Args:
        peak_table: pandas dataframe containing sequence matrix
        motif_table: pandas dataframe containing motif basic information
        background: a numpy array containing ratio of A C G T on the genome
        tmp_output: temporary motifscan output file for a certain motif
    """

    score_c = ctypes.CDLL('%s/score_c.so' % os.path.dirname(os.path.realpath(MotifScan.__file__)))
    score_c.motif_scan_core.restype = ctypes.POINTER(MOTIF_RES)
    score_c.motif_scan_core.argtypes = [ctypes.POINTER(MAT),
                                        ctypes.POINTER(MAT),
                                        ctypes.POINTER(ctypes.c_double*4),
                                        ctypes.c_double,
                                        ctypes.c_double]
    score_c.freeMOTIF_RES.argtypes = [ctypes.POINTER(MOTIF_RES)]

    n_motif = len(motif_table)
    background = np.ctypeslib.as_ctypes(background)
    cnt = 0
    for idx, motif_record in motif_table.iterrows():
        motif_name = motif_record['name']
        motif_len = np.shape(motif_record['matrix'])[1]
        score_cutoff = ctypes.c_double(motif_record['score_cutoff'])
        max_score = ctypes.c_double(motif_record['max_score'])
        mmatrix = arr2MAT(motif_record['matrix'])
        sys.stdout.write("Scan ning for %s...%s%s\
                          \r" % (motif_record['name'], cnt*100/n_motif, "%"))
        sys.stdout.flush()
        ratio, tarnum, tarsite, tarratio = [], [], [], []
        for p_idx, peak_record in peak_table.iterrows():
            smatrix = arr2MAT(peak_record['seq_matrix'])
            r = score_c.motif_scan_core(ctypes.byref(smatrix),
                                        ctypes.byref(mmatrix),
                                        ctypes.byref(background),
                                        max_score,
                                        score_cutoff)

            ratio.append(r.contents.ratio)

            ctypes_tarsite = (ctypes.c_int * r.contents.tarnum).from_address(ctypes.addressof(r.contents.tarsite.contents))
            ts = []
            for i in np.arange(r.contents.tarnum):
                ts.append(ctypes_tarsite[i])

            ctypes_tarratio = (ctypes.c_double * r.contents.tarnum).from_address(ctypes.addressof(r.contents.tarratio.contents))
            tr = []
            for i in np.arange(r.contents.tarnum):
                tr.append(ctypes_tarratio[i])

            # deduplicate overlap motif target sites
            ts, tr = deduplicate_target_site(ts, tr, motif_len)
            tarnum.append(len(tr))
            tarsite.append(ts)
            tarratio.append(tr)

            score_c.freeMOTIF_RES(r)
        peak_result = pd.DataFrame({'%s.ratio' % motif_name: ratio,
                                    '%s.tarnum' % motif_name: tarnum,
                                    '%s.tarsite' % motif_name: tarsite,
                                    '%s.tarratio' % motif_name: tarratio})
        peak_result.to_pickle("%s/%s" % (tmp_output, motif_record['id']))
        cnt += 1
    print 'Scanning...Done!            \r'

def deduplicate_target_site(ts, tr, motif_len):
    """
        if the distance between two neignbor target site is smaller than
        the motif length, we viewed them as duplicate target sites and thus filter out 
        the one with the lower ratio.
    """
    if len(ts) >= 2:
        i_pre = float("-Inf")
        j_pre = float("-Inf")
        for (i, j) in zip(ts,tr):
            if i - i_pre < motif_len:
                if j < j_pre:
                   ts.remove(i)
                   tr.remove(j)
                else:
                   ts.remove(i_pre)
                   tr.remove(j_pre)
                   i_pre = i
                   j_pre = j
            else:
                i_pre = i
                j_pre = j
    return ts, tr
                           
def target_enrichment(peak_table, rnd_table, motif_table, gene_based=False):
    """Perform the enrichment analysis between sample and random.

    Args:
        peak_table: pandas dataframe, motifscan result table on sample
        rnd_table: pandas dataframe, motifscan result table on random
        motif_table: pandas dataframe, motif information table

    Returns:
        motif_table: pandas dataframe, table containing both motif information and
                                       fisher exact test statistics
    """
    n_motif = len(motif_table)
    n_peak = len(peak_table)
    n_rand = len(rnd_table)
    n_samp = int(n_rand/n_peak)

    fold_change = np.zeros(n_motif)
    enrich_pvalue = np.zeros(n_motif)
    deplete_pvalue = np.zeros(n_motif)
    oddsratio = np.ones(n_motif)
    pvalue_corrected = np.ones(n_motif)

    peak_tarnum = np.zeros(n_motif)
    rand_tarnum = np.zeros(n_motif)
    peak_tarnum_table = peak_table[
        pd.Index([i for i in peak_table.columns if re.search(r'\.tarnum', i)])]
    rnd_tarnum_table = rnd_table[
        pd.Index([i for i in rnd_table.columns if re.search(r'\.tarnum', i)])]

    for mti, motif_name in zip(range(n_motif), motif_table['name']):
        if gene_based:
            targeted_peak_idx = peak_table.groupby(['target_gene'])['%s.tarnum' % motif_name].transform(max) == peak_table['%s.tarnum' % motif_name]
            targeted_peak_idx = peak_table.loc[targeted_peak_idx,['target_gene','%s.tarnum' % motif_name]].drop_duplicates().index
            targeted_rand_idx = []
            for i in targeted_peak_idx:
                for j in np.arange(n_samp):
                    targeted_rand_idx.append(i*n_samp+j)

            targeted_rand_idx = pd.Index(targeted_rand_idx)
            peak_tarnum[mti] = len(
                [i for i in peak_tarnum_table.loc[targeted_peak_idx, '%s.tarnum' % motif_name] if i > 0])
            rand_tarnum[mti] = len(
                [i for i in rnd_tarnum_table.loc[targeted_rand_idx, '%s.tarnum' % motif_name] if i > 0])
        else: 
            peak_tarnum[mti] = len(
                [i for i in peak_tarnum_table['%s.tarnum' % motif_name] if i > 0])
            rand_tarnum[mti] = len(
                [i for i in rnd_tarnum_table['%s.tarnum' % motif_name] if i > 0])

        if peak_tarnum[mti] != 0 and rand_tarnum[mti] != 0:
            fold_change[mti] = float(peak_tarnum[mti] * n_rand) / (
                rand_tarnum[mti] * n_peak)
        else:
            fold_change[mti] = 'NaN'
        table = [[peak_tarnum[mti], n_peak - peak_tarnum[mti]],
                 [rand_tarnum[mti], n_rand - rand_tarnum[mti]]]
        oddsratio[mti], enrich_pvalue[mti] = stats.fisher_exact(table, 'greater')
        oddsratio[mti], deplete_pvalue[mti] = stats.fisher_exact(table, 'less')
        pvalue_corrected[mti] = min(min(deplete_pvalue[mti],
                                        enrich_pvalue[mti]) * n_motif, 1)
    motif_table['target_number'] = peak_tarnum
    motif_table['rnd_target_number'] = rand_tarnum
    motif_table['fold_change'] = fold_change
    motif_table['enrich_pvalue'] = enrich_pvalue
    motif_table['deplete_pvalue'] = deplete_pvalue
    motif_table['pvalue_corrected'] = pvalue_corrected

    return motif_table


def target_enrichment_peak2peak(peak1_table, peak2_table, motif_table):
    """Perform the enrichment analysis on two samples
    Args:
        peak_table: pandas dataframe, motifscan result table on sample1
        rnd_table: pandas dataframe, motifscan result table on sample2
        motif_table: pandas dataframe, motif information table

    Returns:
        motif_table: pandas dataframe, table containing both motif information and
                                       fisher exact test statistics
    """

    n_motif = len(motif_table)
    n_peak1 = len(peak1_table)
    n_peak2 = len(peak2_table)

    fold_change = np.zeros(n_motif)
    enrich_pvalue = np.zeros(n_motif)
    deplete_pvalue = np.zeros(n_motif)
    oddsratio = np.ones(n_motif)
    pvalue_corrected = np.ones(n_motif)

    peak1_tarnum = np.zeros(n_motif)
    peak2_tarnum = np.zeros(n_motif)
    # print pd.Index([i for i in peak1_table.columns if re.search(r'\.tarnum',i)])
    peak1_tarnum_table = peak1_table[
        pd.Index([i for i in peak1_table.columns if re.search(r'\.tarnum', i)])]
    peak2_tarnum_table = peak2_table[
        pd.Index([i for i in peak2_table.columns if re.search(r'\.tarnum', i)])]

    for mti, motif_name in zip(range(n_motif), motif_table['name']):
        peak1_tarnum[mti] = len(
            [i for i in peak1_tarnum_table['%s.tarnum' % motif_name] if i > 0])
        peak2_tarnum[mti] = len(
            [i for i in peak2_tarnum_table['%s.tarnum' % motif_name] if i > 0])
        if peak1_tarnum[mti] != 0 and peak2_tarnum[mti] != 0:
            fold_change[mti] = float(peak1_tarnum[mti] * n_peak2) / (
                peak2_tarnum[mti] * n_peak1)
        else:
            fold_change[mti] = 'NaN'
        table = [[peak1_tarnum[mti], n_peak1 - peak1_tarnum[mti]],
                 [peak2_tarnum[mti], n_peak2 - peak2_tarnum[mti]]]
        oddsratio[mti], enrich_pvalue[mti] = stats.fisher_exact(table, 'greater')
        oddsratio[mti], deplete_pvalue[mti] = stats.fisher_exact(table, 'less')
        pvalue_corrected[mti] = min(min(deplete_pvalue[mti],
                                        enrich_pvalue[mti]) * n_motif, 1)
    motif_table['peak1_target_number'] = peak1_tarnum
    motif_table['peak2_target_number'] = peak2_tarnum
    motif_table['fold_change'] = fold_change
    motif_table['enrich_pvalue'] = enrich_pvalue
    motif_table['deplete_pvalue'] = deplete_pvalue
    motif_table['oddsratio'] = oddsratio
    motif_table['pvalue_corrected'] = pvalue_corrected
    motif_table.sort('enrich_pvalue', inplace=True)
    return motif_table


#################################################################################
#  plot functions
#################################################################################
def target_site_distribution(peak_table, motif_table, plot_out_dir, bin_size=5, region_radius=500):
    win_size = region_radius * 2
    half_win_size = win_size / 2

    bin_center = np.arange(-half_win_size, win_size - half_win_size + 1, bin_size)
    bin_edge = bin_center - round(bin_size / 2)
    bin_edge = np.append(bin_edge, bin_center[-1] + round(bin_size / 2))
    bin_edge = map(int, bin_edge)
    n_motif = len(motif_table)

    cnt = 0
    for midx, motif_record in motif_table.iterrows():
        motif_len = np.shape(motif_record['matrix'])[1]
        motif_name = motif_record['name']
        motif_tarsite = []

        for idx, i in enumerate(peak_table['%s.tarsite' % motif_name]):
            if len(i) > 0:
                for j in i:
                    motif_tarsite.append(j + motif_len/2 - half_win_size)

        motif_tarsite_freq = np.histogram(motif_tarsite, bin_edge)
        motif_marker = motif_tarsite_freq[0]
        motif_marker[0] = motif_marker[0] * 2
        motif_marker[-1] = motif_marker[-1] * 2

        motif_marker = motif_marker / float(sum(motif_marker))
        motif_marker = smooth(motif_marker, 20)

        plt.cla()
        #plt.plot(bin_edge[:-1], motif_marker, lw = 2, color='#4169E1', label=motif_name)
        plt.bar(bin_edge[:-1], motif_marker, width=bin_size-1.5,
                color='#4169E1', linewidth=0, label=motif_name)
        plt.legend()
        ax = plt.gca()
        ax.set_xlabel('Distance to Peak Summit', weight='bold')
        ax.set_ylabel('Fraction', weight='bold')
        ax.set_xlim([min(bin_edge), max(bin_edge)])
        ax.set_ylim([0, 1.5 * max(motif_marker)])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        ax.tick_params(axis='x', which='both', top='off')
        ax.tick_params(axis='y', which='both', right='off')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(2)
        plt.savefig('%s/%s_target_site.png' % (plot_out_dir, re.sub(r'[\:\-\.]{1,2}', '_', motif_name)), dpi=600)
        plt.close()

        # sys.stdout.write("\r")
        sys.stdout.write("Target sites distribution plotting for %s...%s%s            \r" % (motif_name, (cnt)*100/n_motif, "%"))
        sys.stdout.flush()
        cnt += 1
    print "Target sites distribution plotting...Done!           \r"


def tarnum_and_tarsite_distribution(peak_table, rand_table, motif_table, plot_out_dir, bin_size=10, region_radius=500):
    # fold change plot parameter
    motif_name_list = extract_motif_name_from_peak_result_table(peak_table)
    npeak = len(peak_table)
    nrand = len(rand_table)
    nsamp = round(nrand / npeak)
    bin = 1000
    half_bin = bin/2
    peak_table.sort('value', ascending=False, inplace=True)
    n_motif = len(motif_name_list)

    # target site paramater
    win_size = region_radius * 2
    half_win_size = win_size / 2

    bin_center = np.arange(-half_win_size, win_size - half_win_size + 1, bin_size)
    bin_edge = bin_center - round(bin_size / 2)
    bin_edge = np.append(bin_edge, bin_center[-1] + round(bin_size / 2))
    bin_edge = map(int, bin_edge)
    n_motif = len(motif_table)

    cnt = 0
    for midx, motif_record in motif_table.iterrows():
        cnt += 1
        motif_len = np.shape(motif_record['matrix'])[1]
        motif_name = motif_record['name']
        # fold change data preparation
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

        tarnum_fc_smooth = pd.Series(tarnum_fc_smooth)
        tarnum_smooth = pd.Series(tarnum_smooth)

        # target site data preparation
        motif_tarsite = []
        for idx, i in enumerate(peak_table['%s.tarsite' % motif_name]):
            if len(i) > 0:
                for j in i:
                    motif_tarsite.append(j + motif_len/2 - half_win_size)
        motif_tarsite_freq = np.histogram(motif_tarsite, bin_edge)
        motif_marker = motif_tarsite_freq[0]
        motif_marker[0] = motif_marker[0] * 2
        motif_marker[-1] = motif_marker[-1] * 2

        motif_marker = motif_marker / float(sum(motif_marker))
        motif_marker = smooth(motif_marker, 20)

        # plot
        plt.cla()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        # fc
        fig.suptitle(motif_name, weight='black')
        ax1.bar(np.arange(npeak), tarnum_fc_smooth, 1,
                color='#4169E1', lw=0, label='fc of %s' % motif_name)
        xtick_step = _get_highest_digit(npeak/6)*10**int(math.log10(npeak/6))
        ax1.set_xticks(np.arange(0, npeak, xtick_step))
        ax1.set_xlim(xmin=1, xmax=npeak)
        ax1.set_ylim(ymin=0, ymax=1.3*max(tarnum_fc_smooth))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(axis='x', which='both', top='off')
        ax1.tick_params(axis='y', which='both', right='off')
        ax1.set_ylabel('Fold Change', weight='bold')
        ax1.set_xlabel('Peak Rank (Sorted by Value in Descending Order)', weight='bold')
        for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontweight('bold')

        # targetsite
        ax2.bar(bin_edge[:-1], motif_marker, width=bin_size-1.5,
                color='#4169E1',  linewidth=0,  label=motif_name)
        ax2.set_xlabel('Distance to Peak Summit', weight='bold')
        ax2.set_ylabel('Fraction', weight='bold')
        ax2.set_xlim([min(bin_edge),  max(bin_edge)])
        ax2.set_ylim([0,  1.5 * max(motif_marker)])
        for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontweight('bold')
        ax2.tick_params(axis='x', which='both', top='off')
        ax2.tick_params(axis='y', which='both', right='off')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax2.spines[axis].set_linewidth(2)
        fig.savefig('%s/%s_%s_tarsite_fc_dist.png' % (plot_out_dir, cnt, re.sub(r'[\:\-\.]{1,2}', '_', motif_name)), dpi=600)
        plt.close()

        sys.stdout.write("Fold change distribution plotting for %s...%s%s           \r" % (motif_name, (cnt)*100/n_motif, "%"))
        sys.stdout.flush()
    print 'Fold change distribution plotting...Done!         \r'


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
    if not window in ['flat',  'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett',  blackman'"

    s = np.r_[2*x[0] - x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
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
