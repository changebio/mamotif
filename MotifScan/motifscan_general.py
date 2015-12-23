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

"""
The core functions of MotifScan.
"""


__authors__ = ['"Jiawei Wang" <jerryeah@gmail.com>']
# -----------------------------
#  modules
# -----------------------------
import os
import os.path
import pandas as pd
import numpy as np
import sys
import MotifScan.core as core
import MotifScan.peak as peak
import MotifScan.motif as motif
import shutil
import re

# ----------------------------------------------------------------------------------------------------------------
# main functions
# ----------------------------------------------------------------------------------------------------------------


def motifscan_general(peak_path, peak_format, motif_path,
                      motif_list_path, genome_name, gene_path,
                      random_times, peak_length, output_dir,
                      region, up, down,
                      is_enrichment, extract_target_site):
    # --------------------------------------------------------------------------
    # Preparing output directory
    # -------------------------------------------------------------------------
    # ---- output directories --------------
    plot_out_dir = output_dir + '/' + 'plot'
    enrichment_csv = output_dir + '/motif_enrichment.csv'
    peak_motif_tarnum = output_dir + '/peak_motif_tarnum.csv'
    peak_motif_score = output_dir + '/peak_motif_score.csv'
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(plot_out_dir):
        os.mkdir(plot_out_dir)
    if extract_target_site:
        motif_tarsite_out_dir = output_dir + '/' + 'motif_target_sites'
        if not os.path.exists(motif_tarsite_out_dir):
            os.mkdir(motif_tarsite_out_dir)
    tmp_dir = output_dir+'/.tmp'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # -------------------------------------------------------------------------
    #  read basic information: genome, motif and gene annotation
    # -------------------------------------------------------------------------
    user_home = "/picb/rsgeno/huangyin"

    genome_db_path = '%s/.MotifScan/genome/%s' % (user_home, genome_name)
    if gene_path is None and os.path.exists('%s/.MotifScan/gene/%s' % (user_home, genome_name)):
        gene_path = '%s/.MotifScan/gene/%s/refSeq.txt' % (user_home, genome_name)

    print '#########################################################'
    print '                 Data Preparation'
    print '#########################################################'
    # -------- loading gene annotation ----------------
    if gene_path:
        gene_table = peak.load_ref_gene(gene_path)
    else:
        gene_table = None
        if region != "genome":
            print 'Error! You can assign region option only when gene annotation is provided.'
            exit()

    # ---------- loading genome ----------------------------
    print 'Loading genome information...',
    background = pd.read_pickle('%s/background' % genome_db_path)['background']
    chromosome_size = pd.read_pickle('%s/chromosome_size' % genome_db_path)
    print 'Done!'

    # ---------- loading motif ------------------
    print 'Loading motif...',
    motif_table = motif.load_motif(motif_list_path, motif_path)

    print 'Done! %d motifs are processed!' % len(motif_table['name'].unique())

    if len(motif_table) == 0:
        print 'Warning: There is no detected motifs in the motif list!'
        exit()
    # -------------------------------------------------------
    #  loading peaks and generate random peaks if necessary
    # -----------------------------------------------------
    print 'Generating peak sequence...',
    peak_table = peak.load_peak(peak_path, genome_db_path, peak_length, peak_format)
    print 'Done! %d peaks are processed!' % len(peak_table)
    

    if region in ['promoter','distal','whole']:
        if not isinstance(gene_table,pd.DataFrame):
            print "Gene annotation file is required."
            exit()
        print 'Split promoter/distal regions...',
        sys.stdout.flush()
        promoter_start = []
        promoter_end = []
        for i, gene_i in gene_table.iterrows():
            if gene_i['strand'] == '+':
                promoter_start.append(gene_i['TSS']-up)
                promoter_end.append(gene_i['TSS']+down)
            else:
                promoter_start.append(gene_i['TSS']-up)
                promoter_end.append(gene_i['TSS']+down)
        gene_table['promoter_start'] = promoter_start
        gene_table['promoter_end'] = promoter_end
        is_promoter = pd.Series(np.zeros(len(peak_table),dtype=bool))
        for chri in peak_table.chr.unique():
            peak_table_chr = peak_table[peak_table.chr==chri]
            gene_table_chr = gene_table[gene_table.chr==chri]
            p_start_chr=np.tile(np.array([peak_table_chr.start]).transpose(),(1,len(gene_table_chr)))
            g_start_chr = np.tile(gene_table_chr.promoter_start,(len(peak_table_chr),1))
            p_end_chr = np.tile(np.array([peak_table_chr.end]).transpose(),(1,len(gene_table_chr)))
            g_end_chr = np.tile(gene_table_chr.promoter_end,(len(peak_table_chr),1))
            start_chr = p_start_chr<g_end_chr
            end_chr = p_end_chr > g_start_chr
            index = start_chr==end_chr
            is_promoter[peak_table.chr==chri] = index.sum(axis=1)!=0
        peak_table['is_promoter']=is_promoter

        if region == 'promoter':
            peak_table = peak_table.ix[peak_table['is_promoter']]
            peak_table.reset_index(inplace=True)
            print '%s promoter peaks extracted!' % len(peak_table)
        elif region == 'distal':
            peak_table = peak_table.ix[~peak_table['is_promoter']]
            peak_table.reset_index(inplace=True)
            print '%s distal peaks extracted!' % len(peak_table)
        elif region == 'whole':
            print '%s whole peaks extracted!' % len(peak_table)
    elif region in ['gene1', 'gene2']:
        if not isinstance(gene_table,pd.DataFrame):
            print "Gene annotation file is required."
            exit()
        print 'Find peak regions target gene...'
        peak_table = peak.extract_target_gene(peak_table,gene_table)
        peak_table = peak_table.loc[peak_table['target_gene'] != 'No Target']
        if region == 'gene1':
            # pick up the peak that is nearest to the target gene
            targeted_peak_idx = peak_table.groupby(['target_gene'])['target_gene_distance'].transform(min) == peak_table['target_gene_distance']
            peak_table = peak_table.loc[targeted_peak_idx]
            peak_table.reset_index(inplace=True)
        else:
            peak_table.reset_index(inplace=True)
        print '%s gene targeted peak was processed!'%len(peak_table)
    # --------------------------------------------------------------------------
    #  Motifscan core
    # --------------------------------------------------------------------------
    print '#########################################################'
    print '                  Motif Scanning'
    print '#########################################################'
    print 'Motif scanning on %s peaks regions...' % region
    peak_result, tarnum_col, score_col = motifscan_on_peaks(peak_table, motif_table, background, tmp_dir)
    peak_result.to_csv(peak_motif_tarnum, index=False, header=True, cols=tarnum_col)
    peak_result.to_csv(peak_motif_score, index=False, header=True, cols=score_col, float_format='%.2f')
    #peak_result.to_pickle("%s/peak_result.pkl" % output_dir) # only for testing
    if region == 'whole':
        peak_result["common_or_unique"] = peak_table.common_or_unique
        peak_result["is_promoter"] = peak_table.is_promoter
    if extract_target_site:
        export_target_site_info(motif_table, peak_result, genome_db_path, tmp_dir, motif_tarsite_out_dir)

    if not is_enrichment or 'value' not in peak_table.columns:
        core.target_site_distribution(peak_result, motif_table, plot_out_dir,region_radius=peak_length/2)

    if not is_enrichment or len(peak_table) < 50:
        if len(peak_table) < 50:
            print 'Motifscan finished! The number of peaks must be greater than 100 if you want the enrichment analysis performed!'
        shutil.rmtree(tmp_dir)
        return peak_result

    print '#########################################################'
    print '                 Enrichment Analysis'
    print '#########################################################'
    print 'Motif scanning on random control sequences...'
    # ------------------------------------------------------------------------------------
    #  generate random sequence
    # ---------------------------------------------------------------------------------
    if is_enrichment and len(peak_table) >= 50:
        print 'Generating random control based on %s %s peaks...' % (len(peak_table), region),
        sys.stdout.flush()
        if isinstance(gene_table,pd.DataFrame):
            rnd_table = peak.generate_random_with_ref2(gene_table, peak_table, genome_db_path, random_times)
        else:
            rnd_table = peak.generate_random_without_ref(peak_table, genome_db_path,  chromosome_size, random_times)
        print '%d random sequences are processed!' % len(rnd_table)
    # ------------------------------------------------------------------------------------
    #  motif scan on random control
    # -------------------------------------------------------------------------------
    core.motif_scan(rnd_table, motif_table, background, tmp_dir)
    rnd_result = {}
    for idx, motif_record in motif_table.iterrows():
        name = motif_record['name']
        tmp_table = pd.read_pickle('%s/%s' % (tmp_dir, idx))
        rnd_result['%s.tarnum' % name] = tmp_table['%s.tarnum' % name]
    rnd_result = pd.DataFrame(rnd_result)
    # ------------------------------------------------------------------------------------
    # motif enrichment analysis
    # -------------------------------------------------------------------------------------
    print 'Doing enrichment...'
    if region == 'gene2':
        gene_based = True
    else:
        gene_based = False
    enrich_result = core.target_enrichment(peak_result, rnd_result, motif_table, gene_based)
    enrich_result.sort(columns=['enrich_pvalue'], inplace=True)
    enrich_result.to_csv(enrichment_csv, index=False, cols=['name', 'target_number', 'rnd_target_number', 'fold_change', 'enrich_pvalue', 'deplete_pvalue', 'pvalue_corrected'])
    if 'value' in peak_table.columns:
        core.tarnum_and_tarsite_distribution(peak_result, rnd_result, enrich_result, plot_out_dir,region_radius=peak_length/2)

    shutil.rmtree(tmp_dir)
    print '############## Finished! ##################'
    return peak_result, rnd_result, enrich_result


def motifscan_on_peaks(peak_table, motif_table, background, tmp_dir):
    core.motif_scan(peak_table, motif_table, background, tmp_dir)
    # motifscan result output
    peak_result = {}
    peak_result['chr'] = peak_table['chr']
    peak_result['start'] = peak_table['start']
    peak_result['end'] = peak_table['end']
    peak_result['seq_start'] = peak_table['seq_start']
    peak_result['seq_end'] = peak_table['seq_end']
    peak_result['summit'] = peak_table['summit']
    if 'target_gene' in peak_table.columns:
        peak_result['target_gene'] = peak_table['target_gene']
        peak_result['target_gene_distance'] = peak_table['target_gene_distance']
        tarnum_col = ['chr', 'start', 'end', 'summit', 'target_gene','target_gene_distance']
        score_col = ['chr', 'start', 'end', 'summit', 'target_gene','target_gene_distance']
    else:
        tarnum_col = ['chr', 'start', 'end', 'summit']
        score_col = ['chr', 'start', 'end', 'summit']
    if 'value' in peak_table.columns:
        peak_result['value'] = peak_table['value']
        tarnum_col.append('value')
        score_col.append('value')
    for idx, motif_record in motif_table.iterrows():
        name = motif_record['name']
        tmp_table = pd.read_pickle('%s/%s' % (tmp_dir, motif_record['id']))
        peak_result['%s.tarnum' % name] = tmp_table['%s.tarnum' % name]
        peak_result['%s.ratio' % name] = tmp_table['%s.ratio' % name]
        peak_result['%s.taridx' % name] = tmp_table['%s.tarsite' % name]
        peak_result['%s.tarsite' % name] = tmp_table['%s.tarsite' % name]
        peak_result['%s.tarratio' % name] = tmp_table['%s.tarratio' % name]
        tarnum_col.append('%s.tarnum' % name)
        score_col.append('%s.ratio' % name)

    peak_result = pd.DataFrame(peak_result)
    return peak_result, tarnum_col, score_col


def export_target_site_info(motif_table, peak_result, genome_db_path,
                            tmp_dir, motif_tarsite_out_dir):
    n_motif = len(motif_table)
    cnt = 0
    for idx, motif_record in motif_table.iterrows():
        name = motif_record['name']
        sys.stdout.write("Collecting motif target site information for %s...%s%s   \r" % (name, cnt*100/n_motif, "%"))
        sys.stdout.flush()
        motif_len = np.shape([motif_record['matrix']])[2]
        tmp_table = pd.read_pickle('%s/%s' % (tmp_dir, motif_record['id']))
        # output details
        fo = open('%s/%s_target_site.bed' % (motif_tarsite_out_dir, re.sub(r'[\:\-\.]{1,2}', '_', name)), 'w')
        for p_idx, taridx in tmp_table['%s.tarsite' % name].iteritems():
            if len(taridx) > 0:
                for i, start in enumerate(taridx):
                    tar_chr = peak_result.iloc[p_idx]['chr']
                    tar_start = int(peak_result.iloc[p_idx]['seq_start'] + start)
                    tar_end = int(tar_start + motif_len)
                    tar_seq = peak.extract_sequence(genome_db_path, tar_chr, tar_start, tar_end)
                    tar_ratio = tmp_table['%s.tarratio' % name].iloc[p_idx][i]
                    fo.write('%s\t%s\t%s\t%s\t%s\n' % (tar_chr, tar_start, tar_end, tar_seq, tar_ratio))
        fo.close()
        cnt += 1
    sys.stdout.write("Collecting motif target site information...Done!\n")

def merge_two_results(peak_result_1, peak_result_2, rnd_result_1, rnd_result_2, motif_table):
    merged_peak_result = pd.concat([peak_result_1, peak_result_2], ignore_index=True)
    merged_rnd_result = pd.concat([rnd_result_1, rnd_result_2], ignore_index=True)
    merged_enrich_result = core.target_enrichment(merged_peak_result, merged_rnd_result, motif_table)
    return [merged_peak_result, merged_rnd_result, merged_enrich_result]


def merge_two_results_no_rnd(peak_result_1, peak_result_2):
    merged_peak_result = pd.concat([peak_result_1, peak_result_2], ignore_index=True)
    return merged_peak_result
