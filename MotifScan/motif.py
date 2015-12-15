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
# Limitations under the License.

import pandas as pd
import numpy as np
import re
import os
import os.path
from scipy import signal
import peak
import time
import shutil
import sys
import ctypes
import MotifScan


class MAT(ctypes.Structure):
    _fields_ = [("n", ctypes.c_int),
                ("a_arr", ctypes.POINTER(ctypes.c_double)),
                ("c_arr", ctypes.POINTER(ctypes.c_double)),
                ("g_arr", ctypes.POINTER(ctypes.c_double)),
                ("t_arr", ctypes.POINTER(ctypes.c_double))]


def arr2MAT(arr):
    n = np.shape(arr)[1]
    a = np.ctypeslib.as_ctypes(arr[0])
    c = np.ctypeslib.as_ctypes(arr[1])
    g = np.ctypeslib.as_ctypes(arr[2])
    t = np.ctypeslib.as_ctypes(arr[3])
    return MAT(n, a, c, g, t)
'''
motif table:
-----------------------------------------------------------------
name    matrix    max_score   score_cutoff_1    score_cutoff_2 ....
------------------------------------------------------------------
'''


class InvalidMotifMatrixFormat(Exception):
    pass


def load_motif_matrix(matrix):
    motif_id = []
    motif_name = []
    motif_matrix = []
    with open(matrix,'r') as fi:
        matrix_list = []
        for line in fi:
            line = line.strip()
            if line[0] == '>':
                result = re.findall(r'[\w:./()\-_]+',line)
                if len(result) < 2:
                    raise InvalidMotifMatrixFormat("Invalid Motif Name. Motif ID \\t MotifName")
                name = ' '.join(result[1:])
                motif_id.append(result[0])
                motif_name.append(name)
                if matrix_list:
                    # normalize the matrix
                    matrix_list = np.asarray(matrix_list,dtype=float)
                    matrix_list = motif_matrix_normalization(matrix_list)
                    motif_matrix.append(matrix_list)
                    matrix_list = []
            else:
                matrix_list.append([float(i) for i in line.split('\t')])
        if len(set(motif_id)) != len(motif_id): #ensure motif id unique
            raise InvalidMotifMatrixFormat("Motif ID should be unique")
        matrix_list = np.asarray(matrix_list,dtype=float)
        # normalize the matrix
        matrix_list = motif_matrix_normalization(matrix_list)
        motif_matrix.append(matrix_list)

    motif_table = pd.DataFrame({'id':motif_id,'name':motif_name,'matrix':motif_matrix})
    return motif_table


def load_motif(motif_list_path, motif_path):
    '''extract all motif information from the directory to form a motif table('pickle table')'''
    motif_table = read_motif_table(motif_path)
    if motif_list_path == "":
        dup_name_idx = []
        for idx, motif in motif_table.iterrows():
            if len(motif_table.loc[motif_table["name"]==motif['name']])>1:
                dup_name_idx.append(idx)

        for idx in dup_name_idx:
            motif = motif_table.loc[idx]
            motif_table.loc[idx,'name']= "%s_%s"%(motif['name'], motif['id'])
        return motif_table
    motif_list = pd.read_csv(motif_list_path, header=None, names=['motif'])['motif'].tolist()
    motif_list = list(set(motif_list)) #remove duplicate
    motif_idx = []
    motif_db_name = motif_table['name'].tolist()
    motif_table_filter = None
    for motif in motif_list:
        if motif in motif_db_name:
            tmp_motif = motif_table.loc[motif_table['name']==motif]
            if len(tmp_motif) == 0:
                print '%s is not found!'
                continue
            elif len(tmp_motif) > 1:
                for idx, mtf in tmp_motif.iterrows():
                    tmp_motif.loc[idx,'name'] = '%s_%s'%(mtf['name'] ,idx)
            if motif_table_filter:
                motif_table_filter = motif_table_filter.append(tmp_motif)
            else:
                motif_table_filter = pd.DataFrame(tmp_motif)
    return motif_table_filter


def motif_matrix_normalization(matrix):
    normalized_matrix = np.zeros(np.shape(matrix))
    n = np.shape(matrix)[1]
    for i in xrange(n):
        cur_base = matrix[:,i]
        n_zero = len(cur_base[cur_base<0.001])
        if n_zero == 0:
            cur_base = cur_base / cur_base.sum()
        elif n_zero == 1:
            cur_base[cur_base<0.001] = cur_base.sum()/999
            cur_base = cur_base / cur_base.sum()
        elif n_zero == 2:
            cur_base[cur_base<0.001] = cur_base.sum()/998
            cur_base = cur_base / cur_base.sum()
        elif n_zero == 3:
            cur_base[cur_base<0.001] = cur_base.sum()/997
            cur_base = cur_base / cur_base.sum()
        else:
            print 'Invalid position weight matrix! Error: whole column elements are zero!'
            exit(2)
        normalized_matrix[:,i] = cur_base
    return normalized_matrix


def simulation_core(smatrix,M,max_score,B):
    motif_len = np.shape(M)[1]
    my_matirx = smatrix[:,:motif_len]   # only doing motifscan on the left genome
    seq_len = np.shape(my_matirx)[1]
    
    M[M<0.001] = 0.001
    logM = np.log(M)
    revlogM = logM[:,::-1]

    #scan plus strand
    my=np.zeros((1,seq_len),np.int8)
    by=np.zeros((1,seq_len),np.int8)

    for i in range(0,4):
        logbg = np.squeeze(np.log(B[i])*np.ones((1,motif_len),np.int8)) 
        my = my + signal.lfilter(revlogM[i,:],1,my_matirx[i,:])
        by = by + signal.lfilter(logbg,1,my_matirx[i,:])

    ratio_plus = my - by
    ratio_plus = np.squeeze(ratio_plus)
    ratio_plus[0:motif_len-1] = -10000
    
    #scan minus strand
    my=np.zeros((1,seq_len),np.int8)
    by=np.zeros((1,seq_len),np.int8)

    for i in range(0,4):
        logbg = np.squeeze(np.log(B[3-i])*np.ones((1,motif_len),np.int8))
        my = my + signal.lfilter(logM[3-i,:],1,my_matirx[i,:])
        by = by + signal.lfilter(logbg,1,my_matirx[i,:])

    ratio_minus = my - by
    ratio_minus = np.squeeze(ratio_minus)
    ratio_minus[0:motif_len-1] = -10000 

    ratio_both = np.maximum(ratio_plus, ratio_minus)/max_score;
    ratio_both = np.squeeze(ratio_both)
    ratio_both = max(ratio_both)
    return ratio_both


def compute_max_score(motif_table, background):
    logB = np.log(background)
    max_score = []
    for i in motif_table.index:
        matrix = motif_table['matrix'][i]
        max_value = np.amax(matrix, axis=0)
        logM = np.log(max_value)  # the consensus sequence
        max_index = matrix.argmax(axis=0)
        logBm = np.zeros(len(max_index))
        for j in np.arange(len(max_index)):
            logBm[j] = logB[max_index[j]]  # the background sequence
        max_score.append(np.sum(logM-logBm))
    motif_table['max_score'] = max_score
    return motif_table


def generate_consensus(motif_matrix):
    max_value = np.argmax(motif_matrix, axis=0)
    seq = ''
    seq_c = ''
    for i in max_value:
        if i == 0:
            seq += 'A'
            seq_c += 'T'
        elif i == 1:
            seq += 'C'
            seq_c += 'G'
        elif i == 2:
            seq += 'G'
            seq_c += 'C'
        elif i == 3:
            seq += 'T'
            seq_c += 'A'
    seq_c = seq_c[::-1]
    return seq, seq_c


def simulation(motif_table, sample_number, genome_db_path, chrom_size, background):
    # prepare the data to compute the cutoff
    score_c = ctypes.CDLL('%s/score_c.so' % os.path.dirname(os.path.realpath(MotifScan.__file__)))
    score_c.motif_scan_core_simulation.restype = ctypes.c_double
    score_c.motif_scan_core_simulation.argtypes = [ctypes.POINTER(MAT),
                                                   ctypes.POINTER(MAT),
                                                   ctypes.POINTER(ctypes.c_double*4),
                                                   ctypes.c_double]
    motif_max_len = 0
    for i in motif_table.matrix:
        if motif_max_len < np.shape(i)[1]:
            motif_max_len = np.shape(i)[1]

    tmp_dir = '.motif_tmp_%s' % int(time.time())
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    try:
        ############################################################################
        print 'Sequence sampling on the whole genome...',
        sys.stdout.flush()
        random_peak_table = sampling_on_genome(sample_number, motif_max_len, chrom_size)
        n_peak = len(random_peak_table)
        genome_iter = [genome_db_path]*n_peak
        random_peak_table['seq'] = map(peak.extract_sequence,
                                       genome_iter,
                                       random_peak_table['chr'],
                                       random_peak_table['seq_start'],
                                       random_peak_table['seq_end'])
        random_peak_table['seq_matrix'] = map(peak.construct_sequence_matrix, random_peak_table['seq'])
        print 'Done!'
        #######################################################################################

        #######################################################################################
        print 'Start simulating...',
        sys.stdout.flush()
        background = np.ctypeslib.as_ctypes(background)

        n_motif = len(motif_table)
        cnt = 0
        for idx, motif_record in motif_table.iterrows():
            motif_id = motif_record['id']
            max_score = ctypes.c_double(motif_record['max_score'])
            mmatrix = arr2MAT(motif_record['matrix'])
            sys.stdout.write("Scanning for %s...%s%s\
                          \r" % (motif_record['name'], cnt*100/n_motif, "%"))
            sys.stdout.flush()
            ratio = []
            for p_idx, peak_record in random_peak_table.iterrows():
                smatrix = arr2MAT(peak_record['seq_matrix'])
                c = score_c.motif_scan_core_simulation(ctypes.byref(smatrix),
                                                       ctypes.byref(mmatrix),
                                                       ctypes.byref(background),
                                                       max_score)
                ratio.append(c)
            ratio_df = pd.DataFrame({'%s' % motif_id: ratio})
            ratio_df.to_pickle("%s/%s" % (tmp_dir, motif_id))
            cnt += 1

        print 'Done!'
        #######################################################################################

        peak_ratio = {}
        cnt = 0
        for motif_id in motif_table['id']:
            cnt += 1
            tmp_table = pd.read_pickle('%s/%s' % (tmp_dir, motif_id))
            tmp_series = tmp_table[motif_id].copy()
            tmp_series.sort(ascending=False)
            tmp_series = np.array(tmp_series)
            peak_ratio[cnt] = tmp_series
        peak_ratio = pd.DataFrame(peak_ratio)

        cutoff = sample_number
        cnt = 0
        while cutoff > 1:
            cutoff = int(cutoff*0.1)
            cnt += 1
            motif_table['score_cutoff_%s' % cnt] = np.array(peak_ratio.ix[cutoff-1])
        motif_table['score_cutoff'] = motif_table['score_cutoff_4']

        motif_table = motif_table.set_index(motif_table['id'])
        return motif_table
    except Exception, e:
        print 'Simulation failed!', 'Exception: ', e
    finally:
        shutil.rmtree(tmp_dir)


def sampling_on_genome(nsamp, samp_len, chrom_size):
    """
    @param nsamp: sampling times
    @param samp_len: sequence length of each time sampling
    @param chrom_size: chromosome size of the genome
    """
    chr_ls = []
    start_ls = []
    end_ls = []

    n_chr = len(chrom_size)
    for i in np.arange(nsamp):
        chr_idx = np.random.randint(n_chr)
        start = np.random.randint(chrom_size['size'][chr_idx]-samp_len)
        end = start + samp_len
        chr_ls.append(chrom_size['chr'][chr_idx])
        start_ls.append(start)
        end_ls.append(end)
    return pd.DataFrame({'chr': chr_ls, 'seq_start': start_ls, 'seq_end': end_ls})


def write_motif_table(motif_table, output_file_prefix):
    for i in [2, 3, 4, 5, 6]:
        f = open("%s_1e-%s.txt" % (output_file_prefix, i), 'w')
        for j, rec in motif_table.iterrows():
            f.write(">%s\t%s\tmax_score:%s\tscore_cutoff:%s\n" %
                    (rec['id'], rec['name'], rec['max_score'], rec['score_cutoff_%s' % i]))
            f.write('%s\n' % '\t'.join(map(str, rec['matrix'][0])))
            f.write('%s\n' % '\t'.join(map(str, rec['matrix'][1])))
            f.write('%s\n' % '\t'.join(map(str, rec['matrix'][2])))
            f.write('%s\n' % '\t'.join(map(str, rec['matrix'][3])))
        f.close()


def read_motif_table(motif_path):
    f = open(motif_path, 'r')
    motif_id = []
    name = []
    max_score = []
    matrix = []
    score_cutoff = []
    cnt = 1
    for line in f:
        line = line.strip('')
        line = line.strip('\n')
        if cnt % 5 == 1:  # score line
            scores = line.split('\t')
            motif_id.append(scores[0][1:])
            name.append(scores[1])
            max_score.append(float(scores[2].split(":")[1]))
            score_cutoff.append(float(scores[3].split(":")[1]))
            cnt += 1
        else:
            pr = line.split('\t')
            pr = map(float, pr)
            if cnt % 5 == 2:
                a = np.array(pr)
            elif cnt % 5 == 3:
                c = np.array(pr)
            elif cnt % 5 == 4:
                g = np.array(pr)
            elif cnt % 5 == 0:
                t = np.array(pr)
                matrix.append(np.array([a, c, g, t]))
            cnt += 1
    f.close()
    motif_table = pd.DataFrame({'id': motif_id, 'name': name, 'max_score': max_score,
                                'matrix': matrix, 'score_cutoff': score_cutoff})
    motif_table.set_index(motif_table['id'], inplace=True, drop=True)
    return motif_table
