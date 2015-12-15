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
import numpy as np
import re
import pandas as pd
import math
import os
'''
read a genome file(.fa) and compute the chromosome size and genome background

'''

ROME_CHR_NAME_1 = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", "chrXVIII", "chrXIX", "chrXX", "chrXXI", "chrXXII", "chrXXIII")
ROME_CHR_NAME_2 = ("ChrI", "ChrII", "ChrIII", "ChrIV", "ChrV", "ChrVI", "ChrVII", "ChrVIII", "ChrIX", "ChrXI", "ChrXII", "ChrXIII", "ChrXIV", "ChrXV", "ChrXVI", "ChrXVII", "ChrXVIII", "ChrXIX", "ChrXX", "ChrXXI", "ChrXXII", "ChrXXIII")


def compute_genome_background(genome_path):
    '''
    this function will solely take autosome chromosome into account(chromosomes such as chrX and chrY will be excluded).

    return an array as following:
    [ 0.29485  0.20491163  0.20500062  0.29523775 ]
    each float represent the ratio of the base,  which is in the order of 'A', 'C', 'G', 'T'
    '''
    a, c, g, t = 0, 0, 0, 0
    is_autosome = 1
    is_chrX_autosome = 0
    with open(genome_path) as fi:
        for line in fi:
            line = line.strip()

            if line[0] == '>':
                line = line[1:]
                if re.search(r'[cC]hr[0-9]+\b', line) or line in ROME_CHR_NAME_1 or line in ROME_CHR_NAME_2 or ((line == "chrX" or line == "ChrX") and is_chrX_autosome == 1):
                    if line == "ChrIX" or line == "chrIX":
                        is_chrX_autosome = 1
                    print 'Processing %s...' % line
                    is_autosome = 1
                else:
                    print 'Processing %s...Skip' % line
                    is_autosome = 0
            elif is_autosome == 1:
                for base in line:
                    base = base.upper()
                    if base == 'A':
                        a += 1
                    if base == 'C':
                        c += 1
                    if base == 'G':
                        g += 1
                    if base == 'T':
                        t += 1
    background = np.array([a, c, g, t])/float(np.sum([a, c, g, t]))
    return pd.Series({'background': background})


def compute_genome_size(genome_path):
    '''
    return a table with two columns,  chr represent the chromosome name and size represent the size the corresponding chromosome,
    demo of hg19:
          chr       size
0    chr1  249250621
1    chr2  243199373
2    chr3  198022430
3    chr4  191154276
4    chr5  180915260
    '''
    chrom_id_list = []
    seq_size_list = []
    with open(genome_path, 'r') as fi:
        seq = ''
        isFirst = 1
        for line in fi:
            if line[0] == '>':  # chromosome id line
                if isFirst == 1:
                    chrom_id = line.split('\n')[0][1:]
                    isFirst = 0
                else:
                    seq_len = len(seq)
                    chrom_id_list.append(chrom_id)
                    seq_size_list.append(seq_len)
                    chrom_id = line.split('\n')[0][1:]
                    seq = ''
            else:  # sequence line
                seq += line.split('\n')[0]
        seq_len = len(seq)
        chrom_id_list.append(chrom_id)
        seq_size_list.append(seq_len)
        chromosome_size = pd.DataFrame({'chr': chrom_id_list, 'size': seq_size_list})
        return chromosome_size

def split_genome(genome_path, output):
    # split the genome by chromosom
    chr_list = []
    with open(genome_path) as fi:
        seq_buff = ''
        for line in fi:
            if line[0] == '>': # the chromosome id line
                if seq_buff != '':
                    with open('%s/%s'%(output, current_chr), 'w') as fo:
                        fo.write(seq_buff)
                seq_buff = ''
                current_chr = line[1:-1]
                chr_list.append(current_chr)
                seq_buff += line
            else: # the sequence line
                seq_buff += line
        with open('%s/%s'%(output, current_chr), 'w') as fo:
            fo.write(seq_buff)

    # format the genome
    for current_chr in chr_list:
        fi = open('%s/%s' % (output, current_chr))
        fo = open('%s/%s_new' % (output, current_chr), 'w')
        chr_name = fi.readline()
        fo.write(chr_name)
        seq = fi.read().replace('\n', '')
        nseq = len(seq)
        for i in range(int(math.ceil(nseq/50.0))):
            s_idx = 50*i
            e_idx = 50*(i+1)
            if e_idx > nseq:
                segment_to_write = seq[s_idx:e_idx]
                fo.write('%s\n'%segment_to_write)
            else:
                segment_to_write = seq[s_idx:e_idx]
                fo.write('%s\n'%segment_to_write)
        fi.close()
        fo.close()
        os.rename('%s/%s_new'%(output, current_chr), '%s/%s'%(output, current_chr))