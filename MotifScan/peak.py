import numpy as np
import random
import pandas as pd
import re
import os.path
import sys

'''
This set of functions will extract the peak information and generate corresponding peak sequence and sequence matrix,
while they also generate the random control sequence

input peak file
----------------------------------------------------------------------------------------
chr    start    end    summit(optional)    m_value(optional)
-------------------------------------------------------------------------------------------

gene_table
-----------------------------------------------------------------------------------------------
name    id      chr     strand    start     end     TSS   
-----------------------------------------------------------------------------------------------

peak_table
---------------------------------------------------------------------------------------
chr    start    end    summit    m_value    a_value    p_value    seq     seq_matrix
-----------------------------------------------------------------------------------------

rnd_table 
---------------------------------------------------------------
chr  start   seq     seq_matrix
-------------------------------------------------------------
'''

__all__ = ["load_peak"]
def load_ref_gene(gene_path):
    '''
    load the gene annotation, default: refSeq Gene
    '''
    gene_table = pd.read_csv(gene_path,sep='\t',usecols=[1,2,3,4,5])
    gene_table.columns = ['id','chr','strand','start','end']

    def judge_tss(strand,start,end):
        if strand=='+':
            return start
        else:
            return end
    gene_table['TSS'] = map(judge_tss,gene_table['strand'],gene_table['start'],gene_table['end'])
    return gene_table

def extract_sequence(genome,chr,bpstart,bpend):
    # 0-based: bpstart is included while bpend is not included.
    gf = file('%s/%s'%(genome,chr))
    gf.seek(0,0)
    gf.readline()  #read the first line; the pointer is at the second line
    nbp = bpend - bpstart
    offset = bpstart + np.floor(bpstart/50) #assuming each line contains 50 characters; add 1 offset per line
    gf.seek(offset,1)
    seq_tmp = gf.read(nbp+int(np.floor(nbp/50))+1)
    seq_tmp = seq_tmp.replace('\n','')
    gf.close()
    return seq_tmp[0:nbp].upper()


def construct_sequence_matrix(seq):
    matrix = np.zeros((4,len(seq)))
    for base,idx in zip(seq,np.arange(len(seq))):
        if base == 'A':
            matrix[0][idx] = 1
        elif base == 'C':
            matrix[1][idx] = 1
        elif base == 'G':
            matrix[2][idx] = 1
        elif base == 'T':
            matrix[3][idx] = 1
    return matrix


def load_peak(peaks_path,genome_path,peak_len,peak_type):
    '''
    load the peak file, and extract the from peak sequence from the genome
    '''
    skip = _get_comment_row_cnt(peaks_path)  # skip the command lines 
    # load the peak table
    if peak_type == "macs": # macs peaks
        print '%s is a xls file generated by MACS.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', skiprows=skip, 
        names=['chr','start','end','length','summit','tags','value','fold_change','FDR'],
        dtype={'chr':'string'})
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
        peak['summit'] = peak['start']+peak['summit']
    elif peak_type == "manorm": # MAnorm peaks in xls format
        print '%s is a xls file generated by MAnorm.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', usecols=[0,1,2,3,4],names=['chr','start','end','summit','value'], skiprows=skip)
        peak['summit'] = peak['start']+peak['summit']
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
    elif peak_type == "MAnorm_bed": # macs peak in bed format
        print '%s is a bed file generated by MAnorm.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', usecols=[0,1,2,4],names=['chr','start','end','value'], skiprows=skip)
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
        peak['summit'] = peak['start'] + (peak['end']-peak['start'])/2
    elif peak_type == "bed3col":# general peak bed file, 3~5 columns
        print '%s is a 3-column bed file.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', skiprows=skip,names=['chr','start','end'])
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
        peak['summit'] = peak['start'] + (peak['end']-peak['start'])/2 
    elif peak_type == "bed4col":# general peak bed file, 3~5 columns
        print '%s is a 4-column bed file.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', skiprows=skip,names=['chr','start','end','summit'])
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
        if peak['summit'].dtype !='int64':
                print 'Invalid summit! Please use option -f to specify your peak format!'
                exit(0)
        peak['summit'] = peak['start'] + peak['summit'] 
        if (peak['summit']>peak['end']).any():
            print 'Invalid summit Please use option -f to specify your peak format!'
            exit(0) 
    elif peak_type == "bed5col":# general peak bed file, 3~5 columns
        print '%s is a 5-column bed file.'%os.path.basename(peaks_path)
        peak = pd.read_csv(peaks_path, '\t', skiprows=skip,names=['chr','start','end','summit','value'])
        if len(peak)==0:
            print 'Warning: No peaks detected in peak file. Exiting...'
            exit(2)
        if peak['summit'].dtype != 'int64':
                print 'Invalid summit! Please use option -f to specify your peak format!'
                exit(0)
        peak['summit'] = peak['start'] + peak['summit'] 
        if (peak['summit']>peak['end']).any():
            print 'Invalid summit! Please use option -f to specify your peak format!'
            exit(0)

    else:
        print 'Invalid peak file format! Please use option -f to specify your peak format!'
        exit(1)
    peak['seq_start'] = peak['summit'] - peak_len/2
    peak.ix[peak['seq_start']<0,'seq_start'] = 0
    peak['seq_end'] = peak['summit'] + peak_len/2
    n_peak = len(peak)
    genome_iter = [genome_path]*n_peak
    peak['seq'] = map(extract_sequence,genome_iter,peak['chr'],peak['seq_start'],peak['seq_end'])
    peak['seq_matrix'] = map(construct_sequence_matrix,peak['seq'])

    return peak

def _get_comment_row_cnt(txt_file):
    """
    @return: the count of comment lines in the header
    last modified: 2014-10-13
    """
    cnt = 0 
    with open(txt_file) as fi:
        for line in fi:
            if line.strip() == '' or line.strip()[0] == '#' or re.search(r'[Cc][Hh][Rr]\Z',line.strip().split('\t')[0]): # the comment lines
                cnt += 1
            else:
                break 
    return cnt

def generate_random_without_ref(peak_table,genome_path,chromosome_size,random_times):
    # generate the random contral
    rnd_chr = []
    rnd_start = []
    rnd_end = []
    peak_length = peak_table['seq_end'].iloc[0]-peak_table['seq_start'].iloc[0]

    npeak = len(peak_table)
    for i in np.arange(npeak):
        for j in np.arange(random_times):
            cur_chr = peak_table['chr'].iloc[i]
            cur_start = random.randint(0,chromosome_size[chromosome_size['chr']==cur_chr]['size']-peak_length)
            cur_end = cur_start + peak_length
            rnd_chr.append(cur_chr)
            rnd_start.append(cur_start)
            rnd_end.append(cur_end)
    rnd_table = pd.DataFrame({'chr':rnd_chr,'start':rnd_start,'end':rnd_end}) 
    n_rnd = len(rnd_table) 
    genome_iter = [genome_path]*n_rnd
    rnd_table['seq'] = map(extract_sequence,genome_iter,rnd_table['chr'],rnd_table['start'],rnd_table['end'])
    rnd_table['seq_matrix'] = map(construct_sequence_matrix,rnd_table['seq'])
    rnd_table.to_pickle('rand_table.pkl')
    return rnd_table

def generate_random_with_ref(gene_table,peak_table,genome_path,chromosome_size,random_times):

    npeak = len(peak_table)
    ngene = len(gene_table)

    a = np.array([np.array(peak_table['summit'])]*ngene)
    b = np.array([np.array(gene_table['TSS'])]*npeak).transpose()
    distance_matrix = a - b
    abs_distance_matrix = np.absolute(distance_matrix)
    min_distance_idx = np.argmin(abs_distance_matrix,axis=0)

    d_list = []
    for i in np.arange(npeak):
        d_list.append(distance_matrix[min_distance_idx[i],i])

    rnd_chr = []
    rnd_start = []
    rnd_end = []
    chr_size_dict = {}
    chrs = gene_table['chr'].unique()
    for chr_id in chrs:
        if chr_id in peak_table['chr']:
            chr_size_dict[chr_id] = max(peak_table[peak_table['chr']==chr_id]['end'].max(),gene_table[gene_table['chr']==chr_id]['end'].max())
        else:
            chr_size_dict[chr_id] = gene_table[gene_table['chr']==chr_id]['end'].max()

    for d in d_list:
        rnd_idx = random.sample(np.arange(len(gene_table)),random_times) # for each peak, generate 5 samples randomly
        for i in rnd_idx:
            chr_name = gene_table.iloc[i]['chr']
            chr_size = chr_size_dict[chr_name]
            rnd_chr.append(gene_table.iloc[i]['chr'])
            if gene_table.iloc[i]['strand'] == '+':
                summit = gene_table.ix[i]['TSS'] + d
                summit = summit if summit <= chr_size else chr_size
                summit = summit if summit >= 500 else 500
            else:
                summit = gene_table.ix[i]['TSS'] - d
                summit = summit if summit <= chr_size else chr_size
                summit = summit if summit >= 500 else 500
            start = summit-500+1
            rnd_start.append(summit-500+1)
            end = summit+500 if summit+500 <= chr_size else chr_size
            rnd_end.append(end)
        
    rnd_table = pd.DataFrame({'chr':rnd_chr,'start':rnd_start,'end':rnd_end}) 
    #print rnd_table
    n_rnd = len(rnd_table) 
    genome_iter = [genome_path]*n_rnd
    rnd_table['seq'] = map(extract_sequence,genome_iter,rnd_table['chr'],rnd_table['start'],rnd_table['end'])
    rnd_table['seq_matrix'] = map(construct_sequence_matrix,rnd_table['seq'])
    return rnd_table


def generate_random_with_ref2(gene_table,peak_table,genome_path,random_times):
    '''
    generate the random sequence according the peak summit distance relative to the target gene TSS
    '''
    peak_gene = merge_peak_gene(peak_table, gene_table)
    d_list = []
    chrom_set = peak_gene["chr"].unique()
    for chrom in chrom_set:
        d_list += process_long_distance(distance2target_gene(peak_gene.loc[peak_gene["chr"]==chrom]))
    chromosome_size = chr_margin(peak_table,gene_table)
    rnd_table = peak_random(d_list, gene_table, chromosome_size, random_times)
    
    # generate the random sequence
    n_rnd = len(rnd_table) 
    genome_iter = [genome_path]*n_rnd
    rnd_table['seq'] = map(extract_sequence,genome_iter,rnd_table['chr'],rnd_table['start'],rnd_table['end'])
    rnd_table['seq_matrix'] = map(construct_sequence_matrix,rnd_table['seq'])

    return rnd_table


def merge_peak_gene(peak_table, gene_table):
    """merge peak table and gene table and sorting by their position


    Args:
        peak_table: pandas dataframe
        gene_table: pandas dataframe
    Returns:
        merged table: pandas dataframe with 4 columns: chr, position, strand, type
        Note: type "gene" is for record from gene table, type "peak" is for record from peak table.
              positon for gene is TSS; position for peak is peak summit.
              the index of dataframe are reset by the sorting.

    """
    peak_part = peak_table[['chr','summit']]
    peak_part['strand'] = 'unknown'
    peak_part['type'] = 'peak'
    peak_part.columns = ['chr','position','strand','type']
    n_peak = len(peak_part)

    gene_part = gene_table[['chr','TSS','strand']]
    gene_part['type'] = 'gene'
    gene_part.columns = ['chr','position','strand','type']
    
    peak_gene = peak_part.append(gene_part,ignore_index=True)
    peak_gene.sort(['chr','position'],inplace=True)
    peak_gene.reset_index(inplace=True,drop=True)
    return peak_gene

def distance2target_gene(peak_gene, distance_cutoff_abs = 10000):
    """Compute each peak's distance to its target gene

        Append the column named after distance_to_target_gene in the peak_gene indicating the distance.

        All features should be on the same chromosome

    Args:
        peak_gene: dataframe  by merge_peak_gene's output
        distance_cutoff_abs: when distance to target gene > 10k, the target gene is considered not exist.
    Returns:
        the peak_gene with the column, distance_to_target_gene, that indicates the distance between peak and its target gene
    """

    # initialize the distance to target gene
    peak_gene['distance_to_target_gene'] = -1
    peak_gene.reset_index(inplace=True, drop=True)
    peak_idx = peak_gene.loc[peak_gene['type'] == 'peak'].index
    d_list = []
    for i in peak_idx:
        cur_peak = peak_gene.iloc[i]
        # search for the nearest gene in upstream
        dis_up = float('inf')
        if i > 0: # do not need to search up for the first element
            for j in np.arange(i-1,-1,-1):
                if i-j > distance_cutoff_abs:
                    break
                cur_gene = peak_gene.iloc[j] 
                if cur_gene['type'] == 'peak': # skip peak record
                    continue
                dis_up = cur_peak['position'] - cur_gene['position']
                dis_up = dis_up if cur_gene['strand'] == '+' else - dis_up
                break
        # search for the nearest gene in downstream
        dis_down = float('-inf')
        if i < len(peak_gene):
            for j in np.arange(i+1,len(peak_gene),1):
                if j-i > distance_cutoff_abs:
                    break
                cur_gene = peak_gene.iloc[j] 
                if cur_gene['type'] == 'peak': # skip peak record
                    continue
                dis_down = cur_peak['position'] - cur_gene['position']
                dis_down = dis_down if cur_gene['strand'] == '+' else - dis_down
                break

        d = dis_up if abs(dis_up) < abs(dis_down) else dis_down
        
        d_list.append(d)
    return d_list

def process_long_distance(d_list, distance_cutoff_abs = 10000):
    # distance that larger than this will be replaced by a random distance between it and 10*it
    for i,d in enumerate(d_list):
        if d == float('inf'):
            d_list[i] = np.random.randint(distance_cutoff_abs, 10*distance_cutoff_abs)
        elif d == float('-inf'):
            d_list[i] = -np.random.randint(distance_cutoff_abs, 10*distance_cutoff_abs)
        else:
            d_list[i] = (d/abs(d))*np.random.randint(distance_cutoff_abs, 10*distance_cutoff_abs) if abs(d) >= distance_cutoff_abs else d
    return d_list 

def chr_margin(peak_table,gene_table):
    chr_size_dict = {}
    chrs = gene_table['chr'].unique()
 
    for chr_id in chrs:
        if chr_id in peak_table['chr']:
            chr_size_dict[chr_id] = max(peak_table[peak_table['chr']==chr_id]['end'].max(),gene_table[gene_table['chr']==chr_id]['end'].max())
        else:
            chr_size_dict[chr_id] = gene_table[gene_table['chr']==chr_id]['end'].max()
    return chr_size_dict

def extract_target_gene(peak_table, gene_table, distance_cutoff_abs = 10000):
    target_gene = []
    target_gene_distance = []
    for p_idx, peak_record in peak_table.iterrows():
        distance = distance_cutoff_abs
        gene = 'No Target'
        for g_idx, gene_record in gene_table.loc[gene_table['chr']==peak_record['chr']].iterrows():
            
            if gene_record['strand'] == '+':
                if gene_record['TSS'] - peak_record['summit'] < distance and gene_record['TSS'] - peak_record['summit'] > 0:
                    distance = gene_record['TSS'] - peak_record['summit']
                    gene = gene_record['id']
            else:
                if peak_record['summit'] - gene_record['TSS']  < distance and  peak_record['summit'] - gene_record['TSS'] > 0:
                    distance = peak_record['summit'] - gene_record['TSS'] 
                    gene = gene_record['id']
        target_gene.append(gene)
        target_gene_distance.append(distance)
    peak_table['target_gene'] = target_gene
    peak_table['target_gene_distance'] = target_gene_distance
    return peak_table

def peak_random(d_list, gene_table,chr_size_dict, rnd_times=5):
    """generate rnd peak table

    Args:
        d_list: contain all distances
        gene_table: pandas dataframe, gene table
        chr_size_dict: chromosome margin of each chromosome
        rnd_times: random times for each distance

    Returns:
        pandas dataframe containing random times
    """

    gene_part = gene_table[['chr','TSS','strand']]
    gene_part['type'] = 'gene'
    gene_part.columns = ['chr','position','strand','type']
    
    rnd_chr = []
    rnd_start = []
    rnd_end = []
    for d in d_list:
        rnd_idx = random.sample(np.arange(len(gene_part)),rnd_times) # for each peak, generate 5 samples randomly
        for i in rnd_idx:
            cur_gene = gene_part.iloc[i]
            chr_name = cur_gene['chr']
            
            chr_size = chr_size_dict[chr_name]
            rnd_chr.append(cur_gene['chr'])
            if cur_gene['strand'] == '+':
                summit = cur_gene['position'] + d
                summit = summit if summit <= chr_size else chr_size
                summit = summit if summit >= 500 else 500
            else:
                summit = cur_gene['position'] - d
                summit = summit if summit <= chr_size else chr_size
                summit = summit if summit >= 500 else 500
            start = summit-500+1
            rnd_start.append(summit-500+1)
            end = summit+500 if summit+500 <= chr_size else chr_size
            rnd_end.append(end)
        
    rnd_table = pd.DataFrame({'chr':rnd_chr,'start':rnd_start,'end':rnd_end}) 
    return rnd_table