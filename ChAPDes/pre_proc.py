# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:22:30 2018

@author: Link
"""


#%%
import csv
import json
import glob
import os
import random
from pprint import pprint
from bisect import bisect
from collections import namedtuple
from bs4 import BeautifulSoup
#%%
def read_sd_csv(filename):
    """Read segmental duplication information, then return a segmental duplication
    list."""
    with open(filename) as fh:
        sd_csv = csv.reader(fh)
        sd_lst = []
        header = next(sd_csv)
        Row = namedtuple('Row', header)
        for line in sd_csv:
            row = Row(*line)
            sd_lst.append(row)
    return sd_lst


def sequence_overlapped_times(sd_lst, chrom):
    """Return a list of sequence overlapped times by specific chromosome. (Count Problem)"""
    intervals = []
    for row in sd_lst:
        if row.chrom == chrom:
            intervals.append((int(row.chromStart), int(row.chromEnd)))
    max_end = 0
    for row in intervals:
        if row[1] > max_end:
            max_end = row[1]
    times = [0] * (max_end+1)  # 0-based
    for row in intervals:
        for i in range(row[0], row[1]+1):
            times[i] += 1
    return times


def read_chromosome(filename):
    """Read chromosome sequence in a fasta file, then reutrn a sequence dict."""
    seq = []
    count = 0
    with open(filename) as fh:
        for line in fh:
            if '>' in line:
                count += 1
                continue
            seq.append(line.rstrip())
    seq = ''.join(seq)
    assert count == 1, 'sequence number greater than 1.'
    return seq

    
def extract_non_overlapped_sequences(seq, chrom_start, chrom_end, times, chrom, templ_len, rc=False):
    """Return a list of non-overlapped sequences, given chromosome sequnce seq (str),
    sequence overlapped times (list), reverse complementary rc(boolean), and PCR amplified template length (int). (Filter Problem)
    """
    pos_lst = []
    s = 0
    e = 0
    for i in range(len(times)):
        if times[i] != 1:
            e = i - 1  
            if e - s >= templ_len:
                pos_lst.append((s, e))
            s = i + 1  
        if i == len(times)-1 and times[i] == 1 and i-s >= templ_len-1:
            pos_lst.append((s, i))
    
    if pos_lst == [] and sum(times) == len(times):
        pos_lst.append((0, len(times)))
    seqs = []
    if chrom_start == 1 and chrom_end == len(seq) + 1:
        for p in pos_lst:
            name = '>' + chrom + ':' + str(p[0]+1) + '-' + str(p[1]+1)
            subseq = seq[p[0]:p[1]+1].upper()
            if rc:
                subseq = reverse_complement(subseq)
            seqs.append((name, subseq))
    else:
        name = '>' + chrom + ':' + str(chrom_start) + '-' + str(chrom_end)
        subseq = seq[chrom_start-1:chrom_end].upper()
        if rc:
            subseq = reverse_complement(subseq)
        seqs.append((name, subseq))
    return seqs


def fast_primer_pairs(seqs, rcseqs, min_k, max_k, min_len, max_len, rc=False):
    """Fast generate a dict of primer pairs given amplification product length k, and primer length len."""
    pps = {}
    step = step_gen(seqs)
    for k in range(min_k, max_k+1):
        for n in range(len(seqs)):
            chrom = seqs[n][0].split(':')[0][1:]
            pos = int(seqs[n][0].split(':')[1].split('-')[0])
            for i in range(0, len(seqs[n][1])-k+1, step):  
                for j in range(min_len, max_len+1):
                    fp = seqs[n][1][i:i+j]
                    if rc:  
                        rp = rcseqs[n][1][-i-k:-i-(k-j)]
                    else:
                        rp = seqs[n][1][i+k-j:i+k]
                    pps[(chrom, pos, k, j)] = (fp, rp)
                pos += 1
    return pps

def reverse_complement(dna):
    table = str.maketrans('ATCGN', 'TAGCN')
    dna = dna.translate(table)
    return dna[::-1]
    
def step_gen(seqs):
    length = sum([len(e[1]) for e in seqs])
    if length > 1000000:
        step = length // 1000000
    else:
        step = 1
    return step    
    
    
def write_non_overlapped_sequences(seqs, filename):
    """Write sequences seqs to a fasta file named filename."""
    with open(filename, 'w') as fh:
        for s in seqs:
            fh.writelines(s[0]+'\n')
            fh.writelines(s[1]+'\n')


def write_primer_pairs(pps, path_prefix):
    """Write primer pairs to file separately."""
    with open(path_prefix+'.1.fa', 'w') as fh:
        for k, v in pps.items():
            if isinstance(k, tuple):
                fh.writelines('>'+k[0]+':'+str(k[1])+':'+str(k[2])+':'+str(k[3])+'\n')  
            elif isinstance(k, str):
                fh.writelines('>'+k+'\n')
            else:
                raise ValueError('pps keys are neither tuple nor string')
            fh.writelines(v[0]+'\n')
    with open(path_prefix+'.2.fa', 'w') as fh:
        for k, v in pps.items():
            if isinstance(k, tuple):
                fh.writelines('>'+k[0]+':'+str(k[1])+':'+str(k[2])+':'+str(k[3])+'\n')
            elif isinstance(k, str):
                fh.writelines('>'+k+'\n')
            else:
                raise ValueError('pps keys are neither tuple nor string')
            fh.writelines(v[1]+'\n')


def random_select(pps, num=100000, seed=None):
    """Return a dict of sample given primer pairs pps (dict) and population num."""
    random.seed(seed)
    return dict(random.sample(list(pps.items()), k=num))


def write_json(lst, path, filename):
    if not os.path.exists(path):
        os.mkdir(path)
    with open(path+filename+'.json', 'w') as fh:
        json.dump(lst, fh, indent=4)


def read_json(path, chrom):
    for fn in glob.glob(path+chrom+'*.json'):
        with open(fn) as fh:
            lst = json.load(fh)
    return lst


def filter_primer_tandem_repeat(pps, chrom):
    """ Filter the primer pairs with TR """
    tandem_repeats = read_json('raw_data/tandem_repeats/', chrom)
    tr_l = [e[0] for e in tandem_repeats]
    tr_r = [e[1] for e in tandem_repeats]
    f_pps = {}
    for k, v in pps.items():
        i = bisect(tr_l, k[1]) - 1
        if k[1] < tr_r[i] and k[1] + k[3] - 1 <= tr_r[i]:
            continue
        start = k[1] + k[2] - k[3]
        i = bisect(tr_l, start) - 1
        if start < tr_r[i] and start + k[3] - 1 <= tr_r[i]:
            continue
        f_pps[k] = v
    return f_pps


def filter_primer_snp(pps, popul, chrom):
    """Discard the primer pairs obtaining SNP"""
    snp_info = read_json('raw_data/snp_info/'+popul+'/', chrom)
    snp_dict = dict([(int(e[1]), e[0]) for e in snp_info])
    f_pps = {}
    for k, v in pps.items():
        has_snp = False
        for i in range(k[1], k[1]+k[3]):
            if i in snp_dict:
                has_snp = True
                break
        if not has_snp:
            for i in range(k[1]+k[2]-k[3], k[1]+k[2]):
                if i in snp_dict:
                    has_snp = True
                    break
        if not has_snp:
            f_pps[k] = v
    return f_pps


def find_tandem_repeats(seq):
    """ Find the TR positon list in chromosome."""
    pos = []
    i = 0
    while i < len(seq):
        if seq[i].islower():
            start = i + 1
            for j in range(start, len(seq)):
                if seq[j].isupper():
                    pos.append((start, j))
                    i = j
                    break               
        i += 1
    return pos


def extract_tandem_repeats(path='raw_data/hg19/'):
    """ Write the TR positon list to the json file."""
    for fn in glob.glob(path+'*.fa'):
        print('Processing file:'+fn+'...')
        seq = read_chromosome(fn)
        pos = find_tandem_repeats(seq)
        file = fn.split('/')[-1].split('.')[0]
        print(file)
        print('Writing file: '+file+'.json...')
        write_json(pos, 'raw_data/tandem_repeats/', file)

    
def extract_snp_info(path='raw_data/1000genomes/', popul='EAS', rf=0.1):
    """ Write the SNP information to the json file."""
    for fn in glob.glob(path+'*.vcf'):
        print('Processing file:'+fn+'...')
        snp_abs = filter_snp(fn, popul, rf)
        file = 'chr'+snp_abs[0][0]
        print('Writing file: '+file+'.json...')
        write_json(snp_abs, 'raw_data/snp_info/'+popul+'/', file)
 
 
def filter_snp(filename, popul, rf):
    """Return the SNP positon based on different type of population and rf value. """
    with open(filename) as fh:
        tsv_file = csv.reader(fh, delimiter='\t')
        snp_abs = []
        for row in tsv_file:
            if len(row) < 8:
                continue
            info = row[7].split(';')
            if len(info) < 5:
                continue
            temp = info[[e.startswith(popul) for e in info].index(True)].split('=')[1]
            try:
                af = float(temp)
            except:
                af = float(temp.split(',')[1])
            if af >= rf:
                snp_abs.append((row[0], row[1]))
    return snp_abs



def main():
    pass

    
if __name__ == "__main__":
    main()