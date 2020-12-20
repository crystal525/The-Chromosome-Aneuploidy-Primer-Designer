# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 17:33:56 2018

@author: Link
"""


#%%
from pre_proc import *
from glob import glob
import subprocess
import os
#%%



def alignment(task, primer_type):
    """Do alignment by Bowtie2 and BLAST Tool."""
    os.chdir('spec_check/'+task)
    
    bowtie2_cmd = "time bowtie2 -p 8 -f --very-sensitive --mp 1,1 --rdg \
    99,99 --rfg 99,99 --score-min L,-0.35,-0.35 -k 3 -x ../ChAPDes/raw_data/bowtie2_index/hg19 \
    -1 primer_pairs.1.fa -2 primer_pairs.2.fa -S pair_alns.sam &> pairs.log" ### change your code according to your path
    with subprocess.Popen(bowtie2_cmd, shell=True) as proc:
        proc.communicate()
        
    sam_lst = read_tsv('pair_alns.sam')
    sam_lst = sam_lst[sam_seq_index(sam_lst):]
    sam_dict = lst_to_dict(sam_lst)
    align_pairs = sam_filter(sam_dict, primer_type)
    write_filtered_pairs(align_pairs)
        
    wc_cmd = "wc -l primer_pairs.1.filt.fa"
    with subprocess.Popen(wc_cmd, stdout=subprocess.PIPE, shell=True) as proc:
        line_num = proc.communicate()
        line_num = str(line_num[0]).split(' ')[0].split("'")[1]
        
    pad_cmd = "awk 'BEGIN {for (i = 1; i <= %s; ++i) \
    {if (i %% 2 == 0) print \"NNNNNNNNNNNNNNNNNNNN\" ; else print \"\"}}'  > padding.txt" % line_num
    with subprocess.Popen(pad_cmd, shell=True) as proc:
        proc.communicate()
        
    trim_end_cmd = "cat primer_pairs.2.filt.fa | awk '{if (match($1, />/)) print \"\" ; else print $1}' \
    > primer_pairs.2.filt.trim.fa"
    with subprocess.Popen(trim_end_cmd, shell=True) as proc:
        proc.communicate()
        
    paste_cmd = "paste -d'\b' primer_pairs.1.filt.fa padding.txt | paste -d'\b' \
    - primer_pairs.2.filt.trim.fa > merge.fa"
    with subprocess.Popen(paste_cmd, shell=True) as proc:
        proc.communicate()
         
    blast_cmd = "time blastn -db ../ChAPDes/raw_data/blastdb/hg19New -query merge.fa \
    -out blast.out -evalue 30000 -max_hsps 100 -num_threads 16 -outfmt '6 std qseq sseq' \
    -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -perc_identity 65 &> blast.log" ###change your code according to your path
    with subprocess.Popen(blast_cmd, shell=True) as proc:
        proc.communicate()
        
    os.chdir('../../')
 
 
def read_tsv(filename):
    """Read the .sam file."""
    with open(filename) as fh:
        tsv_file = csv.reader(fh, delimiter='\t')
        lst = []
        for row in tsv_file:
            lst.append(row)
    return lst 
    
    
def sam_seq_index(sam_lst):
    """Get the position of useful infomation of the .sam file."""
    for i in range(len(sam_lst)):
        if sam_lst[i][0] == '@PG':
            return i+1

 
def lst_to_dict(lst):
    """Reverse the list format into the dictionary format."""
    res = {}
    for e in lst:
        res.setdefault(e[0], []).append(e[1:])
    return res

    
def are_in_wrong_chroms(chrom1, chrom2, primer_type):
    if primer_type == 'a':
        return chrom1 in ['chrX', 'chrY'] or chrom2 in ['chrX', 'chrY']
    elif primer_type == 's':
        return not (chrom1 in ['chrX', 'chrY'] and chrom2 in ['chrX', 'chrY'])
    elif primer_type == 'xa':
        return chrom1 == 'chrY' or chrom2 == 'chrY' \
    or (chrom1 != 'chrX' and chrom2 != 'chrX')
    elif primer_type == 'ya':
        return chrom1 == 'chrX' or chrom2 == 'chrX' \
    or (chrom1 != 'chrY' and chrom2 != 'chrY')
    else:
        raise ValueError('value of primer_type is wrong: %s' % primer_type)
    
 
def sam_filter(sam_dict, primer_type):
    """Filter SAM files based on rules. """
    
    paths = glob('../ChAPDes/raw_data/hg19/*.fa')###change your code according to your path 
    
    chroms = {}
    for fn in paths:
        chrom = fn.split('/')[-1].split('.')[0]
        chroms[chrom] = read_chromosome(fn)
    
    align_pairs = {} 
    for k, v in sam_dict.items():
        prim_len = int(k.split(':')[-1])
        pps = paired_sam_primers(v)
        if len(pps) != 2: # make sure primers align into only two chromosome
            continue
        if pps[0][0][8] != pps[1][0][8]: 
            continue
        chrom1 = pps[0][0][1]
        chrom2 = pps[1][0][1]
        if chrom1 == chrom2:  # make sure primers in the different chromosome
            continue
        if are_in_wrong_chroms(chrom1, chrom2, primer_type):  # choose different models
            continue
        if [e for e in pps[0][0] if e[:2] == 'NM'][0].split(':')[-1] != '0' \
        or [e for e in pps[0][1] if e[:2] == 'NM'][0].split(':')[-1] != '0' \
        or [e for e in pps[1][0] if e[:2] == 'NM'][0].split(':')[-1] != '0' \
        or [e for e in pps[1][1] if e[:2] == 'NM'][0].split(':')[-1] != '0':  
            continue
        if chrom1 == k.split(':')[0]:  
            i = 0
            j = 1
        else:
            i = 1
            j = 0
        templ_len = int(k.split(':')[-2])
        templ_start1 = int(pps[i][0][2])
        templ_end1 = int(pps[i][1][2]) + prim_len - 1
        templ_start2 = int(pps[j][0][2])
        templ_end2 = int(pps[j][1][2]) + prim_len - 1
        if templ_end1 - templ_start1 + 1 != templ_len \
        or templ_end2 - templ_start2 + 1 != templ_len:  
            continue
        templ1 = chroms[pps[i][0][1]][templ_start1-1:templ_end1].upper()  
        templ2 = chroms[pps[j][0][1]][templ_start2-1:templ_end2].upper()
        if templ1[prim_len:-prim_len] == templ2[prim_len:-prim_len]:  # make sure to have mismatch
            continue
        align_pairs[k] = pps
    return align_pairs
    
def paired_sam_primers(aligns):
    """ Find primer pairs. """
    left_aligns = [e for e in aligns if int(e[7]) > 0]
    right_aligns = [e for e in aligns if int(e[7]) < 0]
    pps = []
    while len(left_aligns) > 0:
        for i in range(len(right_aligns)):
            if left_aligns[0][2] == right_aligns[i][6]:
                pps.append((left_aligns[0], right_aligns[i]))
        del(left_aligns[0])
    return pps


def write_filtered_pairs(align_pairs):
    with open('primer_pairs.1.filt.fa', 'w') as fh:
        for k, v in align_pairs.items():
            if v[0][0][1] == k.split(':')[0]:
                i = 0
            else:
                i = 1
            fh.writelines('>'+k+'\n')
            fh.writelines(v[i][0][8]+'\n')
    with open('primer_pairs.2.filt.fa', 'w') as fh:
        for k, v in align_pairs.items():
            if v[0][1][1] == k.split(':')[0]:
                i = 0
            else:
                i = 1
            fh.writelines('>'+k+'\n')
            fh.writelines(v[i][1][8]+'\n')
            
            
            
def spec_check(blast_dict, primer_type):
    """Check the specificity again of the file output from BLAST tool."""
     paths = glob('../ChAPDes/raw_data/hg19/*.fa')###change your code according to your path
    
    chroms = {}
    for fn in paths:
        chrom = fn.split('/')[-1].split('.')[0]
        chroms[chrom] = read_chromosome(fn)
    
    align_pairs = {}
    align_templs = {}  
    for k, v in blast_dict.items():
        prim_len = int(k.split(':')[-1])
        f_hsps = filter_hsps(v, prim_len)
        pps = paired_primers(f_hsps, prim_len)
         if len(pps) != 2:
            continue
        is_reverse = False
        for e in pps:
            if int(e[1][8]) - int(e[0][7]) < 0:
                is_reverse = True
        if is_reverse:
            continue
        chrom1 = pps[0][0][0]
        chrom2 = pps[1][0][0]
        if chrom1 == chrom2:  # in the same chromosome
            continue
        if are_in_wrong_chroms(chrom1, chrom2, primer_type):
            continue
        if chrom1 == k.split(':')[0]:  
            i = 0
            j = 1
        else:
            i = 1
            j = 0
        templ_len = int(k.split(':')[-2])
        templ_start1 = int(pps[i][0][7])
        templ_end1 = int(pps[i][1][8])
        templ_start2 = int(pps[j][0][7])
        templ_end2 = int(pps[j][1][8])
        
        if templ_end1 - templ_start1 + 1 != templ_len \
        or templ_end2 - templ_start2 + 1 != templ_len:
            continue
        templ1 = chroms[pps[i][0][0]][templ_start1-1:templ_end1].upper()  
        templ2 = chroms[pps[j][0][0]][templ_start2-1:templ_end2].upper()
        if templ1[prim_len:-prim_len] == templ2[prim_len:-prim_len]:
            continue
        align_pairs[k] = pps
        align_templs[k] = (templ1, templ2)
    return align_pairs, align_templs


def filter_hsps(hsps, length):
    filtered_hsps = []
    for e in hsps:
        q_end = int(e[6])
        if q_end <= length - 2:  # experienced rules
            continue
        if q_end > length and q_end <= length*2 + 18:  # experienced rules
            continue
        if '-' in e[-2] or '-' in e[-1]:  # have gap
            continue
        if mismatch_number(e[-2][-5:], e[-1][-5:]) >= 2:  # mismatch in 3'end is more than or equal to 2
            continue
        filtered_hsps.append(e)
    return filtered_hsps   
    
    
def mismatch_number(query, subject):
    counter = 0
    for i in range(len(query)):
        if query[i] != subject[i]:
            counter += 1
    return counter

    
def paired_primers(hsps, length):
    left_hsps = [e for e in hsps if int(e[5]) <= length]
    right_hsps = [e for e in hsps if int(e[5]) > length]
    pps = []
    while len(left_hsps) > 0:
        left_start = int(left_hsps[0][7])  
        for i in range(len(right_hsps)):
            right_end = int(right_hsps[i][8])
            if left_hsps[0][0] == right_hsps[i][0] and abs(right_end - left_start) < 4000:
                pps.append((left_hsps[0], right_hsps[i]))
        del(left_hsps[0])
    return pps    
       

def query_summary(blast_dict, filtered=True):
    blast_counter = []
    for k, v in blast_dict.items():
        length = int(k.split(':')[-1])
        if filtered:
            f_hsps = filter_hsps(v, length)
            pps = paired_primers(f_hsps, length)
        else:
            pps = paired_primers(v, length)
        blast_counter.append((k, len(pps)))
    summary = {}
    for e in blast_counter:
        summary[e[1]] = summary.get(e[1], 0) + 1
    return summary, blast_counter

#%%
def alignment_memo():
    os.chdir('spec_check/')
    
    for file in glob('../pre_proc/*.fa'):
        blast_cmd = "time blastn -db ../ChAPDes/raw_data/blastdb/hg19New -query %s \
        -out %s.out -evalue 30000 -max_hsps 100 -num_threads 16 -outfmt '6 std qseq sseq' \
        -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -perc_identity 65 &> blast.log" % (file, file.split('/')[-1])
        ### change your code according to your path
        with subprocess.Popen(blast_cmd, shell=True) as proc:
            proc.communicate()
        
    os.chdir('../')
#%%
def main():
    alignment_memo()


if __name__ == '__main__':
    main()