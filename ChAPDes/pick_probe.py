# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 14:14:28 2018

@author: Link
"""

#%%
from primer3 import *
from spec_check import *
import csv
#%%


def pick_probe_task(align_templs, task, **kwargs):
    param_lst = primer3_params(kwargs)
    templ_lst = templ_seqs(align_templs, kwargs['primer_internal_min_size'])
    path = 'pick_probe/'+task+'/pick_probe_task.bio'
    with open(path, 'w') as fh:
        for templ in templ_lst:
            for e in templ:
                fh.writelines(e+'\n')
            for e in param_lst:
                fh.writelines(e+'\n')
    primer3_cmd = 'primer3_core'
    subprocess.call(primer3_cmd, stdin=open(path), stdout=open('pick_probe/'+task+'/pick_probe_out.bio', 'w'))


def templ_seqs(align_templs, probe_len):
    templ_lst = []
    for k, v in align_templs.items():
        templ = []
        templ.append('SEQUENCE_ID='+k)
        templ.append('SEQUENCE_TEMPLATE='+v[0])
        templ.append('PRIMER_TASK=pick_hyb_probe_only')
        mm_sites = mismatch_sites(v)
        prim_len = int(k.split(':')[-1])
        templ_len = int(k.split(':')[-2])
        regions = excluded_regions(mm_sites, prim_len, probe_len, templ_len)
        region_str = ''
        for t in regions:
            region_str = region_str + str(t[0]) + ',' + str(t[1]) + ' '
        templ.append('SEQUENCE_INTERNAL_EXCLUDED_REGION='+region_str)
        templ_lst.append(templ)
    return templ_lst
    
    
def mismatch_sites(templs):
    mm_sites = []
    for i in range(len(templs[0])):
        if templs[0][i] != templs[1][i]:
            mm_sites.append(i)
    return mm_sites


def excluded_regions(mm_sites, prim_len, probe_len, templ_len):
    regions = []
    l_start = mm_sites[0] - (probe_len - 1)
    if l_start <= prim_len - 1:  
        regions.append((0, prim_len))
    else:
        regions.append((0, l_start))  
    for i in range(1, len(mm_sites)):
        r_start = mm_sites[i] - (probe_len - 1)
        l_end = mm_sites[i-1] + probe_len - 1
        if r_start > l_end + 1:  
            regions.append((l_end+1, r_start-(l_end+1))) 
    r_end = mm_sites[-1] + probe_len - 1
    if r_end + 1 >= templ_len - prim_len:
        regions.append((templ_len-prim_len, prim_len))
    else:
        regions.append((r_end+1, templ_len-(r_end+1)))  
    return regions


def check_all_task(align_pairs, task, **kwargs):
    templ_dict = read_bio('pick_probe/'+task+'/pick_probe_out.bio')
    param_lst = primer3_params(kwargs)
    seqs_lst = all_seqs(templ_dict, align_pairs)
    path = 'pick_probe/'+task+'/check_all_task.bio'
    with open(path, 'w') as fh:
        for seqs in seqs_lst:
            for e in seqs:
                fh.writelines(e+'\n')
            for e in param_lst:
                fh.writelines(e+'\n')
    primer3_cmd = 'primer3_core'
    subprocess.call(primer3_cmd, stdin=open(path), stdout=open('pick_probe/'+task+'/check_all_out.bio', 'w'))
    

def all_seqs(templ_dict, align_pairs):
    seqs_lst = []
    for k, v in templ_dict.items():
        for i in range(int(v['PRIMER_INTERNAL_NUM_RETURNED'])):
            seqs = []
            seqs.append('SEQUENCE_ID='+k+':'+str(i))
            seqs.append('SEQUENCE_TEMPLATE='+v['SEQUENCE_TEMPLATE'].upper())
            oligo_name = 'PRIMER_INTERNAL_'+str(i)+'_SEQUENCE'
            seqs.append('SEQUENCE_INTERNAL_OLIGO'+'='+v[oligo_name].upper())
            for j in range(2):
                if not align_pairs[k][j][0][0] == k.split(':')[0]:
                    continue
                seqs.append('SEQUENCE_PRIMER='+align_pairs[k][j][0][11])
                seqs.append('SEQUENCE_PRIMER_REVCOMP='+reverse_complement(align_pairs[k][j][1][11]))
            seqs.append('PRIMER_TASK=check_primers')
            seqs.append('PRIMER_INTERNAL_MAX_SIZE=35')
            seqs_lst.append(seqs)
    return seqs_lst

    
#%%
def process(task):
    print('Reading blast file...')
    blast_lst = read_tsv('spec_check/'+task+'/blast.out')
    blast_dict = lst_to_dict(blast_lst)
    print('Filtering blast result...')
    align_pairs, align_templs = spec_check(blast_dict)
    print('Checking primers and probes...')
    check_all_task(align_pairs, task, primer_min_tm=57.0, primer_opt_tm=60.0, primer_max_tm=63.0, 
             primer_pair_max_diff_tm=3.0, primer_min_gc=20.0, primer_opt_gc_percent=50.0,
             primer_max_gc=80.0, primer_internal_min_size=30, primer_internal_opt_size=33, 
             primer_internal_max_size=35, primer_internal_min_tm=67.0, primer_internal_opt_tm=70.0,
             primer_internal_max_tm=73.0, primer_internal_min_gc=35.0, primer_internal_opt_gc_percent=50.0,
             primer_internal_max_gc=65.0)

    
def main():
    process(task)
    
    
    
if __name__ == '__main__':
    main()