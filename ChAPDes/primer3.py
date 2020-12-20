# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 22:58:41 2018

@author: Link
"""


#%%
import subprocess
from datetime import datetime
from pre_proc import *
#%%


def primer_check_task(task, **kwargs):
    """Run Primer3 based on parameter and primer pairs. """
    param_lst = primer3_params(kwargs)
    left_primer_file = 'pre_proc/'+task+'/primer_pairs.1.fa'
    right_primer_file = 'pre_proc/'+task+'/primer_pairs.2.fa'
    seq_lst = primer_seqs(left_primer_file, right_primer_file)
    for i in range(0, len(seq_lst), 999999):  # step length
        path = 'primer3_out/'+task+'/primer_check_task.'+str(i)+'.bio'
        with open(path, 'w') as fh:
            for primer in seq_lst[i:i+999999]:
                for e in primer:
                    fh.writelines(e+'\n')
                for e in param_lst:
                    fh.writelines(e+'\n')
        primer3_cmd = 'primer3_core'
        subprocess.call(primer3_cmd, stdin=open(path), stdout=open('primer3_out/'+task+'/primer_check_out.bio', 'a'))


def primer3_params(param_dict):
    """Load parameters into the format of Primer3."""
    param_lst = []
    for k, v in param_dict.items():
        param_lst.append(k.upper()+'='+str(v))
    param_lst.append('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=../primer3-2.4.0/src/primer3_config/')## change your code according to your path
    param_lst.append('=')
    return param_lst


def primer_seqs(left_primer_file, right_primer_file):
  """Load primer pairs sequences into the format of Primer3."""
    seqs = {}
    with open(left_primer_file) as fh:
        while True:
            header = fh.readline().rstrip()
            seq = fh.readline().rstrip()
            if len(header) == 0:
                break
            seqs[header] = [seq]
    with open(right_primer_file) as fh:
        while True:
            header = fh.readline().rstrip()
            seq = fh.readline().rstrip()
            if len(header) == 0:
                break
            seqs[header].append(seq)
    seq_lst = []
    for k, v in seqs.items():
        primer = []
        primer.append('SEQUENCE_ID='+k[1:])
        primer.append('SEQUENCE_PRIMER='+v[0])
        primer.append('SEQUENCE_PRIMER_REVCOMP='+v[1])
        primer.append('PRIMER_TASK=check_primers')
        seq_lst.append(primer)
    return seq_lst
    

def read_bio(path):
    """Read .bio file."""
    primer_dict = {}
    with open(path) as fh:
        seq_info = {}
        while True:
            line = fh.readline().rstrip()
            if len(line) == 0:
                break
            elif line == '=':
                primer_dict[seq_info['SEQUENCE_ID']] = seq_info
                seq_info = {}
            else:
                k, v = line.split('=')
                seq_info[k] = v
    return primer_dict


def primer_dict_to_pps(primer_dict, rc=False):
    """Get sequences of primers and reversed primers."""
    pps = {}
    for k, v in primer_dict.items():
        if int(v['PRIMER_PAIR_NUM_RETURNED']) != 0:
            if rc:
                pps[k] = (v['SEQUENCE_PRIMER'], v['SEQUENCE_PRIMER_REVCOMP'])
            else:
                pps[k] = (v['SEQUENCE_PRIMER'], reverse_complement(v['SEQUENCE_PRIMER_REVCOMP']))
    return pps



def main():
    primer_check_task(primer_min_tm=57.0, primer_opt_tm=60.0, primer_max_tm=63.0, 
                primer_pair_max_diff_tm=3.0, primer_min_gc=20.0, primer_opt_gc_percent=50.0,
                primer_max_gc=80.0)
    
    
if __name__ == '__main__':
    main()