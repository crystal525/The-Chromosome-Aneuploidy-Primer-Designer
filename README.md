# The Chromosome Aneuploidy Primer Designer (ChAPDes)



## Introduction

#### Implementation of ChAPDes object detection pipeline using Python. The ChAPDes is a primers and probes design method for aneuploidy detection. This method introduced Bowtie2 tool and improved present Primer-BLAST tool that designs target-specific PCR primers. Compared with the Primer-BLAST tool, this method sped up the primer specificity detection greatly.



##   Description

#### The ChAPDes consists of 4 modules, and performs its evaluation in the following order: data preprocessing, primer design, primer specificity checking and probe design. 



## How to Use?

#### Create new folders named raw_data, pre_proc, primer3_out, spec_check and pick_probe in the current environment. Make sure that all files are present in the corresponding directory. 

- #### By default,it runs as follows.

```python
$python ChAPDes.py
```



## Arguments

- ### Data Preprocessing Module

  - #### chrom :  Number of chromosome for aneuploidy detiction. 

  - #### popul :  Population for detection. Possible types are 

    - ##### "EAS" :   East Asian

    - ##### "SAS" :  South Asian

    - ##### "EUR" :  European

    - ##### "AFR" : African

      	- ##### "AMR" : Ad Mixed American

  - #### chrom_start : Start site of chosen chromosome.

  - #### chrom_end : End site of chosen chromosome. 

  - #### templ_min_len : Minimum length of templates.

  - #### templ_max_len : Maxmum length of templates.

  - #### primer_min_len : Minimum length of primers.

  - #### primer_max_len : Maxmum length of primers.

  - #### primer_min_tm :  Minimum acceptable melting temperature(Celsius) for a primer oligo.

  - #### primer_opt_tm : Optimum melting temperature(Celsius) for a primer oligo.

  - #### primer_max_tm : Maximum acceptable melting temperature(Celsius) for a primer oligo.

  - #### primer_pair_max_diff_tm :  Maximum acceptable (unsigned) difference between the melting temperatures of the left and right primers.

  - #### primer_min_gc :  Minimum allowable percentage of Gs and Cs in any primer.

  - #### primer_opt_gc_percent : Optimum GC percent. 

  - #### primer_max_gc : Maxmum allowable percentage of Gs and Cs in any primer.

  - #### primer_internal_min_size : Minimum acceptable length of a primer. Must be greater than 0.

  - #### primer_internal_opt_size : Optimum length (in bases) of a primer.

  - #### primer_internal_max_size : Maximum acceptable length (in bases) of a primer. This parameter cannot be larger than 35.

  - #### primer_internal_min_tm : Equivalent parameter of primer_min_tm for the internal oligo.

  - #### primer_internal_opt_tm : Equivalent parameter of primer_opt_tm for the internal oligo.

  - #### primer_internal_max_tm : Equivalent parameter of primer_max_tm for the internal oligo.

  - #### primer_internal_min_gc : Equivalent parameter of primer_min_gc for the internal oligo.

  - #### primer_internal_opt_gc_percent : Equivalent parameter of primer_opt_gc_percent for the internal oligo.

  - #### primer_internal_max_gc : Equivalent parameter of primer_max_gc for the internal oligo.





## TODO List

- [ ]  Complete the Linux environment

- [ ]  Prepare the [Segmental Duplication file](http://humanparalogy.gs.washington.edu/build37/build37.htm) in the folder named "/raw_data"

- [ ]  Prepare the [chromosome sequences ](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/)in the folder named"/raw_data/hg19"

- [ ]  Prepare the [index of human genome reference sequence](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for Bowtie2 in the folder named "/raw_data/bowtie2_index"

- [ ]  Prepare the message of [SNP](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) in the folder named "/raw_data/1000genomes"

- [ ]  Create new folders named tandem_repeats and snp_info in "/raw_data".  Run the code as follows. Users could get infomation of TR and SNP.

```python
from pre_proc import *
extract_tandem_repeats(path='raw_data/hg19/')
extract_snp_info(path='raw_data/1000genomes/', popul='EAS', rf=0.1) #obtain SNP infomation of East Asian, users can change parameters to get result.

```

- [ ]  Change your code in the program named "primer3.py" and "spec_check.py" according to your real path.



## Requirements

 - #### Linux 

 - #### Python 3

 - #### Primer3

 - #### BLAST 2.3.0+

 - #### Bowtie2 2.3.3.1
