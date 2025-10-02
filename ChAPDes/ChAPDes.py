from pre_proc import *
from primer3 import *
from spec_check import *
from pick_probe import *
from datetime import datetime


def pipeline(chrom, popul, primer_type='a', chrom_start=None, chrom_end=None, templ_min_len=100, templ_max_len=120,
             primer_min_len=20, primer_max_len=25, primer_min_tm=57.0, primer_opt_tm=60.0, primer_max_tm=63.0, 
             primer_pair_max_diff_tm=3.0, primer_min_gc=20.0, primer_opt_gc_percent=50.0,
             primer_max_gc=80.0, primer_internal_min_size=30, primer_internal_opt_size=33, 
             primer_internal_max_size=35, primer_internal_min_tm=67.0, primer_internal_opt_tm=70.0,
             primer_internal_max_tm=73.0, primer_internal_min_gc=35.0, primer_internal_opt_gc_percent=50.0,
             primer_internal_max_gc=65.0):
             
"""   Generated the primer pairs which discarded the TR and SNP,
      checked the specificity by Bowtie2 and BLAST tool and designed the probes step by step.""" 
             
    assert templ_min_len >= 50 and templ_max_len <= 300 and templ_min_len <= templ_max_len, \
    'Template length must be >=50 and <=300!'
    assert primer_min_len >= 18 and primer_max_len <=30 and primer_min_len <= primer_max_len, \
    'Primer length must be >=18 and <=30!'
    assert primer_min_tm >= 0 and primer_max_tm >= 0, \
    'Primer tm must be greater than 0!'
    assert primer_opt_tm >= primer_min_tm and primer_opt_tm <= primer_max_tm, \
    'Primer optimal tm must be between primer min tm and max tm!'
    assert primer_pair_max_diff_tm >= 0, \
    'Primer pair different length must be greater than 0!'
    assert primer_min_gc >= 0 and primer_min_gc <= 100 and primer_max_gc >= 0 and primer_max_gc <= 100, \
    'Primer GC% must be >=0 and <=100!'
    assert primer_opt_gc_percent >= primer_min_gc and primer_opt_gc_percent <= primer_max_gc, \
    'Primer optimal GC% must be between primer min GC% and max GC%!'
    assert primer_internal_min_size >= 0 and primer_internal_max_size >= 0, \
    'Primer internal size must be greater than 0!'
    assert primer_internal_opt_size >= primer_internal_min_size and primer_internal_opt_size <= primer_internal_max_size, \
    'Primer optimal internal size must be between primer min size and max size!'
    assert primer_internal_min_tm >= 0 and primer_internal_max_tm >= 0, \
    'Primer internal tm must be greater than 0!'
    assert primer_internal_opt_tm >= primer_internal_min_tm and primer_internal_opt_tm <= primer_internal_max_tm, \
    'Primer optimal internal tm must be between primer min tm and max tm!'
    assert primer_internal_min_gc >= 0 and primer_internal_min_gc < 100 and primer_internal_max_gc >= 0 and primer_internal_max_gc <= 100, \
    'Primer internal GC% must be >=0 and <=100!'
    assert primer_internal_opt_gc_percent >= primer_internal_min_gc and primer_internal_opt_gc_percent <= primer_internal_max_gc, \
    'Primer optimal internal GC% must be between primer min GC% and max GC%!'
    
    print('Generating primer pairs...')
    templ_len = '%d-%d' % (templ_min_len, templ_max_len)
    primer_len = '%d-%d' % (primer_min_len, primer_max_len)
    task = chrom+'_'+popul+'_'+templ_len+'_'+primer_len+'_'+datetime.strftime(datetime.now(), "%Y-%m-%d_%H:%M:%S")
    for path in ['pre_proc/', 'primer3_out/', 'spec_check/', 'pick_probe/']:
        os.mkdir(path+task)
    sd_lst = read_sd_csv('raw_data/build37.csv')
    times = sequence_overlapped_times(sd_lst, chrom)
    seq = read_chromosome('raw_data/hg19/'+chrom+'.fa')
    if chrom_start is None:
        chrom_start = 1
    if chrom_end is None:
        chrom_end = len(seq) + 1
    seqs = extract_non_overlapped_sequences(seq, chrom_start, chrom_end, times, chrom, templ_len=templ_max_len)
    rcseqs = extract_non_overlapped_sequences(seq, chrom_start, chrom_end, times, chrom, templ_len=templ_max_len, rc=True)
    del(sd_lst)
    del(times)
    del(seq)
    pps = fast_primer_pairs(seqs, rcseqs, min_k=templ_min_len, max_k=templ_max_len, min_len=primer_min_len, max_len=primer_max_len)
    rc_pps = fast_primer_pairs(seqs, rcseqs, min_k=templ_min_len, max_k=templ_max_len, min_len=primer_min_len, max_len=primer_max_len, rc=True)
    del(seqs)
    del(rcseqs)
    
    print('Filtering SNPs and tandem repeats...')
    f_pps = filter_primer_tandem_repeat(pps, chrom)   
    del(pps)
    f_pps = filter_primer_snp(f_pps, popul, chrom)
    f_rc_pps = filter_primer_tandem_repeat(rc_pps, chrom)
    del(rc_pps)
    f_rc_pps = filter_primer_snp(f_rc_pps, popul, chrom)
    
    print('Writing primer pairs...')
    write_primer_pairs(f_pps, 'pre_proc/'+task+'/primer')
    del(f_pps)
    write_primer_pairs(f_rc_pps, 'pre_proc/'+task+'/primer_pairs')
    del(f_rc_pps)
    
    print('Filtering primer pairs...')
    primer_check_task(task, primer_min_tm=primer_min_tm, primer_opt_tm=primer_opt_tm, primer_max_tm=primer_max_tm, 
                      primer_pair_max_diff_tm=primer_pair_max_diff_tm, primer_min_gc=primer_min_gc, 
                      primer_opt_gc_percent=primer_opt_gc_percent, primer_max_gc=primer_max_gc)
    primer_dict = read_bio('primer3_out/'+task+'/primer_check_out.bio')
    pps = primer_dict_to_pps(primer_dict)
    write_primer_pairs(pps, 'spec_check/'+task+'/primer')
    del(pps)
    rc_pps = primer_dict_to_pps(primer_dict, rc=True)
    write_primer_pairs(rc_pps, 'spec_check/'+task+'/primer_pairs')
    del(rc_pps)
    del(primer_dict)
    
    print('Checking specificity...')
    alignment(task, primer_type)
    
    print('Reading blast file...')
    blast_lst = read_tsv('spec_check/'+task+'/blast.out')
    blast_dict = lst_to_dict(blast_lst)
    del(blast_lst)
    
    print('Filtering blast result...')
    align_pairs, align_templs = spec_check(blast_dict, primer_type)
    del(blast_dict)
    
    print('Picking probes...')
    pick_probe_task(align_templs, task, primer_internal_min_size=primer_internal_min_size, 
                    primer_internal_opt_size=primer_internal_opt_size, primer_internal_max_size=primer_internal_max_size, 
                    primer_internal_min_tm=primer_internal_min_tm, primer_internal_opt_tm=primer_internal_opt_tm,
                    primer_internal_max_tm=primer_internal_max_tm, primer_internal_min_gc=primer_internal_min_gc, 
                    primer_internal_opt_gc_percent=primer_internal_opt_gc_percent,
                    primer_internal_max_gc=primer_internal_max_gc)
    del(align_templs)
    
    print('Checking primers and probes...')
    check_all_task(align_pairs, task)
    
def main():
"""    main function of ChAPDes """
    pipeline(chrom='chr21', popul='EAS') # an example to design the probes on chromosome 21 of Asian


def ml_integration_example():
    """
    Example of integrating machine learning for enhanced primer design decisions.
    This demonstrates how to use the imbalanced classification toolkit with ChAPDes results.
    """
    try:
        from imbalanced_classifier import ImbalancedBinaryClassifier
        from genomic_imbalanced_example import create_genomic_like_dataset

        print("Machine Learning Integration Example")
        print("=" * 50)

        # Create example genomic data (normally this would come from primer design results)
        X, y, feature_names = create_genomic_like_dataset(
            n_samples=1000, n_features=20, imbalance_ratio=0.2
        )

        print(f"Example dataset: {X.shape[0]} samples, {X.shape[1]} features")
        print(f"Class distribution: {y.value_counts().to_dict()}")

        # Train a classifier on primer/genomic features
        clf = ImbalancedBinaryClassifier(method='smote')
        clf.fit(X, y, validation_split=0.2)

        # Example prediction (normally on new primer candidates)
        test_predictions = clf.predict(X[:10])
        test_probabilities = clf.predict_proba(X[:10])

        print(f"\nExample predictions: {test_predictions}")
        print(f"Prediction probabilities shape: {test_probabilities.shape}")

        # Evaluate performance
        metrics = clf.evaluate()
        print(f"\nModel F1-Score: {metrics['f1_score']".4f"}")
        print(f"ROC AUC: {metrics['roc_auc']".4f"}")

        print("\nML integration completed successfully!")
        print("You can now use this classifier to predict primer performance,")
        print("classify genomic variants, or enhance other ChAPDes decisions.")

    except ImportError as e:
        print(f"ML modules not available: {e}")
        print("To use ML features, install required packages:")
        print("pip install -r requirements_ml.txt")


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--ml-example':
        ml_integration_example()
    else:
        main() 