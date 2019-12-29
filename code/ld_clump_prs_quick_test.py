import argparse
parser = argparse.ArgumentParser(prog='ld_clump_prs_quick_test.py', description='''
    This script calculates LD-clumping based PRS using hail.
''')


## Genotype inputs
parser.add_argument('--mt-input', default=None, help='''
    Input of pre-subsetted genotype
''')

## subset and GWAS list 
parser.add_argument('--gwas-tsv', required=True, help='''
    GWAS table (LD clump variants)       
''')
parser.add_argument('--trait-names', default=None, help='''
    List of trait name, separaed by ','
''')

## PRS parameter
parser.add_argument('--pval-thresholds', default='5e-8,0.05', help='''
    GWAS p-value cutoffs used to calculate PRS.
    It should be separated by ',' 
''')

## output 
parser.add_argument('--output', required=True, help='''
    Output file
''')

## hail log 
parser.add_argument('--hail-log', default=None, help='''
    Path of hail log file.
    The default is args.output_prefix + '.log'
''')


args = parser.parse_args()
# print(args.dont_overwrite)
import hail as hl
import logging, os, time, sys
import prs_helper


# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# input sanity check
if args.hail_log is None:
    args.hail_log = args.output_prefix + '.log'
    
pval_thresholds = []
try:
    for i in args.pval_thresholds.split(','):
        pval_thresholds.append(float(i))
except:
    print('--args.pval-thresholds {} is not floats separated by ",". Exit!'.format(args.pval_thresholds))
    sys.exit()   
    


# some hail environment logging before run
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')


# initialize hail
logging.info('Initialize hail')
hl.init(log = args.hail_log)

# load genotype
logging.info('Start reading genotype file')
mt_sub_name = args.mt_input
mt_subset = hl.read_matrix_table(mt_sub_name)


# load gwas table
logging.info('Start read GWAS table')
columns = [ 'variant' ]
# columns = []
for trait in args.trait_names.split(','):
    columns.append('beta.' + trait)
    columns.append('pval.' + trait)
type_dic = {
    i : hl.tfloat for i in columns
}
type_dic['variant'] = hl.tstr
gwas_tsv = prs_helper.read_gwas_table(args.gwas_tsv, type_dic = type_dic)
logging.info('--> Start caching GWAS table')
gwas_tsv = gwas_tsv.repartition(40)
gwas_tsv = gwas_tsv.cache()


# prs calculation
logging.info('Start PRS calculation')
# mt_subset = mt_subset.annotate_rows(**annot_gwas)
logging.info('--> Start annotating genotype with GWAS information')
annot_gwas = {}
for trait in args.trait_names.split(','):
    annot_gwas[f'beta_{trait}'] = gwas_tsv[mt_subset.locus, mt_subset.alleles][f'beta.{trait}']
    annot_gwas[f'pval_{trait}'] = gwas_tsv[mt_subset.locus, mt_subset.alleles][f'pval.{trait}']
mt_subset = mt_subset.annotate_rows(**annot_gwas)
logging.info('--> Start PRS annotation')
prs = {}
for trait in args.trait_names.split(','):
    for pval in pval_thresholds:
        prs[f'{trait}_pval_thres_{pval}'] = hl.agg.sum(mt_subset[f'beta_{trait}'] * mt_subset.dosage * hl.int(mt_subset[f'pval_{trait}'] < pval))
mt_subset = mt_subset.annotate_cols(**prs)
logging.info('--> Start export')
mt_subset.col.export(args.output)
        
