import argparse
parser = argparse.ArgumentParser(prog='ld_clump_prs.py', description='''
    This script calculates LD-clumping based PRS using hail.
    It takes a list of sample subsets.
    For each subset, it calculates PRS of one or multiple GWAS. 
''')

## Genotype inputs
parser.add_argument('--bgen-path', required=True, help='''
    Template of bgen file, need to contain {chr_num} in place of chromosome 
    number
''')
all_chrs = ','.join([ str(i) for i in range(1, 23) ])
parser.add_argument('--chrs', default=all_chrs, type=str, help='''
    The chromosome indexes separated by ',' to include. 
    If not specified, all autosomes will be used.
''')
parser.add_argument('--bgen-sample', required=True, help='''
    Sample file of bgen
''')
parser.add_argument('--bgen-index', default=None, help='''
    If not specified, assume the idx2 files are in the same directory
    and same name as bgen files.
    Otherwise, specify with {chr_num} as well
''')

## subset and GWAS list 
parser.add_argument('--subset-and-gwas', required=True, help='''
    A YAML file taking the structure:
    subset1:
        indiv_lists: 'list1.txt,list2.txt,...'
        GWASs:
            gwas1: 
                sum_stat: 'gwas_path1.txt'
                ld_clump: 'clump_path1.txt'
            gwas2: 
                sum_stat: 'gwas_path2.txt'
                ld_clump: 'clump_path2.txt'
    subset2:
        indiv_lists: 'list3.txt,list4.txt,...'
        GWASs:
            gwas1: 
                sum_stat: 'gwas_path3.txt'
                ld_clump: 'clump_path3.txt'
            gwas2: 
                sum_stat: 'gwas_path4.txt'
                ld_clump: 'clump_path4.txt'        
''')
parser.add_argument('--variant-pool', required=True, help='''
    A TSV GWAS summary statistics file (in Neale's lab format).
    The variants listed there are used as the starting pool of variants.
''')

## PRS parameter
parser.add_argument('--pval-thresholds', required=True, help='''
    GWAS p-value cutoffs used to calculate PRS.
    It should be separated by ',' 
''')

## output 
parser.add_argument('--output-prefix', required=True, help='''
    Prefix of output file
''')

## hail log 
parser.add_argument('--hail-log', default=None, help='''
    Path of hail log file.
    The default is args.output_prefix + '.log'
''')


args = parser.parse_args()

import hail as hl
import logging, os, time, sys
import gwas_helper
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
if '{chr_num}' not in args.bgen_path:
    logging.info('Wrong --bgen-path! It should contain {chr_num}! Exit')
    sys.exit()
if '{chr_num}' not in args.bgen_index:
    logging.info('Wrong --bgen-index! It should contain {chr_num}! Exit')
    sys.exit()
    
pval_thresholds = []
try:
    for i in args.pval_thresholds.split(','):
        pval_thresholds.append(float(i))
except:
    print('--args.pval-thresholds {} is not floats separated by ",". Exit!'.format(args.pval_thresholds))
    sys.exit()
    
    
# more on args
if args.bgen_index is None:
    args.bgen_index = args.bgen_path + '.idx2'


# some hail environment logging before run
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')


# initialize hail
logging.info('Initialize hail')
hl.init(log = args.hail_log)


# read in GWAS sum stats and clumped variants
logging.info('Read subset and GWAS / LD-clumping YAML')
myinputs = gwas_helper.read_yaml(args.subset_and_gwas)

## collect clump variant
clump_var_files = []
for i in list(myinputs.keys()):
    for j in list(myinputs[i]['GWASs'].keys()):
        clump_var_files.append(myinputs[i]['GWASs'][j]['ld_clump'])

## load variant loop (limiting to clump variant)
logging.info('--> Start loading variant pool')
ht_var_pool = prs_helper.read_gwas_table_with_varlist(args.variant_pool, clump_var_files, type_dic = {'beta' : hl.tfloat, 'pval' : hl.tfloat})
logging.info('--> Loading variant pool FINISHED!')



# load bgen
logging.info('Start to load genotype files')
bgen_path = args.bgen_path.format(
    chr_num = '{' + args.chrs + '}'
)
bgen_sample = args.bgen_sample
bgen_idx_dict = {}
for i in range(1, 23):
    bgen = args.bgen_path.format(chr_num = i)
    bgen_idx = args.bgen_index.format(chr_num = i)
    bgen_idx_dict[bgen] = bgen_idx
tstart = time.time()
mt = hl.import_bgen(
    path = bgen_path, 
    sample_file = bgen_sample, 
    n_partitions = None, 
    index_file_map = bgen_idx_dict, 
    entry_fields = ['dosage'],
    variants = ht_var_pool
)
tend = time.time()
logging.info('Loading genotype FINISHED! {} seconds elapsed'.format(tend - tstart))


# loop over all subsets
logging.info('Start looping over all subsets')
for subset in list(myinputs.keys()):
    logging.info('--> Working on subset = {}'.format(subset))
    indiv_files = myinputs[subset]['indiv_lists'].split(',')
    ht_indiv = hl.import_table(indiv_files, key = ['f0'], no_header = True, delimiter = ' ')
    ## subset by individual
    logging.info('--> Start subsetting genotype')
    tstart = time.time()
    mt_subset = mt.filter_cols(hl.is_defined(ht_indiv[mt.s]))
    mt_subset = mt_subset.repartition(200)
    mt_subset = mt_subset.cache()
    tend = time.time()
    logging.info('--> Subsetting genotype FINISHED! {} seconds elapsed'.format(tend - tstart))
    for gwas in list(myinputs[subset]['GWASs']):
        logging.info('----> Start subset = {} and gwas = {}'.format(subset, gwas))
        gwas_file = myinputs[subset]['GWASs'][gwas]['sum_stat']
        clump_file = myinputs[subset]['GWASs'][gwas]['ld_clump']
        logging.info('----> Start loading GWAS TSV'.format(subset, gwas))
        tstart = time.time()
        gwas_tsv = prs_helper.read_gwas_table_with_varlist(gwas_file, clump_file, type_dic = {'beta' : hl.tfloat, 'pval' : hl.tfloat})
        tend = time.time()
        logging.info('--> Loading GWAS TSV FINISHED! {} seconds elapsed'.format(tend - tstart))
        ## subset by variant
        logging.info('----> Start subsetting GWAS-specific clumping variant'.format(subset, gwas))
        tstart = time.time()
        mt_this = mt_subset.filter_rows(hl.is_defined(gwas_tsv[mt_subset.locus, mt_subset.alleles]))
        tend = time.time()
        logging.info('--> Subsetting GWAS-specific clumping variant FINISHED! {} seconds elapsed'.format(tend - tstart))
        ## annotate variant with gwas sum stat
        logging.info('----> Start calculating PRS'.format(subset, gwas))
        tstart = time.time()
        mt_this = mt_this.annotate_rows(
            gwas_beta = gwas_tsv[mt_this.locus, mt_this.alleles].beta, 
            gwas_pval = gwas_tsv[mt_this.locus, mt_this.alleles].pval
        )
        prs = {
            'pval_thres_' + str(i) : hl.agg.sum(mt_this.gwas_beta * mt_this.dosage * hl.int(mt_this.gwas_pval < i)) for i in pval_thresholds
        }
        mt_this = mt_this.annotate_cols(**prs)
        tend = time.time()
        logging.info('--> Calculating PRS FINISHED! {} seconds elapsed'.format(tend - tstart))
        ## FIXME: this is temporary! the output format should by tsv.bgz once everything gets settled down 
        logging.info('----> Start writing to disk'.format(subset, gwas))
        tstart = time.time()
        mt_this.write('{prefix}_x_{subset}_x_{gwas}.prs.ht'.format(prefix = args.output_prefix, subset = subset, gwas = gwas))
        mt_this = mt_this.annotate_cols(**prs)
        tend = time.time()
        logging.info('--> Writing to disk FINISHED! {} seconds elapsed'.format(tend - tstart))
        
