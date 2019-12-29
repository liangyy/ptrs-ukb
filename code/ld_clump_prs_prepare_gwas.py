import argparse
parser = argparse.ArgumentParser(prog='ld_clump_prs_prepare_gwas.py', description='''
    This script prepare GWAS sum stats (subset to LD clumpped variants)
''')

## If using Google Cloud 
parser.add_argument('--google-cloud-project', default=None, help="Add project name if using Google Cloud (default: None)")

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
parser.add_argument('--gwas-ht', default=None, help='''
    It should be hail Table containing nested GWAS results.
    If it is specified, the GWAS sum_stats should be 
    the list name and inner list name of the results 
    inside hail Table.
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
# print(args.dont_overwrite)
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


# some hail environment logging before run
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')


# initialize hail
logging.info('Initialize hail')
hl.init(log = args.hail_log)


# read in GWAS sum stats and clumped variants
logging.info('Read subset and GWAS / LD-clumping YAML')
myinputs = gwas_helper.read_yaml(args.subset_and_gwas, args.google_cloud_project)
logging.info('Read GWAS hail Table')
gwas_ht = hl.read_table(args.gwas_ht)
gwas_ht = gwas_ht.key_by('rsid')
gwas_ht = gwas_ht.repartition(40)
gwas_ht = gwas_ht.cache()


# loop over all subsets
logging.info('Start looping over all subsets')
for subset in list(myinputs.keys()):
    for gwas in list(myinputs[subset]['GWASs']):
        logging.info('----> Start subset = {} and gwas = {}'.format(subset, gwas))
        gwas_file = f'{subset}_x_{gwas}'  # myinputs[subset]['GWASs'][gwas]['sum_stat']
        clump_file = myinputs[subset]['GWASs'][gwas]['ld_clump']
        logging.info('----> Start working GWAS TSV'.format(subset, gwas))
        tstart = time.time()
        gwas_tmp_ht = args.output_prefix + '_x_checkpoint_x_' + subset + '_x_' + gwas + '.ht'
        prs_helper.read_gwas_table_with_varlist(gwas_file, clump_file, type_dic = {'beta' : hl.tfloat, 'pval' : hl.tfloat}, checkpoint_path = gwas_tmp_ht, gwas_ht = gwas_ht, no_return = True)
        tend = time.time()
        logging.info('----> Working GWAS TSV FINISHED! {} seconds elapsed'.format(tend - tstart))
        
        
