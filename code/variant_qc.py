import argparse
parser = argparse.ArgumentParser(prog='variant_qc.py', description='''
    This script reads in genotype files in bgen
    and calls variant_qc in hail.
    The resulting file is saved in hail MatrixTable format.
''')
parser.add_argument('--bgen-path', required=True, help='''
    Template of bgen file, need to contain {chr_num} in place of chromosome number
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
parser.add_argument('--hail-block-size', type=int, default=128, help='''
    The block_size in hail.import_bgen
''')
parser.add_argument('--output-filename', required=True, help='''
    Filename of output. 
    If the directory does not exist, it will be created
''')
args = parser.parse_args()

import hail as hl
import logging, time, sys, os

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# input sanity check
if '{chr_num}' not in args.bgen_path:
    logging.info('Wrong --bgen-path! It should contain {chr_num}! Exit')
    sys.exit()
if '{chr_num}' not in args.bgen_index:
    logging.info('Wrong --bgen-index! It should contain {chr_num}! Exit')
    sys.exit()

# more on args
if args.bgen_index is None:
    args.bgen_index = args.bgen_path + '.idx2'

# init export pyspark setup
logging.info('Echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')

# init hail
logging.info('Initialize hail')
hl.init()

# load genotype files
logging.info('Start to load genotype files with block_size = {}'.format(args.hail_block_size))
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
    entry_fields = ['GT'],  # to make variant QC work, I can only use GT here 
    block_size = args.hail_block_size
)
tend = time.time()
logging.info('Loading genotype finished! {} seconds elapsed'.format(tend - tstart))

# show n_partitions
logging.info('n_partitions of the loaded data is {}'.format(mt.n_partitions()))

# do variant QC
logging.info('Perform variant QC')
tstart = time.time()
mt = hl.variant_qc(mt)
tend = time.time()
logging.info('hail.variant_qc finished! {} seconds elapsed'.format(tend - tstart))

# and write it to disk
logging.info('Export variant QC')
target_folder = os.path.dirname(args.output_filename)
if not os.path.exists(target_folder) and target_folder is not '':
    os.makedirs(target_folder)
filename = args.output_filename + '.mt'
tstart = time.time()
# mt.variant_qc.export(filename)
tend = time.time()
logging.info('hail.export finished! {} seconds elapsed'.format(tend - tstart))

