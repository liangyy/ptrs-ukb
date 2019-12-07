import argparse
parser = argparse.ArgumentParser(prog='index_bgen_for_hail.py', description='''
    Index BGEN file for HAIL
''')
parser.add_argument('--bgen', help='''
    input BGEN file
''')
parser.add_argument('--index_file', type = int, default = 0, help='''
    output index file path
''')
args = parser.parse_args()

import hail as hl
import logging, os, time, sys

logging.basicConfig(level = logging.INFO, stream = sys.stderr)
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')

bgen_file = args.bgen  
index_file = args.index_file

logging.info('Start indexing {file}'.format(file = bgen_file))
tstart = time.time()
# hl.index_bgen(bgen_file, index_file_map = {bgen_file : index_file})
logging.info('Finished! {time} seconds elapsed'.format(time = time.time() - tstart))
logging.info('Index file saved as {file}'.format(index_file))
