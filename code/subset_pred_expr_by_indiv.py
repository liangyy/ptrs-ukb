import argparse
parser = argparse.ArgumentParser(prog='subset_pred_expr_by_indiv.py', description='''
    Input1: predicted expression matrix generated from predixcan_prediction.py
    Input2: a list of individual ID
    Output: the predicted expression matrix (as TSV) for the subset of individuals.
    Also, if a gene as no variation, it will be removed.
''')

parser.add_argument('--pred-expr', required=True, help='''
    Predicted expression in HDF5 format (generated by predixcan_prediction.py)
''')
parser.add_argument('--indiv-list', required=True, help='''
    The list of individuals (it can have several columns but the first one 
    will be treated as individual ID)
''')
parser.add_argument('--output', required=True, help='''
    Predicted expression on subset individuals in TSV format
''')

args = parser.parse_args()

import pandas as pd
import numpy as np
import h5py
import logging, os, time, sys
import my_hail_helper as myhelper
import gcta_helper as ghelper


# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# read individual list
logging.info('Read individual list')
indiv_list = myhelper.read_indiv_list(args.indiv_list)

# load predicted expression
logging.info('Load predicted expression')
f = h5py.File(args.pred_expr, 'r')

# get individual and gene ids from predicted expression
logging.info('Extracting individual IDs from predicted expression and individual list')
col_names = f['samples'][:]
row_names = f['genes'][:]
target_names = np.array(indiv_list).astype('S10')  # to match the string type to the one saved in hdf5
logging.info('OK, it\'s done!')

# main
logging.info('Main function starts')
tstart = time.time()
sub_pred_expr, sub_indivs = ghelper.extract_cols(f['pred_expr'], col_names, target_names)
tend = time.time()
logging.info('Main function FINISHED! {} seconds elapsed'.format(tend - tstart))

# remove genes with no variation
logging.info('Removing genes with constant predicted expression')
sub_pred_expr, sub_genes = ghelper.remove_constant_row(sub_pred_expr, row_names)

# convert subset predicted expression to pandas data.frame
logging.info('Convert predicted expression to data.frame')
tstart = time.time()
df = pd.DataFrame(sub_pred_expr)
df.columns = sub_indivs.astype('str')
df.index = sub_genes.astype('str')
df['gene'] = df.index
df = df.reset_index(drop = True)
tend = time.time()
logging.info('Conversion FINISHED! {} seconds elapsed'.format(tend - tstart))

# write to disk
logging.info('Writing to disk as TSV file')
tstart = time.time()
df.to_csv(args.output, index = None, header = True, sep = '\t')
tend = time.time()
logging.info('Writing to disk FINISHED! {} seconds elapsed'.format(tend - tstart))


