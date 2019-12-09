import argparse
parser = argparse.ArgumentParser(prog='extract_var_list_from_predictdb.py', description='''
    Extract the variant ID's as a one column list from predictdb
''')
parser.add_argument('--predictdb', help='''
    input predictdb file in DB format
''')
parser.add_argument('--varcol', help='''
    column name of variant ID's you'd like to extract.
    If you are not sure what to put here, you can use 'SHOW'
    and the script to print all available columns in predictdb
''')
parser.add_argument('--output', help='''
    output file name
''')
args = parser.parse_args()

import sqlite3 
import logging, time, sys
from tabulate import tabulate


# configing util
logging.basicConfig(level = logging.INFO, stream = sys.stderr)

conn = sqlite3.connect(args.predictdb)
cur = conn.cursor()

logging.info('Showing tables')
cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
tab_nested = cur.fetchall()
tab = [ i[0] for i in tab_nested ]
weight_header = None
for i in tab:
    sql = 'PRAGMA table_info({})'.format(i)
    print('schema for {}'.format(i))
    tmp = cur.execute(sql).fetchall()
    print(tabulate(tmp), file = sys.stderr)
    if i == 'weights':
        weight_header = [ i[1] for i in tmp ]

# check if weights is table in predictdb
if 'weights' not in tab:
    logging.info('`weights` is not a table in predictdb! Exit')
    sys.exit()

# work with weights table
logging.info('--varcol is {}'.format(args.varcol))
if args.varcol == 'SHOW':
    logging.info('Show first 3 rows in weights and exit')
    sql_weights = 'select * from weights limit 3'
    content = cur.execute(sql_weights).fetchall()
    print(tabulate(
        content,
        headers = weight_header
    ), file = sys.stderr)
elif args.varcol not in weight_header:
    logging.info('Wrong --varcol is specified! Exit')
    sys.exit()
else:
    sql_weights = 'select distinct({}) from weights'.format(args.varcol)
    
    logging.info('Start to query --varcol from weights using SQL:')
    print(sql_weights, file = sys.stderr)
    tstart = time.time()
    content = cur.execute(sql_weights).fetchall()
    tend = time.time()
    logging.info('Query finished! {} seconds elapsed'.format(tend - tstart))
    
    logging.info('Start to write to disk')
    tstart = time.time()
    f = open(args.output, 'w')
    for i in content:
        f.write(i[0] + '\n')
    tend = time.time()
    logging.info('Query finished! {} seconds elapsed'.format(tend - tstart))
            
    
