import yaml

def read_yaml(filename):
    with open(filename, 'r') as f:
        mydic = yaml.safe_load(f)
    return mydic
def get_variant_qc_filter(qc_table, name, value):
    if name is 'maf':
        return [qc_table['AF'][0] > value, qc_table['AF'][1] > value]
    if name is 'hwe':
        return [qc_table['p_value_hwe'] > value]
