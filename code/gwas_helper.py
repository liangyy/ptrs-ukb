import yaml, re, gcsfs

def read_yaml(filename, gcp_project = None):
    if gcp_project is None:
        with open(filename, 'r') as f:
            mydic = yaml.safe_load(f)
    else:
        fs = gcsfs.GCSFileSystem(project = gcp_project)
        filename = re.sub('gs://', '', filename)
        with fs.open(filename, 'r') as f:
            mydic = yaml.safe_load(f)
    return mydic

def get_variant_qc_filter(qc_table, name, value):
    if name == 'maf':
        return [qc_table['AF'][0] > value, qc_table['AF'][1] > value]
    if name == 'hwe':
        return [qc_table['p_value_hwe'] > value]
def apply_variant_qc_filter(qc_ht, name, value):
    if name == 'maf':
        value = float(value)
        qc_ht = qc_ht.filter(qc_ht.variant_qc.AF[0] > value)
        qc_ht = qc_ht.filter(qc_ht.variant_qc.AF[1] > value)
        return qc_ht
    if name == 'hwe':
        value = float(value)
        qc_ht = qc_ht.filter(qc_ht.variant_qc.p_value_hwe > value)
        return qc_ht
    return qc_ht
