import os 
import datetime
from pharmvip_guideline import pack_path
path_dir = pack_path + '/../resources/allele_definitions/'

now = datetime.datetime.now()
value = [now.year,now.month,now.day]
string_value = [str(int) for int in value]
now_format = '_'.join(string_value)
path_dir = path_dir + now_format + '/table'
try:
    os.makedirs(path_dir)
except:
    print('directory already exist')

gene_set = {'CACNA1S', 'CFTR', 'CYP2B6', 
    'CYP2C9', 'CYP2C19', 'CYP2D6' , 'CYP3A5',
    'CYP4F2', 'DPYD', 'G6PD', 'IFNL3', 'NUDT15',
    'RYR1', 'SLCO1B1', 'TPMT', 'UGT1A1', 'VKORC1'}
    
os.chdir(path_dir)

def update_allele_table():
    for gene in gene_set:

        path = f'https://api.pharmgkb.org/v1/download/file/attachment/{gene}_allele_definition_table.xlsx'
        os.system(f'wget {path}')

if __name__ == "__main__":   
    update_allele_table()
