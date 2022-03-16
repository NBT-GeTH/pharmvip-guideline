import pandas as pd
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))

def  default_gen():
    '''
    emaxple input data to gen hla dipolotype.tsv
    '''
    temp =[
        {
        'sampleid': "tester",
        'gene': "HLA-B",
        'guide_diplotype' : "['HLA-B*15:02:01_negative/HLA-B*15:02:01_negative,HLA-B*57:01:01_positive/HLA-B*57:01:01_positive,HLA-B*58:01_negative/HLA-B*58:01_negative', 'HLA-B*15:02:01_negative/HLA-B*15:02:01_negative,HLA-B*57:01:01_positive/HLA-B*57:01:01_positive,HLA-B*58:01_positive/HLA-B*58:01_positive']",	
        'print_diplotype' : "['Other/Other,*57:01:01/Other,Other/Other', 'Other/Other,*57:01:01/Other,*58:01/Other']",	
        'tool' : "['ATHLATES', 'HLA-HD,KOURAMI']"
        },{
        'sampleid': "tester",
        'gene': "HLA-B",
        'guide_diplotype' : "['HLA-B*15:02:01_positive/HLA-B*15:02:01_positive,HLA-B*57:01:01_negative/HLA-B*57:01:01_negative,HLA-B*58:01_negative/HLA-B*58:01_negative']",	
        'print_diplotype' : "['*15:02:01/Other,Other/Other,Other/Other']",	
        'tool' : "['ATHLATES,HLA-HD,KOURAMI']"
        }
        ]
    file_path_set = gen_hls_tsv(temp)
    return file_path_set
   		

def  gen_hls_tsv(data:list[dict]):
    file_path_set = []
    for inx,data in enumerate(data):
        df = pd.DataFrame([data])
        file_path = f'{dir_path}/hla_diplotype_tester{inx}.tsv'
        df.to_csv(file_path, sep='\t', encoding='utf-8',index=False)
        file_path_set.append(file_path)
    return file_path_set