import pandas as pd
import ast

'''
['*1x2', '*1≥3', '*2x2', '*2≥3', '*3x2', 
'*4≥2', '*6x2', '*9x2','*10x2', '*17x2', 
'*29x2' , '*35x2', '*36x2', '*41x2', '*41x3', 
'*43x2', '*45x2']

'''
def  transform_cyp2d6(df:pd.DataFrame):
    guide_dip = df.iloc[0]['guide_diplotype']
    guide_dip = ast.literal_eval(guide_dip)
    diplotype = guide_dip[0]
    diplotype = diplotype.split('/')
    allele = diplotype[0]
    allele, num_copy = extract_copy(allele=allele)
    group = [{'*1', '*2'}, {'*4'}]
    num_copy_int = int(num_copy) if num_copy.isdigit() else 0
    if allele in group[0]:
        if num_copy_int > 2 :
            allele = f'{allele}≥{num_copy}'
        elif num_copy_int == 2:
            f'{allele}x{num_copy}'
    elif allele in group[1]:
        if num_copy_int > 1 :
            allele = f'{allele}≥{num_copy}'
    else:
        allele =  f'{allele}x{num_copy}' if num_copy else allele
        pass
    
    return allele

'''
*1_*4_*13
*1/*10
*1/*10x2
*1/*1x2
*1/*36+*10
*1/*36+*36+*10
*10_*10_*36_*36_*4.013;*10_*36_*36_*36_*4
*10_*10_*36_*83
*10_*36_*36_*36_*36
*10/*106;*1/*52
'''
def  extract_copy(allele:str):
    extraction = allele.split('x')
    allele = extraction[0]
    num_copy = extraction[1] if len(extraction) > 1 else ''
    return allele, num_copy
