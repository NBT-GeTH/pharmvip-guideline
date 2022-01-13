import pandas as pd
import ast
import copy

'''
['*1x2', '*1≥3', '*2x2', '*2≥3', '*3x2', 
'*4≥2', '*6x2', '*9x2','*10x2', '*17x2', 
'*29x2' , '*35x2', '*36x2', '*41x2', '*41x3', 
'*43x2', '*45x2']
'''

def  transform_cyp2d6(diplotypes:list[str]):
    group = [{'*1', '*2'}, {'*4'}]
    for ix, diplotype in enumerate(diplotypes):    
        diplotype = diplotype.split('/')
        for inx, allele in enumerate(diplotype):
            allele, num_copy = extract_copy(allele=allele)
            
            num_copy_int = int(num_copy) if num_copy.isdigit() else 0
            if allele in group[0]:
                if num_copy_int > 2 :
                    allele = f'{allele}≥{num_copy}'
                elif num_copy_int == 2:
                    allele = f'{allele}x{num_copy}'
            elif allele in group[1]:
                if num_copy_int > 1 :
                    allele = f'{allele}≥2'
            else:
                allele =  f'{allele}x{num_copy}' if num_copy else allele
                pass
            diplotype[inx] = allele
        diplotype = '/'.join(diplotype)
        diplotypes[ix] = diplotype
    return diplotypes

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

def  detect_or(diplotypes:list[str]):
    diplotypes = [i.split(';') for i in diplotypes]
    diplotypes = [item for sublist in diplotypes for item in sublist]
    return diplotypes

def  cyp2d6_hanlder(diplotypes:list[str]):
    diplotypes = detect_or(diplotypes)
    diplotypes = transform_cyp2d6(diplotypes)
    return diplotypes
