import sys, os

from pandas.core.frame import DataFrame
import pandas  as pd

def disable_print():
    sys.stdout = open(os.devnull, "w")

def enable_print():
    sys.stdout = sys.__stdout__

def  write_exel(df:DataFrame,path:str='./'):
    writer = pd.ExcelWriter(f'{path}comparing.xlsx', engine='xlsxwriter')
    df.to_excel(writer,index=None)
    writer.save()