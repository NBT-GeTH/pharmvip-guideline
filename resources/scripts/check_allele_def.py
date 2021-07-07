#%%
import glob
import json
import os
from sys import path 
from pharmvip_guideline.utils.functional import import_allele_definition_set
import os
pack_path = os.path.join(os.path.dirname(__file__))

class  AnalzeAlleleDef():
    def  __init__(self,ph_cat_def:dict,ph_vip_def:dict) -> None:
        self.set_none_rsid(ph_cat_def)
        self.ph_cat_def = ph_cat_def
        self.ph_vip_def = ph_vip_def
        

    def clear_log(self):
        for log in glob.glob(pack_path + '/*.log'):
            os.remove(log)


    def  write_log(self,title:str,logger:str,file_name:str='name_log',):
        t_bar = '+'*80 + '\n'
        leaft = 15*'%'
        name_log = leaft + title.center(50,'-') + leaft +'\n'
        starter_log = '\n' + t_bar + name_log + t_bar +'\n'
        ender_log =  ('*'*80 + '\n')*2
        logger = starter_log + logger + ender_log
        with open(file_name,'a') as f :
            f.write(logger)

    def set_none_rsid(self,ph_cat_def_set):
        for gene in ph_cat_def_set:
            for inx,val in enumerate(ph_cat_def_set[gene]['variants']) :
                if val['rsid'] == None:
                    ph_cat_def_set[gene]['variants'][inx]['rsid'] = ''
                    # print(ph_cat_def_set[gene]['variants'][inx]['rsid'])
            # break
            
    def find_conflict(self,a_set:set,b_set:set):
        a_set_conflict = []
        in_both = []
        for element in a_set:
            if element in b_set: 
                b_set.remove(element)
                in_both.append(element)
            else:
                a_set_conflict.append(element)

        b_set_conflict = list(b_set)
        # in_both.sort()
        # a_set_conflict.sort()
        return in_both, a_set_conflict, b_set_conflict


    def  compare_allele_name(self,gene:str):
        cat_allele = set()
        vip_allele = set()
        for allele in self.ph_cat_def[gene]["namedAlleles"]:
            allele_name = allele['name']
            cat_allele.add(allele_name)
        for allele in self.ph_vip_def[gene]["name_relation_to_hgvs"] :
            allele_name = allele['name']
            vip_allele.add(allele_name)

        if vip_allele == cat_allele:
            logger = "there is no conflict in allele name\n"
            logger = logger + "allele_name_list : " + str(vip_allele) + '\n\n'
        else:
            in_both, cat_conflict, vip_conflict = self.find_conflict(cat_allele,vip_allele)
            logger = self.get_conflict_log(in_both=in_both,
                cat_conflict=cat_conflict,
                vip_conflict=vip_conflict,
                which_conflict="GENE")

        self.write_log(title=gene, logger=logger, file_name='compare_allele.log')
        
    def  compare_position(self,gene:str,is_write:bool=True):
        
        cat_pos = {}
        cat_pos['chromosome'] = self.ph_cat_def[gene]["variants"][0]["chromosome"]
        cat_pos['variants'] = []
        for variant in self.ph_cat_def[gene]["variants"]:
            data_pack = dict()
            needed_info = {'position','rsid','type'}
            for key,val in variant.items():
                if key in needed_info:
                    data_pack[key] = val
            cat_pos['variants'].append(data_pack)

        #%%
        vip_pos = {}
        vip_pos['chromosome'] = self.ph_vip_def[gene]["haplotypes"][0]["chromosome"]
        vip_pos['variants'] = []
        logger = ''
        for variant in self.ph_vip_def[gene]["haplotypes"][0]["variants"]:
            data_pack = {}
            data_pack['position'] = int(variant['start'])
            data_pack['rsid'] = variant['rsid']
            data_pack['type'] = variant['hgvs_type']
            vip_pos['variants'].append(data_pack)

        if vip_pos == cat_pos:
            logger = "there is no conflict in posision\n"
            logger = logger + "position_list : " + str(vip_pos) + '\n\n'
        else:
            if vip_pos['chromosome'] != cat_pos['chromosome']:
                logger = logger + 'Conflict in  chromosome possition have found\n'
            
            if vip_pos['variants'] != cat_pos['variants']:
                cat_variants = cat_pos['variants']
                vip_variants = vip_pos['variants']
                in_both, cat_conflict, vip_conflict  = self.find_conflict(cat_variants,vip_variants)
                logger = logger + self.get_conflict_log(in_both=in_both,
                    cat_conflict=cat_conflict,
                    vip_conflict=vip_conflict,
                    which_conflict="variant")

        if is_write : self.write_log(title=gene, logger=logger, file_name='compare_possition.log')
        return cat_pos, vip_pos

    def get_conflict_log(self,in_both:set, 
            cat_conflict:set, 
            vip_conflict:set,
            which_conflict:str=''):
        in_both, cat_conflict, vip_conflict = ('{∅}' if str(element).__contains__('[]') else element for element in (in_both, cat_conflict, vip_conflict) )
        logger = f'Conflict in {which_conflict} possition have found\n'
        star = '*'*10
        logger = logger + 'there is ' + str(in_both) + f' {which_conflict} in both system.\n'
        logger = logger + f'{star}CAT*{star}\n'
        logger = logger + 'there is ' + str(cat_conflict) + f' {which_conflict} in PharmCAT system but do not have in PharmVIP system.\n'
        logger = logger + f'{star}VIP{star}\n'
        logger = logger + 'there is ' + str(vip_conflict) + f' {which_conflict} in PharmVIP system but do not have in PharmCAT system.\n'
        return logger


    def  check_gene(self,ilog:bool=True):
        cat_gene = {gene for gene in self.ph_cat_def}
        vip_gene = {gene for gene in self.ph_vip_def}
        in_both = cat_conflict = vip_conflict = set()
        if cat_gene == vip_gene :
            logger = "there is no conflict in GENE list\n"
            logger = logger + "GENE list : " + str(cat_gene) + '\n\n'
            in_both = cat_gene
        else:
            in_both, cat_conflict, vip_conflict = self.find_conflict(cat_gene,vip_gene)
            logger = self.get_conflict_log(in_both=in_both,
                cat_conflict=cat_conflict,
                vip_conflict=vip_conflict,
                which_conflict="GENE")
            
        if ilog: self.write_log(title="Compare Gene",logger=logger,file_name='compare_gene.log')
        return in_both, cat_conflict, vip_conflict
    
#%%
ph_cat_path = '/home/xixe/phCat/0.8/PharmCAT-0.8.0/build/resources/main/org/pharmgkb/pharmcat/definition/alleles'
ph_vip_def_set = import_allele_definition_set()
ph_cat_def_set = import_allele_definition_set(ph_cat_path)
tooler = AnalzeAlleleDef(ph_cat_def=ph_cat_def_set,ph_vip_def=ph_vip_def_set)
tooler.clear_log()
in_both, cat_conflict, vip_conflict = tooler.check_gene()
for gene in in_both:
    tooler.compare_allele_name(gene=gene)
    _, _ = tooler.compare_position(gene=gene)
    
#%%
ph_cat_path = '/home/xixe/phCat/0.8/PharmCAT-0.8.0/build/resources/main/org/pharmgkb/pharmcat/definition/alleles'
ph_vip_def_set = import_allele_definition_set()
ph_cat_def_set = import_allele_definition_set(ph_cat_path)
#%%
gene = 'CACNA1S'
ph_cat_def_set[gene]
# %%
ph_cat_def_set[gene]['namedAlleles']

# %%
