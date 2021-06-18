import glob
import json
import random
from pharmvip_guideline import *
from pharmvip_guideline.allele_matcher.query import convert_allele

class Generator:

    def __init__(self,allele_definition_set=None, num_each_gene=2,missing_rate=0.25,false_rate_in_combine = 0.5,dot_flag=True):
        self.num_each_gene = num_each_gene
        self.allele_definition_set = allele_definition_set
        self.missing_rate = missing_rate
        self.false_rate_in_combine = false_rate_in_combine
        self.__sample_collector = []


    def random_missing_genome(self,genome):
        return random.choices([genome,'.'],[1 - self.missing_rate,self.missing_rate])[0]
    
    def random_missing_pair_genome(self,genome1,genome2):
        return genome1,genome2 if random.choices([True,False],[1 - self.missing_rate,self.missing_rate])[0] else '.','.'

    def random_gt_phase(self):
        return random.choices([True,False],[1 - self.false_rate_in_combine,self.false_rate_in_combine])[0]
   
    def random_haplotype_pair(self,gene_name):
       return random.choices(self.allele_definition_set[gene_name]["haplotypes"],k=2)
   
    def sample_query_generator(self,gene_name,gene_phase, sample_id=0):
        sample_query = {
            "sample_id": sample_id,
            "gene": self.allele_definition_set[gene_name]["gene"],
            "variants": [],
            "gene_phases": gene_phase
        }
        haplotype1, haplotype2 = self.random_haplotype_pair(gene_name=gene_name)
        haplotype1_name = haplotype1['name']
        haplotype2_name = haplotype2['name']
        sample_query["haplotype_allele_name"] = f"{haplotype1_name}/{haplotype2_name}"
        number_of_variant = len(haplotype1['variants'])
        missing_phase_count = 0
        for var_inx in range(number_of_variant):
            allele1 = self.random_missing_genome(genome=haplotype1['variants'][var_inx]["allele"])
            allele2 = self.random_missing_genome(genome=haplotype2['variants'][var_inx]["allele"])
            
            
            if allele1 == '.' or allele2 == '.' :
                missing_phase_count = missing_phase_count + 1 
                gt_phase = '.'
            elif allele1 == allele2:
                gt_phase = True
            elif gene_phase == "Combine":
                gt_phase = self.random_gt_phase()
            else:
                gt_phase = gene_phase
           
            gene = {
            "hgvs_type": haplotype1['variants'][var_inx]["hgvs_type"],
            "allele1_convert": allele1,
            "allele2_convert": allele2,
            "gt_phases": gt_phase
            }
            sample_query["variants"].append(gene)

        if missing_phase_count == number_of_variant:
            sample_query["gene_phases"] = '.'

        phase_checker = [var['gt_phases'] for var in  sample_query['variants']]
        if gene_phase == "Combine":
            if (True not in phase_checker) or (False not in phase_checker):
                return self.sample_query_generator(gene_name=gene_name,gene_phase=gene_phase,sample_id=sample_id)
        elif gene_phase == False:  
            if True in phase_checker : return self.sample_query_generator(gene_name=gene_name,gene_phase=gene_phase,sample_id=sample_id)

        return sample_query

    def massive_generation(self,gene_name,gene_phase,id_prefix='A'):
        for i in range(self.num_each_gene):
            sample_query = self.sample_query_generator(gene_name=gene_name,gene_phase=gene_phase,sample_id=f'sample_{id_prefix}{i}_{gene_name}')
            self.__sample_collector.append(sample_query)

    def overall_generation(self,gene_phase,id_prefix='A'):
        for allele_definition in self.allele_definition_set:
            self.massive_generation(gene_name=allele_definition,gene_phase=gene_phase,id_prefix=id_prefix)
        

    def get_sample_collector(self):
        return self.__sample_collector

    def clear_sample_query_collector(self):
        self.__sample_collector = []

'''''''''''''''''''''''''''''''''
=================================
========= Just for test =========
=================================    
'''''''''''''''''''''''''''''''''
# # %%
# allele_definition_set = {}
# for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
#         allele_definition = json.load(open(allele_definition_json))
#         allele_definition_set[allele_definition["gene"]] = allele_definition
# #%%
# gener = Generator(allele_definition_set)

# # %%
# gener.num_each_gene = 20
# gener.massive_generation(gene_name='CACNA1S',gene_phase="Combine",id_prefix="TA")
# # %%
# gener.get_sample_collector()

# # %%
# # %%
# test = gener.sample_query_generator(gene_name='CACNA1S',gene_phase="combine")
# test
# #%%
