# %%
from ast import ExtSlice
import random
import glob
import json
import pickle
import re 
import os
from typing import Type
from pharmvip_guideline import *
from tests.query_testing.gen_query_testset import QueryGenerator
from pharmvip_guideline.utils.functional import import_allele_definition_set
path_tomodule = os.path.join(os.path.dirname(__file__))

#%%
class  VCFsFileGenerator(QueryGenerator):
    def  __init__(self,
            allele_definition_set=None, 
            num_each_gene=2,
            missing_rate=0.25,
            false_rate_in_combine = 0.5):
        super().__init__(allele_definition_set=allele_definition_set, 
            num_each_gene=num_each_gene, 
            missing_rate=missing_rate, 
            false_rate_in_combine=false_rate_in_combine)
        self.possition_collector = self.load_possition_collector()
        self.body = ''
        self.header = ''

    def  sample_collector_to_VCFs(self):
        for inx,sample in enumerate(self.get_sample_collector()):
            self.from_query_set_to_VCFs(sample)

            # self.get_sample_collector()[inx]['header'] = self.header
            # self.get_sample_collector()[inx]['body'] = self.body
            self.clear_body()

    def  load_possition_collector(self):
        with open(defaults_allele_definitions_dbpmcgenomics + '/gPOS_collector.pickle', 'rb') as f:
            data = pickle.load(f)
        return data

    def  from_hgvs_to_possition(self,hgvs):
        pos_finder = re.match((r'^g\.(\d+)(\_.*|\w.*)$'),hgvs)
        return pos_finder.group(1)

    def  clear_body(self):
        self.body = ''

    def  write_vcf(self,sample_id):
        with open(f'./vcfs/{sample_id}.vcf', "w") as f :
            for key,val in self.header.items():
                f.write(val)
            f.write(self.body)

    def  from_query_set_to_VCFs(self,query):
        sample_id = query['sample_id']
        self.header = self.get_vcfs_header(sample_id)
        for variant in query['variants']:
            self.body = self.body + self.variant_to_VCF_body(variant,chrom=query['chromosome'])
        self.write_vcf(sample_id)

    def  run_vcfs_shell_script(self):
        os.system(f'bash {path_tomodule}/script/vcf_handle_withARG.sh {path_tomodule}/vcfs' )

    def  get_vcfs_header(self,sample_id):
        header = {}
        header['fileformat'] = "##fileformat=VCFv4.2\n"
        header['source'] = '##source=SelectVariants\n'
        header['filter'] = '##FILTER=<ID=PASS,Description="All filters passed">\n'
        header['info'] = '##INFO=<ID=PX,Number=.,Type=String,Description="PGX">\n'
        # header['info'] = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n'
        header['format'] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n'
        header['format2'] = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n'
        header['table_head'] = f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n'
        return header

    def  find_genome_inx(self,genome_list,allele) :
        allele = allele.replace('.','')
        try:
            genome_inx = genome_list.index(allele)
        except:
            genome_inx = '.'
        return genome_inx

    def  find_phase_pipe(self,gt_bases):
        phase_pipe = '/' if '/' in gt_bases else '|'
        return phase_pipe

    def  variant_to_VCF_body(self,variant,chrom):
        hgvs = variant['hgvs']
        hgvs_start_pos = self.from_hgvs_to_possition(hgvs)
        dp = variant['dp']
        hgvs_type = variant['hgvs_type']
        ref = None
        alt = self.possition_collector[hgvs_start_pos]['alt']
        phase = self.find_phase_pipe(variant['gt_bases'])
        if hgvs_type == "SNP" or hgvs_type == "CNV": #test by using CACNA1S
            ref = self.possition_collector[hgvs_start_pos]['ref']
            alt_fomat =  ','.join(alt)
            genome_list = [ref] + list(alt)
            genome1 = self.find_genome_inx(genome_list,variant['allele1'])
            genome2 = self.find_genome_inx(genome_list,variant['allele2'])
            # phase = self.find_phase_pipe(variant['gt_bases'])
            body = f'{chrom}\t{hgvs_start_pos}\t.\t{ref}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'

        if hgvs_type == 'INS': #test by using CYP3A5
            allele_list = []
            for inx,allele in enumerate([variant['allele1'],variant['allele2']]):
                if allele != 'del' and allele != '.' :
                    ref = allele
                allele_list.append(allele)

            if ref == None:
                ref = self.random_gt_base()

            allele_list = [ref if allele == 'del' else allele for allele in allele_list ]
            alt = [genome.replace('ins',ref) for genome in alt]
            alt_fomat =  ','.join(alt)
            genome_list = [ref] + list(alt)
            genome1 = self.find_genome_inx(genome_list,allele_list[0])
            genome2 = self.find_genome_inx(genome_list,allele_list[1])
            # phase = self.find_phase_pipe(variant['gt_bases'])
            body = f'{chrom}\t{hgvs_start_pos}\t.\t{ref}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'

            
        
        if hgvs_type == 'DEL':
            ref = self.possition_collector[hgvs_start_pos]['ref']
            genome_list = [ref] + list(alt)
            genome1 = self.find_genome_inx(genome_list,variant['allele1_convert'])
            genome2 = self.find_genome_inx(genome_list,variant['allele2_convert'])

            find_which_del = re.match((r'^g\.(\d+)\_(\d+)del([A-Z]+)$'),hgvs)
            if find_which_del == None :
                find_which_del = re.match((r'^g\.(\d+)del([A-Z]+)$'),hgvs)
                pos_list = int(find_which_del.group(1))
                allele_list = find_which_del.group(2)
                pos_befor = int(pos_list) - 1
            else :
                pos_list = [i for i in range(int(find_which_del.group(1)),int(find_which_del.group(2)) + 1)]
                allele_list = find_which_del.group(3)
                pos_befor = int(pos_list[0]) - 1

            ref_base = None
            for inx,allele in enumerate([variant['allele1'],variant['allele2']]):
                if allele in self.basic_base:
                    ref_base = allele
            
            if ref_base == None:
                ref_base = self.random_gt_base()
            
            # for_extended_pos = genome_list[inx] 
            for inx,allele in enumerate(genome_list):
                if allele.__contains__('del'):
                    genome_list[inx] = ref_base
                else:
                    genome_list[inx] = ref_base + genome_list[inx] 

            ref_format = genome_list.pop(0)
            alt_fomat =  ','.join(genome_list)
            body = f'{chrom}\t{pos_befor}\t.\t{ref_format}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'
            alt = [ '*' * len(ref)  if allele.__contains__('del') else allele for allele in alt]
            extend_alt =  ['*'] * len(ref)
            extend_alt_inx1 = extend_alt_inx2 = [0] * len(ref)

            for ainx,allele in enumerate(alt):
                for binx,base in enumerate(allele):
                    if (base not in ref[binx]) and (base not in extend_alt[binx]):
                        extend_alt[binx] = extend_alt[binx] + base
                        if genome1 - 1 == ainx : 
                            extend_alt_inx1[binx] = len(extend_alt[binx]) - 1 
                        if genome2 - 1 == ainx : 
                            extend_alt_inx2[binx] = len(extend_alt[binx]) - 1 

            if type(pos_list) == int :
                inx = 0
                body = body + f'{chrom}\t{pos_list}\t.\t{ref}\t{extend_alt[inx]}\t.\tPASS\tPX=;\tGT:DP\t{extend_alt_inx1[inx]}{phase}{extend_alt_inx1[inx]}:{dp}\n'
            else:
                for inx,poss in enumerate(pos_list):
                    extend_alt[inx] = ','.join(extend_alt[inx])
                    try:
                        body = body + f'{chrom}\t{poss}\t.\t{ref[inx]}\t{extend_alt[inx]}\t.\tPASS\tPX=;\tGT:DP\t{extend_alt_inx1[inx]}{phase}{extend_alt_inx1[inx]}:{dp}\n'
                    except:
                        print(inx,ref,extend_alt,extend_alt_inx1,extend_alt_inx2)
            
        return body


#%%
#%%
