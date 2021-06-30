# %%
import random
import glob
import json
import pickle
import re 
import os
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

    def  load_possition_collector(self):
        with open(defaults_allele_definitions_dbpmcgenomics + '/gPOS_collector.pickle', 'rb') as f:
            data = pickle.load(f)
        return data

    def  from_hgvs_to_possition(self,hgvs):
        pos_finder = re.match((f'^g\.(\d+)(\_.*|\w.*)$'),hgvs)
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

    def  run_vcds_shell_script(self):
        os.system(f'bash {path_tomodule}/script/vcf_handle_withARG.sh {path_tomodule}/vcfs' )

    def  get_vcfs_header(self,sample_id):
        header = {}
        header['fileformat'] = "##fileformat=VCFv4.2\n"
        header['source'] = '##source=SelectVariants\n'
        header['filter'] = '##FILTER=<ID=PASS,Description="All filters passed">\n'
        header['info'] = '##INFO=<ID=PX,Number=.,Type=String,Description="PGX">\n'
        # header['info'] = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n'
        header['format'] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n'
        header['format2'] = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
        header['table_head'] = f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n'
        return header

    def  find_genome_inx(self,genome_list,allele) :
        try:
            genome_inx = genome_list.index(allele)
        except:
            genome_inx = '.'
        return genome_inx

    def  find_phase_pipe(self,gt_bases):
        phase_pipe = '/' if '/' in gt_bases else '|'
        return phase_pipe

    def  variant_to_VCF_body(self,variant,chrom):
        pos = self.from_hgvs_to_possition(variant['hgvs'])
        dp = variant['dp']
        hgvs_type = variant['hgvs_type']
        ref = None
        alt = self.possition_collector[pos]['alt']
        if hgvs_type == "SNP" or hgvs_type == "CNV": #test by using CACNA1S
            ref = self.possition_collector[pos]['ref']
            alt_fomat =  ','.join(alt)
            genome_list = [ref] + list(alt)
            genome1 = self.find_genome_inx(genome_list,variant['allele1'])
            genome2 = self.find_genome_inx(genome_list,variant['allele2'])
            phase = self.find_phase_pipe(variant['gt_bases'])
            body = f'{chrom}\t{pos}\t.\t{ref}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'

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
            phase = self.find_phase_pipe(variant['gt_bases'])
            body = f'{chrom}\t{pos}\t.\t{ref}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'
        
        return body

    # def  handle_body(self,chrom,pos,ref,alt,variant):
    #     alt_fomat =  ','.join(alt)
    #     genome_list = [ref] + list(alt)
    #     genome1 = self.find_genome_inx(genome_list,variant['allele1'])
    #     genome2 = self.find_genome_inx(genome_list,variant['allele2'])
    #     phase = self.find_phase_pipe(variant['gt_bases'])
    #     body = f'{chrom}\t{pos}\t.\t{ref}\t{alt_fomat}\t.\tPASS\tPX=;\tGT:DP\t{genome1}{phase}{genome2}:{dp}\n'
