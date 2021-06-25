import random
from tests.matcher_testing.gen_matcher_testset import MatcherGenerator

class QueryGenerator(MatcherGenerator):

    def  __init__(self,allele_definition_set=None, num_each_gene=2,missing_rate=0.25,false_rate_in_combine = 0.5):
        super().__init__(allele_definition_set=allele_definition_set, num_each_gene=num_each_gene, missing_rate=missing_rate, false_rate_in_combine=false_rate_in_combine)

    # def  handle_gt_bases(self,gene, varian_pack, allele1, allele2, gt_phase):
    #     if gt_phase :
    #         phase = random.choice(['/','|']) if allele1 == allele2 else '|'
    #         gene["gt_bases"] = f'{allele1}{phase}{allele2}'
    #         gene["allele1"] = allele1
    #         gene["allele2"] = allele2

    def  gen_gene_variant(self, varian_pack, allele1, allele2, gt_phase):
        haplotype = varian_pack[0]
        var_inx = varian_pack[1]
        dp_min = 0
        dp_max = 1000
        hgvs = haplotype["variants"][var_inx]["hgvs"]
        hgvs_type = haplotype["variants"][var_inx]["hgvs_type"]
        gene = {
            "hgvs": hgvs,
            "hgvs_type": hgvs_type,
            "start": haplotype["variants"][var_inx]["start"],
            "end": haplotype["variants"][var_inx]["end"],
            "rsid": haplotype["variants"][var_inx]["rsid"],
            "dp": int(random.uniform(dp_min,dp_max)),
            "allele1_convert": allele1,
            "allele2_convert": allele2,
            "gt_phases": gt_phase
        }
        allele = allele1,allele2
        phase = '/'
        # phase = random.choice(['/','|']) if allele1 == allele2 else '|'
        # gt_bases = f'{mini_a1}{phase}{mini_a2}'
        # basic_base = ['A','T','C','G']
        
        if hgvs_type == 'INS' :
            for i in (allele):
                if i == 'del' : 
            # mini_a1 = random.choice(basic_base)
            # mini_a2 = random.choice(basic_base)
            # gt_bases = f'{mini_a1}{allele1}{phase}{mini_a2}{allele2}'
            
        if gt_phase :
            gene["gt_bases"] = 
            gene["allele1"] = allele[1]
            gene["allele2"] = allele[2]
        

        return gene

    def  sample_query_generator(self, gene_name, gene_phase, sample_id=0):
        sample =  super().sample_query_generator(gene_name, gene_phase, sample_id=sample_id)
        sample["chromosome"] = self.allele_definition_set[gene_name]["haplotypes"][0]["chromosome"]
        del sample["haplotype_allele_name"]
        return sample

