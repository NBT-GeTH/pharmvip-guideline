import random
import re
from tests.matcher_testing.gen_matcher_testset import MatcherGenerator

class QueryGenerator(MatcherGenerator):

    def  __init__(self,allele_definition_set=None, num_each_gene=2,missing_rate=0.25,false_rate_in_combine = 0.5):
        super().__init__(allele_definition_set=allele_definition_set, num_each_gene=num_each_gene, missing_rate=missing_rate, false_rate_in_combine=false_rate_in_combine)
        self.basic_base = ['A','T','C','G']

    def  random_gt_base(self) :
        
        new_base = random.choice(self.basic_base)
        return new_base

    def  which_del(self,hgvs):
        find_which_del = re.match((r'^.+del([A-Z]+)$'),hgvs)
        return find_which_del.group(1)

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
        allele = [allele1,allele2]
        
        if hgvs_type == 'INS' :
            new_base = self.random_gt_base()
            for inx,val in enumerate(allele):
                
                if val == 'del':
                    allele[inx] = new_base
                else :
                    allele[inx] = str(allele[inx]).replace('ins',new_base)
        
        elif hgvs_type == 'DEL' : 
            new_base = self.random_gt_base()
            for inx,val in enumerate(allele):
                if val.__contains__('del'):
                    allele[inx] = new_base
                elif val == '.' :
                    find_del = self.which_del(hgvs)
                    allele[inx] = '.' + '.'*len(find_del)
                else:
                    allele[inx] = '.' + allele[inx] if allele[inx] != '.' else allele[inx]
            
        if gt_phase :
            phase = random.choice(['/','|']) if allele1 == allele2 else '|'
        else :
            phase = '/'

        gene["gt_bases"] = allele[0] + phase + allele[1]
        gene["allele1"] = allele[0]
        gene["allele2"] = allele[1]
        
        

        return gene

    def  sample_query_generator(self, gene_name, gene_phase, sample_id=0):
        sample =  super().sample_query_generator(gene_name, gene_phase, sample_id=sample_id)
        sample["chromosome"] = self.allele_definition_set[gene_name]["haplotypes"][0]["chromosome"]
        sample["ana_user_id"] = 1
        sample["ana_id"] = 1
        try :
            del sample["haplotype_allele_name"]
        except:
            pass
        return sample

#%%
# %%
