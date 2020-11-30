import re

def match_gene(gene_cell):
    match_gene = re.match(r"^GENE:\s*(\w+)\s*$", str(gene_cell))
    if match_gene:
        """
        gene_cell           match_gene
        'GENE: CYP2C19 '    1. 'CYP2C19'
        'GENE: CYP2D6'      1. 'CYP2D6'
        'GENE:CYP3A5'       1. 'CYP3A5'
        """
        return match_gene.group(1)
    else:
        print(f"error match gene with: {str(gene_cell)}")
        exit()
