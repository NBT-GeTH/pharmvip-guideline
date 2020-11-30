def extract_allele(haplotype_cell):
    haplotype_cell = haplotype_cell.values.tolist()
    haplotype = []
    allele = []
    for cell in range(len(haplotype_cell)):
        haplotype.append(haplotype_cell[cell][0])
        allele.append(haplotype_cell[cell][1:])
    return haplotype, allele
