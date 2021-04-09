import re
import copy

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

def manual_customize(allele_definition_df, gene):
    if gene == "G6PD":
        allele_definition_df.iloc[1, 0] = allele_definition_df.iloc[1, 1]
        allele_definition_df.iloc[2, 0] = allele_definition_df.iloc[2, 1]
        allele_definition_df.iloc[3, 0] = allele_definition_df.iloc[3, 1]
        allele_definition_df.iloc[4, 0] = allele_definition_df.iloc[4, 1]
        allele_definition_df.iloc[5, 0] = allele_definition_df.iloc[5, 1]

        allele_definition_df.drop([1], axis=1, inplace=True)
        allele_definition_df.columns = range(allele_definition_df.shape[1])

        return allele_definition_df
    elif gene == "TPMT":
        exceptions = None
        for column in range(allele_definition_df.shape[1]):
            if allele_definition_df.iloc[3, column] ==  "g.18138983G>A":
                exceptions = column
        allele_definition_df.drop(allele_definition_df.columns[[exceptions]], axis=1, inplace=True)
        allele_definition_df.columns = range(allele_definition_df.shape[1])

        return allele_definition_df
    elif gene == "VKORC1":
        for row in range(allele_definition_df.shape[0]):
            if allele_definition_df.iloc[row, 0] ==  "rs9923231 reference (C)":
                allele_definition_df.iloc[row, 0] = "-1639G"
            elif allele_definition_df.iloc[row, 0] ==  "rs9923231 variant (T)":
                allele_definition_df.iloc[row, 0] = "-1639A"
        
        return allele_definition_df
    else:
        return allele_definition_df

def automatic_customize(allele_definition_df):
    allele_definition_df = sort_hgvs(allele_definition_df)

    allele_definition_df = allele_definition_df.iloc[:find_end_row(allele_definition_df), :find_end_col(allele_definition_df)]

    return allele_definition_df

def sort_hgvs(allele_definition_df):
    allele_definition_df.sort_values(axis=1, by=3, inplace=True)
    allele_definition_df.columns = range(allele_definition_df.shape[1])
    return allele_definition_df

def find_end_row(allele_definition_df):
    for index, cell in allele_definition_df.iloc[:, 0].items():
        if str(cell) == "nan" or re.match(r"^NOTES:\s*$", str(cell)):
            return index

def find_end_col(allele_definition_df):
    for index, cell in allele_definition_df.iloc[3, :].items():
        if str(cell) == "nan":
            return index

def clean_allele_cell(allele_definition_df):
    allele_definition_df = clean_nan_allele_cell(allele_definition_df)

    allele_definition_df = clean_whitespace_allele_cell(allele_definition_df)

    return allele_definition_df

def clean_nan_allele_cell(allele_definition_df):
    allele_definition_df_ = copy.deepcopy(allele_definition_df)
    for i in range(allele_definition_df_.shape[0]):
        if i >= 7:
            for j in range(allele_definition_df_.shape[1]):
                if j >= 1:
                    if str(allele_definition_df_.iloc[i, j]) == "nan":
                        allele_definition_df_.iloc[i, j] = allele_definition_df_.iloc[7, j]
    return allele_definition_df_

def clean_whitespace_allele_cell(allele_definition_df):
    allele_definition_df_ = copy.deepcopy(allele_definition_df)
    for i in range(allele_definition_df_.shape[0]):
        if i >= 7:
            for j in range(allele_definition_df_.shape[1]):
                if j >= 1:
                    if " " in str(allele_definition_df_.iloc[i, j]):
                        allele_definition_df_.iloc[i, j] = allele_definition_df_.iloc[i, j].replace(" ", "")
    return allele_definition_df_

def search_chromosome(position_cell):
    search_chromosome = re.search(r"chromosome\s(\w+)", str(position_cell))
    if search_chromosome:
        """
        position_cell                                                       search_chromosome
        'Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)'   1. '1'
        'Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)'   1. '7'
        """
        return position_cell, "chr" + search_chromosome.group(1)
    else:
        print(f"error search chromosome with: {str(position_cell)}")
        exit()

def match_hgvs(hgvs_cell):
    hgvs_cell = hgvs_cell.values.tolist()
    hgvs_type = []
    start = []
    end = []
    for cell in hgvs_cell:
        match_snp = re.match(r"^g\.((\d+)([A-Z]>([A-Z]|[A-Z](\/[A-Z])*)))*$", str(cell))
        match_ins = re.match(r"^g\.(\d+)_?(\d*)ins.+$", str(cell))
        match_del = re.match(r"^g\.(\d+)_?(\d*)del.+$", str(cell))
        match_cnv = re.match(r"^g\.(\d+)", str(cell))
        if match_snp:
            """
            cell                match_snp
            'g.201060815C>T'    1. '201060815C>T'
                                2. '201060815'
                                3. 'C>T'
                                4. 'T'
                                5.
            'g.40991390C>A/T'   1. 40991390C>A/T
                                2. 40991390
                                3. C>A/T
                                4. A/T
                                5. /T
            'g.41004406G>A/C/T' 1. 41004406G>A/C/T
                                2. 41004406
                                3. G>A/C/T
                                4. A/C/T
                                5. /T
            """
            hgvs_type.append("SNP")
            start.append(match_snp.group(2))
            end.append(match_snp.group(2))
        elif match_ins:
            """
            cell                                                match_ins
            'g.42126666_42126667insAGTGGGCAC'                   1. '42126666'
                                                                2. '42126667'
            'g.42128936_42128937insGGGGCGAAA/insGGGGCGAAAGGGGC' 1. '42128936'
                                                                2. '42128937'
            """
            hgvs_type.append("INS")
            start.append(match_ins.group(1))
            end.append(match_ins.group(2))
        elif match_del:
            """
            cell                                match_del
            'g.94942213_94942222delAGAAATGGAA'  1. '94942213'
                                                2. '94942222'
            'g.94949283delA'                    1. '94949283'
                                                2. ''
            """
            hgvs_type.append("DEL")
            start.append(match_del.group(1))
            if match_del.group(1) and not match_del.group(2):
                end.append(match_del.group(1))
            elif match_del.group(1) and match_del.group(2):
                end.append(match_del.group(2))
        elif match_cnv:
            """
            cell            match_cnv
            'g.233760233'   1. '233760233'
            """
            hgvs_type.append("CNV")
            start.append(match_cnv.group(1))
            end.append(match_cnv.group(1))
        else:
            print(f"error match hgvs with: {str(cell)}")
            exit()
    return hgvs_cell, hgvs_type, start, end

def findall_rsid(rsid_cell):
    rsid_cell = rsid_cell.values.tolist()
    rsid = []
    for cell in rsid_cell:
        findall_rsid = re.findall(r"rs\d+", str(cell))
        if findall_rsid:
            """
            cell                            findall_rsid
            'rs1800559'                     ['rs1800559']
            ';;rs137852347;rs76723693;'     ['rs137852347', 'rs76723693']
            """
            rsid.append(", ".join(findall_rsid))
        else:
            rsid.append("")
    return rsid

def extract_allele(haplotype_cell):
    haplotype_cell = haplotype_cell.values.tolist()
    haplotype = []
    allele = []
    for cell in range(len(haplotype_cell)):
        haplotype.append(haplotype_cell[cell][0])
        allele.append(haplotype_cell[cell][1:])
    return haplotype, allele
