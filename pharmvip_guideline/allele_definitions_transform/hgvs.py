import re

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
