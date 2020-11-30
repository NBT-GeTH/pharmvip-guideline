import re

def all_rsid(findall_rsid):
    """
    G6PD_allele_definition_table:
    findall_rsid                    all_rsid
    ['rs137852347', 'rs76723693']   'rs137852347, rs76723693'
    """
    all_rsid = ""
    for i in range(len(findall_rsid)):
        if i != (len(findall_rsid) - 1):
            all_rsid += findall_rsid[i] + ", "
        else:
            all_rsid += findall_rsid[i]
    return all_rsid

def findall_rsid(rsid_cell):
    rsid_cell = rsid_cell.values.tolist()
    rsid = []
    for cell in rsid_cell:
        findall_rsid = re.findall(r"rs\d+", str(cell))
        if findall_rsid:
            if len(findall_rsid) > 1:
                """
                G6PD_allele_definition_table:
                cell                            findall_rsid
                ';;rs137852347;rs76723693;'     ['rs137852347', 'rs76723693']
                """
                rsid.append(all_rsid(findall_rsid))
            else:
                """
                cell            findall_rsid
                'rs1800559'     ['rs1800559']
                """
                rsid.append(findall_rsid[0])
        else:
            rsid.append("")
    return rsid
