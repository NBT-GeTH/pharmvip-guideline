
def  hla_ingroup(guideID,info):
    guide_map = {
        '100412': {'gene' : 'HLA-B','type' : '15:02'}, 
        '100421': {'gene' : 'HLA-B','type' :'57:01'},
        '100422' : {'gene' : 'HLA-B','type' :'58:01'},}
    target_gene = guide_map[guideID]['gene']
    target_type = guide_map[guideID]['type']
    state = False if  not (target_type in info[target_gene]) else True
    return state


def  hla_specific(info):
    num_gene = len(info)
    a_type = '31:01'
    b_type = '15:02'
    if (num_gene == 1) and not(b_type in info['HLA-B']) : return False
    elif num_gene == 2 :
        a_bool = not(a_type in info['HLA-A'])
        b_bool = not(b_type in info['HLA-B'])
        return not(a_bool or b_bool)


def  hla_subjection(target:list[dict], guideID):
    statment = True
    info = {}
    for i in target:
        info = info | i ['key']
    group = ['100412','100421','100422']
    if guideID in group : statment = hla_ingroup(guideID,info)
    elif guideID == '100423': statment = hla_specific(info)

    return statment
