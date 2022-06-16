import re 


def  hla_subjection(target:list[dict], guideID, hla_relation:dict):
    '''
    return false if hla type didn't match in hla relation, true otherwise
    '''
    statment = True
    if guideID in hla_relation:
        info = {}
        for i in target:
            info = info | i ['key']
        hla_set = hla_relation[guideID]
        for i in hla_set:
            if i in info:
                check_type = hla_set[i]
                htype = info[i]
                hla_finder = re.match('\*((\d{2}:*)+)',htype)
                if hla_finder:
                    htype = hla_finder.group(0)
                    if htype != check_type:
                        statment = False
                        break

    return statment
