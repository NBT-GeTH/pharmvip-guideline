import copy
import re

def locate_max(a):
    biggest = max(a)
    return biggest, [index for index, element in enumerate(a) if biggest == element]

def find_best_candidate(allele_definition, allele_matcher):
    raw_count_diplotype = copy.deepcopy(allele_matcher["count_diplotype"])
    raw_guide_dip = copy.deepcopy(allele_matcher["guide_dip"])
    raw_print_dip = copy.deepcopy(allele_matcher["print_dip"])
    
    missing_position = []
    guide_dip_remove_list = []
    print_dip_remove_list = []
    for variant in allele_matcher["variants"]:
        if re.match(r"^(\.+)(\/|\|)(\.+)$", variant["gt_bases"]):
            missing_position.append(variant["hgvs"])
            for relation in allele_definition["hgvs_relation_to_name"]:
                if relation["hgvs"] == variant["hgvs"]:
                    for name in relation["name"]:
                        for guide_dip in allele_matcher["guide_dip"]:
                            if name in guide_dip:
                                guide_dip_remove_list.append(guide_dip)
                        for print_dip in allele_matcher["print_dip"]:
                            if name in print_dip:
                                print_dip_remove_list.append(print_dip)
    for guide_dip_remove in guide_dip_remove_list:
        if guide_dip_remove in allele_matcher["guide_dip"]:
            allele_matcher["guide_dip"].remove(guide_dip_remove)
    for print_dip_remove in print_dip_remove_list:
        if print_dip_remove in allele_matcher["print_dip"]:
            allele_matcher["print_dip"].remove(print_dip_remove)

    assert len(allele_matcher["guide_dip"]) == len(allele_matcher["print_dip"])
    if not allele_matcher["guide_dip"] and not allele_matcher["print_dip"]:
        allele_matcher["count_diplotype"] = raw_count_diplotype
        allele_matcher["guide_dip"] = raw_guide_dip
        allele_matcher["print_dip"] = raw_print_dip
    else:
        if len(allele_matcher["guide_dip"]) == 1 and len(allele_matcher["print_dip"]) == 1:
            allele_matcher["count_diplotype"] = len(allele_matcher["guide_dip"])
            return allele_matcher
    
    guide_dip_score = []
    print_dip_score = []
    for guide_dip in allele_matcher["guide_dip"]:
        name1 = guide_dip.split("/")[0]
        name2 = guide_dip.split("/")[1]
        score1 = 0
        score2 = 0
        for relation in allele_definition["name_relation_to_hgvs"]:
            if name1 == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                score1 = len(relation["hgvs"])
                break
        for relation in allele_definition["name_relation_to_hgvs"]:
            if name2 == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                score2 = len(relation["hgvs"])
                break
        guide_dip_score.append(score1 + score2)
    for print_dip in allele_matcher["print_dip"]:
        name1 = print_dip.split("/")[0]
        name2 = print_dip.split("/")[1]
        score1 = 0
        score2 = 0
        for relation in allele_definition["name_relation_to_hgvs"]:
            if name1 == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                score1 = len(relation["hgvs"])
                break
        for relation in allele_definition["name_relation_to_hgvs"]:
            if name2 == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                score2 = len(relation["hgvs"])
                break
        print_dip_score.append(score1 + score2)

    guide_dip_new = []
    print_dip_new = []
    for index in locate_max(guide_dip_score)[1]:
        guide_dip_new.append(allele_matcher["guide_dip"][index])
    for index in locate_max(print_dip_score)[1]:
        print_dip_new.append(allele_matcher["print_dip"][index])

    allele_matcher["count_diplotype"] = len(print_dip_new)
    allele_matcher["guide_dip"] = guide_dip_new
    allele_matcher["print_dip"] = print_dip_new

    return allele_matcher
