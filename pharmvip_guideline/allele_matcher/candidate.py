import copy
import re

<<<<<<< HEAD
def locate_the_most(score):
    the_most = max(score)
    return the_most, [index for index, element in enumerate(score) if the_most == element]
=======
def locate_max(a):
    biggest = max(a)
    return biggest, [index for index, element in enumerate(a) if biggest == element]
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c

def find_best_candidate(allele_definition, allele_matcher):
    raw_count_diplotype = copy.deepcopy(allele_matcher["count_diplotype"])
    raw_guide_dip = copy.deepcopy(allele_matcher["guide_dip"])
    raw_print_dip = copy.deepcopy(allele_matcher["print_dip"])
    
    missing_position = []
    name_relation_to_missing_hgvs = []
<<<<<<< HEAD
=======
    #loop to find allele with missing position
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    for variant in allele_matcher["variants"]:
        if re.match(r"^(\.+)(\/|\|)(\.+)$", variant["gt_bases"]):
            missing_position.append(variant["hgvs"])
            for relation in allele_definition["hgvs_relation_to_name"]:
                if relation["hgvs"] == variant["hgvs"]:
                    for name in relation["name"]:
                        name_relation_to_missing_hgvs.append(name)

<<<<<<< HEAD
=======
    #loop to remove allele name which is missing in prindip and guidip
    got_remove = []#for test 
    looper = []#for test 
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    for relation in allele_definition["name_relation_to_hgvs"]:
        for name in name_relation_to_missing_hgvs:
            if name == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
<<<<<<< HEAD
                        relation["hgvs"].remove(missing_pos)
                if not relation["hgvs"]:
                    for inx, guide_dip in enumerate(allele_matcher["guide_dip"]):
=======
                        #already remove but still looping could we just break
                        relation["hgvs"].remove(missing_pos)
                if not relation["hgvs"]:
                    for inx,guide_dip in enumerate(allele_matcher["guide_dip"]):
                        # if allele_definition["gene"] == "CYP2C9" : print("inx : ",inx," with",guide_dip)
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
                        if name in guide_dip.split('/'):
                            allele_matcher["guide_dip"].remove(guide_dip)
                    for print_dip in allele_matcher["print_dip"]:
                        if name in print_dip.split('/'):
                            allele_matcher["print_dip"].remove(print_dip)

<<<<<<< HEAD
    assert len(allele_matcher["guide_dip"]) == len(allele_matcher["print_dip"])
=======
    # assert len(allele_matcher["guide_dip"]) == len(allele_matcher["print_dip"])
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
        name1, name2 = guide_dip.split("/")
        score1 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name1 == relation["name"])
        score2 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name2 == relation["name"])
        guide_dip_score.append(score1 + score2)
    for print_dip in allele_matcher["print_dip"]:
        name1, name2 = print_dip.split("/")
=======
        name1,name2 = guide_dip.split("/")
        score1 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name1 == relation["name"])
        score2 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name2 == relation["name"])
        guide_dip_score.append(score1 + score2)

    for print_dip in allele_matcher["print_dip"]:
        name1,name2 = print_dip.split("/")
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
        score1 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name1 == relation["name"])
        score2 = next(len(relation["hgvs"]) for relation in allele_definition["name_relation_to_hgvs"] if name2 == relation["name"])
        print_dip_score.append(score1 + score2)

    guide_dip_new = []
    print_dip_new = []
<<<<<<< HEAD
    for index in locate_the_most(guide_dip_score)[1]:
        guide_dip_new.append(allele_matcher["guide_dip"][index])
    for index in locate_the_most(print_dip_score)[1]:
        print_dip_new.append(allele_matcher["print_dip"][index])

    assert len(allele_matcher["guide_dip"]) == len(allele_matcher["print_dip"])
    allele_matcher["count_diplotype"] = len(guide_dip_new)
=======
    for index in locate_max(guide_dip_score)[1]:
        guide_dip_new.append(allele_matcher["guide_dip"][index])
    for index in locate_max(print_dip_score)[1]:
        print_dip_new.append(allele_matcher["print_dip"][index])

    assert len(allele_matcher["guide_dip"]) == len(allele_matcher["print_dip"])
    allele_matcher["count_diplotype"] = len(print_dip_new)
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    allele_matcher["guide_dip"] = guide_dip_new
    allele_matcher["print_dip"] = print_dip_new

    return allele_matcher
