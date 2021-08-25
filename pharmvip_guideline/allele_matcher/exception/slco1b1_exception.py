

def slco1b1_exception_handler(allele_matcher):
    if allele_matcher["guide_dip"] == ["?/?"] and allele_matcher["print_dip"] == ["?/?"]:
        for variant in allele_matcher["variants"]:
            if variant["rsid"] == "rs4149056":
                if variant["gt_bases"] == "T/T" or variant["gt_bases"] == "T|T":
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*1A/*1A"]
                    allele_matcher["print_dip"] = ["rs4149056T/rs4149056T"]
                elif variant["gt_bases"] == "C/C" or variant["gt_bases"] == "C|C":
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*5/*5"]
                    allele_matcher["print_dip"] = ["rs4149056C/rs4149056C"]
                elif variant["gt_bases"] == "T/C" or variant["gt_bases"] == "T|C":
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*1A/*5"]
                    allele_matcher["print_dip"] = ["rs4149056T/rs4149056C"]
                elif variant["gt_bases"] == "C/T" or variant["gt_bases"] == "C|T":
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*5/*1A"]
                    allele_matcher["print_dip"] = ["rs4149056C/rs4149056T"]
    else:
        for variant in allele_matcher["variants"]:
            if variant["rsid"] == "rs4149056":
                if variant["gt_bases"] == "T/T" or variant["gt_bases"] == "T|T":
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*1A/*1A")
                    allele_matcher["print_dip"].append("rs4149056T/rs4149056T")
                elif variant["gt_bases"] == "C/C" or variant["gt_bases"] == "C|C":
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*5/*5")
                    allele_matcher["print_dip"].append("rs4149056C/rs4149056C")
                elif variant["gt_bases"] == "T/C" or variant["gt_bases"] == "T|C":
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*1A/*5")
                    allele_matcher["print_dip"].append("rs4149056T/rs4149056C")
                elif variant["gt_bases"] == "C/T" or variant["gt_bases"] == "C|T":
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*5/*1A")
                    allele_matcher["print_dip"].append("rs4149056C/rs4149056T")    