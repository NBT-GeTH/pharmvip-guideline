
def dpyd_exception_handler(allele_matcher):
    for i in range(len(allele_matcher["print_dip"])):
        if allele_matcher["print_dip"][i] == "Reference/Reference":
            allele_matcher["print_dip"][i] = "No CPIC decreased or no function variant with strong or moderate evidence found"
