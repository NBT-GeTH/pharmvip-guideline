

def cftr_exception_handler(allele_matcher):
    for i in range(len(allele_matcher["print_dip"])):
        if allele_matcher["print_dip"][i] != "No info" and allele_matcher["print_dip"][i] != "?/?":
            incidental = ["G542X", "N1303K", "W1282X", "R553X", "1717-1G->A", "621+1G->T", "2789+5G->A", "3849+10kbC->T", "R1162X", "G85E", "3120+1G->A", "I507", "1898+1G->A", "3659delC", "R347P", "R560T", "R334W", "A455E", "2184delA", "711+1G->T"]
            if allele_matcher["print_dip"][i].split("/")[0] in incidental and allele_matcher["print_dip"][i].split("/")[1] in incidental:
                    allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[0]} (homozygous)"
            elif allele_matcher["print_dip"][i] == "Reference/Reference":
                allele_matcher["print_dip"][i] = "No CPIC variants found"
            elif allele_matcher["print_dip"][i].split("/")[0] == "Reference" or allele_matcher["print_dip"][i].split("/")[1] == "Reference":
                if allele_matcher["print_dip"][i].split("/")[0] == "Reference":
                    allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[1]} (heterozygous)"
                elif allele_matcher["print_dip"][i].split("/")[1] == "Reference":
                    allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[0]} (heterozygous)"
    