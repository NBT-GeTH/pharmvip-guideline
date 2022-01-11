import unittest
from pharmvip_guideline.allele_matcher.query import *

def test_check_null_dp():
    assert check_null_dp(None) == 0
    assert check_null_dp(None) != 1
    assert check_null_dp(0) == 0
    assert check_null_dp(1) == 1
    assert check_null_dp(1.1) == 1.1
    assert check_null_dp("1.11") == "1.11"

def test_check_gt_format():
    assert check_gt_format("A/A") == True
    assert check_gt_format("A/ATCG") == True
    assert check_gt_format("ATCG/A") == True
    assert check_gt_format("A") == False
    assert not check_gt_format("A") == True

def test_extract_genotype():
    assert extract_genotype("CACNA1S", "A/T") == ("A", "T")
    assert extract_genotype("CACNA1S", "C|G") == ("C", "G")
    assert extract_genotype("CACNA1S", "A/.") == ("A", ".")
    assert extract_genotype("CACNA1S", ".|G") == (".", "G")
    assert extract_genotype("CACNA1S", "./.") == (".", ".")
    assert extract_genotype("CACNA1S", ".|.") == (".", ".")

    assert extract_genotype("G6PD", "A/T") == ("A", "T")
    assert extract_genotype("G6PD", "C|G") == ("C", "G")
    assert extract_genotype("G6PD", "A/.") == ("A", ".")
    assert extract_genotype("G6PD", ".|G") == (".", "G")
    assert extract_genotype("G6PD", "./.") == (".", ".")
    assert extract_genotype("G6PD", ".|.") == (".", ".")
    assert extract_genotype("G6PD", "A") == ("A", "A")
    assert extract_genotype("G6PD", ".") == (".", ".")

def test_extract_genotype_in_rage():
    assert extract_genotype_in_rage("CACNA1S", ["A/T", "./.", "C|G", ".|."]) == ".A.C.|.T.G."
    assert extract_genotype("CACNA1S", extract_genotype_in_rage("CACNA1S", ["A/.", "./.", ".|G", ".|."])) == (".A...", "...G.")

    assert extract_genotype_in_rage("G6PD", ["A/T", "./.", "C|G", ".|.", "A", "."]) == ".A.C.A.|.T.G.A."
    assert extract_genotype("G6PD", extract_genotype_in_rage("G6PD", ["./T", "./.", "C|.", ".|.", "A", "."])) == ("...C.A.", ".T...A.")

def test_convert_ins():
    assert convert_ins("A", "A") == "del"
    assert convert_ins("A", "ATCG") == "insTCG"

def test_convert_del():
    assert convert_del("ATCG", "ATCG") == "TCG"
    assert convert_del("ATCG", "A") == "delTCG"

def test_convert_allele():
    assert convert_allele("SNP", "snp", False, "A", "A") == "A"
    assert convert_allele("SNP", ".", False, "A", ".") == "."
    assert convert_allele("SNP", "unknown", False, "A", "A") == "A"
    assert convert_allele("SNP", "indel", False, "A", "A") == "del"
    assert convert_allele("SNP", "indel", False, "A", "ATCG") == "insTCG"
    assert convert_allele("SNP", "indel", True, "ATCG", "ATCG") == "TCG"
    assert convert_allele("SNP", "indel", True, "ATCG", "A") == "delTCG"

    assert convert_allele("INS", "snp", False, "A", "A") == "A"
    assert convert_allele("INS", ".", False, "A", ".") == "."
    assert convert_allele("INS", "unknown", False, "A", "A") == "del"
    assert convert_allele("INS", "indel", False, "A", "A") == "del"
    assert convert_allele("INS", "indel", False, "A", "ATCG") == "insTCG"
    assert convert_allele("INS", "indel", True, "ATCG", "ATCG") == "TCG"
    assert convert_allele("INS", "indel", True, "ATCG", "A") == "delTCG"

    assert convert_allele("DEL", "snp", False, "A", "A") == "A"
    assert convert_allele("DEL", ".", False, "A", ".") == "."
    assert convert_allele("DEL", "unknown", False, "A", "A") == ""
    assert convert_allele("DEL", "indel", False, "A", "A") == "del"
    assert convert_allele("DEL", "indel", False, "A", "ATCG") == "insTCG"
    assert convert_allele("DEL", "indel", True, "ATCG", "ATCG") == "TCG"
    assert convert_allele("DEL", "indel", True, "ATCG", "A") == "delTCG"

    assert convert_allele("CNV", "snp", False, "A", "A") == "A"
    assert convert_allele("CNV", ".", False, "A", ".") == "."
    assert convert_allele("CNV", "unknown", False, "A", "A") == "A"
    assert convert_allele("CNV", "indel", False, "A", "A") == "A"
    assert convert_allele("CNV", "indel", False, "A", "ATCG") == "ATCG"
    assert convert_allele("CNV", "indel", True, "ATCG", "ATCG") == "ATCG"
    assert convert_allele("CNV", "indel", True, "ATCG", "A") == "A"

def test_sum_up_gt_phases():
    assert sum_up_gt_phases(True, ".", ".") == "."
    assert sum_up_gt_phases(False, "A", ".") == "."
    assert sum_up_gt_phases(True, ".", "A") == "."

    assert sum_up_gt_phases(True, "A", "A") == True
    assert sum_up_gt_phases(True, "A", "ATCG") == True
    assert sum_up_gt_phases(True, "ATCG", "A") == True

    assert sum_up_gt_phases(False, "A", "A") == True
    assert sum_up_gt_phases(False, "ATCG", "ATCG") == True
    assert sum_up_gt_phases(False, "A", "ATCG") == False
    assert sum_up_gt_phases(False, "ATCG", "A") == False

def test_count_call_variants():
    assert count_call_variants([{"gt_phases": "."}, {"gt_phases": "."}, {"gt_phases": "."}]) == 0
    assert count_call_variants([{"gt_phases": True}, {"gt_phases": False}, {"gt_phases": "."}]) == 2

def test_count_missing_call_variants():
    assert count_missing_call_variants([{"gt_phases": "."}, {"gt_phases": "."}, {"gt_phases": "."}]) == 3
    assert count_missing_call_variants([{"gt_phases": "."}, {"gt_phases": "."}, {"gt_phases": "Combine"}]) == 2

def test_sum_up_gene_phases():
    assert sum_up_gene_phases([{"gt_phases": "."}, {"gt_phases": "."}, {"gt_phases": "."}]) == "."

    assert sum_up_gene_phases([{"gt_phases": True}, {"gt_phases": "."}, {"gt_phases": "."}]) == True
    assert sum_up_gene_phases([{"gt_phases": True}, {"gt_phases": True}, {"gt_phases": True}]) == True

    assert sum_up_gene_phases([{"gt_phases": False}, {"gt_phases": "."}, {"gt_phases": "."}]) == False
    assert sum_up_gene_phases([{"gt_phases": False}, {"gt_phases": False}, {"gt_phases": False}]) == False

    assert sum_up_gene_phases([{"gt_phases": True}, {"gt_phases": False}, {"gt_phases": "."}]) == "Combine"
    assert sum_up_gene_phases([{"gt_phases": "."}, {"gt_phases": False}, {"gt_phases": True}]) == "Combine"

def test_query():
    return

if __name__ == "__main__":
    test_check_null_dp()
    test_check_gt_format()
    test_extract_genotype()
    test_extract_genotype_in_rage()
    test_convert_ins()
    test_convert_del()
    test_convert_allele()
    test_sum_up_gt_phases()
    test_count_call_variants()
    test_count_missing_call_variants()
    test_sum_up_gene_phases()
    test_query()

    print("test_query pass")
