import unittest

import pandas as pd
from pharmvip_guideline.allele_definitions_transform.rsid import all_rsid, findall_rsid

def test_all_rsid():
    finall_rsid = ["rs137852347", "rs76723693"]
    actual = all_rsid(finall_rsid)
    expect = "rs137852347, rs76723693"
    assert actual == expect

def test_findall_rsid():
    rsid_cell = pd.Series([";;rs137852347;rs76723693;"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs137852347, rs76723693"]
    assert actual == expect

    rsid_cell = pd.Series(["rs1800559"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs1800559"]
    assert actual == expect

    rsid_cell = pd.Series([";;rs137852347;rs76723693;", "rs1800559"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs137852347, rs76723693", "rs1800559"]
    assert actual == expect
