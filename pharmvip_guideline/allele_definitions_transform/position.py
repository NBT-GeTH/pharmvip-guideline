import re

def search_chromosome(position_cell):
    search_chromosome = re.search(r"chromosome\s(\w+)", str(position_cell))
    if search_chromosome:
        """
        position_cell                                                       search_chromosome
        'Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)'   1. '1'
        'Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)'   1. '7'
        """
        return position_cell, "chr" + search_chromosome.group(1)
    else:
        print(f"error search chromosome with: {str(position_cell)}")
        exit()
