import re

def manual_customize(allele_definition_df, gene):
    if gene == "G6PD":
        allele_definition_df.iloc[1, 0] = allele_definition_df.iloc[1, 1]
        allele_definition_df.iloc[2, 0] = allele_definition_df.iloc[2, 1]
        allele_definition_df.iloc[3, 0] = allele_definition_df.iloc[3, 1]
        allele_definition_df.iloc[4, 0] = allele_definition_df.iloc[4, 1]
        allele_definition_df.iloc[5, 0] = allele_definition_df.iloc[5, 1]

        allele_definition_df.drop([1], axis=1, inplace=True)
        allele_definition_df.columns = range(allele_definition_df.shape[1])

        return allele_definition_df
    elif gene == "VKORC1":
        allele_definition_df.iloc[7, 0] = "-1639A"
        allele_definition_df.iloc[8, 0] = "-1639G"
        
        return allele_definition_df
    else:
        return allele_definition_df

def automatic_customize(allele_definition_df):
    allele_definition_df = sort_hgvs(allele_definition_df)

    allele_definition_df = allele_definition_df.iloc[:find_end_row(allele_definition_df), :find_end_col(allele_definition_df)]

    return allele_definition_df

def sort_hgvs(allele_definition_df):
    allele_definition_df.sort_values(axis=1, by=3, inplace=True)
    allele_definition_df.columns = range(allele_definition_df.shape[1])
    return allele_definition_df

def find_end_row(allele_definition_df):
    for index, cell in allele_definition_df.iloc[:, 0].items():
        if str(cell) == "nan" or re.match(r"^NOTES:\s*$", str(cell)):
            return index

def find_end_col(allele_definition_df):
    for index, cell in allele_definition_df.iloc[3, :].items():
        if str(cell) == "nan":
            return index
