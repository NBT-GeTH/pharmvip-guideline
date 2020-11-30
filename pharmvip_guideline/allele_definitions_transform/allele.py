import copy

def clean_allele_cell(allele_definition_df):
    allele_definition_df = clean_nan_allele_cell(allele_definition_df)

    allele_definition_df = clean_whitespace_allele_cell(allele_definition_df)

    return allele_definition_df

def clean_nan_allele_cell(allele_definition_df):
    allele_definition_df_ = copy.deepcopy(allele_definition_df)
    for i in range(allele_definition_df_.shape[0]):
        if i >= 7:
            for j in range(allele_definition_df_.shape[1]):
                if j >= 1:
                    if str(allele_definition_df_.iloc[i, j]) == "nan":
                        allele_definition_df_.iloc[i, j] = allele_definition_df_.iloc[7, j]
    return allele_definition_df_

def clean_whitespace_allele_cell(allele_definition_df):
    allele_definition_df_ = copy.deepcopy(allele_definition_df)
    for i in range(allele_definition_df_.shape[0]):
        if i >= 7:
            for j in range(allele_definition_df_.shape[1]):
                if j >= 1:
                    if " " in str(allele_definition_df_.iloc[i, j]):
                        allele_definition_df_.iloc[i, j] = allele_definition_df_.iloc[i, j].replace(" ", "")
    return allele_definition_df_
