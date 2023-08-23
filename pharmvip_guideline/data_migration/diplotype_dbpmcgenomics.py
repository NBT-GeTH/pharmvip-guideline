def diplotype_gene_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_diplotype_gene.txt", "w")
    text = ""
    for index, row in diplotype_cpic.iterrows():
        text += f"{ana_user_id}\t{ana_id}\t{row['sample_id']}\t{row['gene']}\t{row['missing_call_variants']}\t{row['total_variants']}\t{row['gene_phases']}\t{row['count_diplotype']}\t{row['guide_dip']}\t{row['guide_dip']}\t{', '.join(row['guide_dip'])}\n"
    f.write(text[:-1])
    f.close()

def diplotype_matched_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_diplotype_matched.txt", "w")
    text = ""
    for index, row in diplotype_cpic.iterrows():
        for i in range(len(row["guide_dip"])):
            a_tag = ""
            if row["gene"] != "CYP2D6" :
                # a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i]}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>' #errr

                if row["gene"] == "SLCO1B1" and (len(row["guide_dip"]) > len(row["guide_dip"])):
                    if row["guide_dip"][i] == "rs4149056T/rs4149056T":
                        a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else "*1A/*1A"}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
                    elif row["guide_dip"][i] == "rs4149056C/rs4149056C":
                        a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else "*5/*5"}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
                    elif row["guide_dip"][i] == "rs4149056T/rs4149056C":
                        a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else "*1A/*5"}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
                    elif row["guide_dip"][i] == "rs4149056C/rs4149056T":
                        a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else "*5/*1A"}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
                    else:
                        a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else row["guide_dip"][i]}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
                else:
                    a_tag += f'<a href="/analysis/report_diplotype_details?pk={ana_id}&gene={row["gene"]}&matched_diplotypes={row["guide_dip"][i] if row["guide_dip"][i] == "No info" or row["guide_dip"][i] == "?/?" else row["guide_dip"][i]}" target="_blank"><span style="font-size: 12px;">{row["guide_dip"][i]}</span></a>'
            else:
                a_tag += f'{row["guide_dip"][i]}'
            if i != int(len(row["guide_dip"]) - 1):
                text+= f"{ana_user_id}\t{ana_id}\t{row['gene']}\t{a_tag},\n"
            else:
                text+= f"{ana_user_id}\t{ana_id}\t{row['gene']}\t{a_tag}\n"
    f.write(text[:-1])
    f.close()

def diplotype_gene_detail_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_diplotype_gene_detail.txt", "w")
    text = ""
    for index, row in diplotype_cpic.iterrows():
        assert len(row["dp"]) == len(row["gt_bases"]) == len(row["gt_phases"])
        for i in range(len(row["dp"])):
            text += f"{ana_user_id}\t{ana_id}\t{row['gene']}\t{row['gt_bases'][i]}\t{row['dp'][i]}\t{row['gt_phases'][i]}\n"
    f.write(text[:-1])
    f.close()

def diplotype_gene_guide_diplotype_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_diplotype_gene_print_diplotype.txt", "w")
    text = ""
    for index, row in diplotype_cpic.iterrows():
        for i in range(len(row["guide_dip"])):
            if "/" in row["guide_dip"][i]:
                for name in row["guide_dip"][i].split("/"):
                    text += f"{ana_user_id}\t{ana_id}\t{row['gene']}\t{name}\n"
            else:
                text += f"{ana_user_id}\t{ana_id}\t{row['gene']}\t{row['guide_dip'][i]}\n"
    f.write(text[:-1])
    f.close()

def diplotype_dbpmcgenomics(ana_user_id, ana_id, diplotype_cpic, diplotype_cyp2d6, dbpmcgenomics):
    diplotype_cpic = diplotype_cpic.append(diplotype_cyp2d6)
    diplotype_cpic = diplotype_cpic.sort_values(by=["gene"])
    diplotype_cpic = diplotype_cpic.reset_index(drop=True)

    diplotype_gene_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    diplotype_matched_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    diplotype_gene_detail_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    diplotype_gene_guide_diplotype_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
