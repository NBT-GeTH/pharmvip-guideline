import json

def create_json_file(ana_user_id, ana_id, diplotype_cpic, diplotype_cyp2d6, dbpmcgenomics):
    data = {}
    json.dump(data,dbpmcgenomics + '/analysis_cpic_summary_drug_view.json')
    json.dump(data,dbpmcgenomics + '/analysis_cpic_summary_gene_view.json')
    

    
    # diplotype_cpic = diplotype_cpic.append(diplotype_cyp2d6)
    # diplotype_cpic = diplotype_cpic.sort_values(by=["gene"])
    # diplotype_cpic = diplotype_cpic.reset_index(drop=True)

    # diplotype_gene_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    # diplotype_matched_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    # diplotype_gene_detail_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)
    # diplotype_gene_print_diplotype_text(ana_user_id, ana_id, diplotype_cpic, dbpmcgenomics)