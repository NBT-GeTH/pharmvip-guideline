import argparse
import sys

import time
from pharmvip_guideline import *
from pharmvip_guideline.allele_definitions_transform.transform import transform
from pharmvip_guideline.allele_definitions_transform.transform_dbpmcgenomics import transform_dbpmcgenomics
from pharmvip_guideline.allele_matcher.matcher import matcher
from pharmvip_guideline.allele_matcher.diplotype import create_diplotype_cpic, read_diplotype, read_hla
from cyvcf2 import VCF
from pharmvip_guideline.allele_matcher.diplotype_dbpmcgenomics import diplotype_dbpmcgenomics
from pharmvip_guideline.annotation.guideline_annotation import *
import numpy as np

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"error: {message}\n")
        self.print_help()
        sys.exit(2)

def replace_blank(df):
    # df = df.fillna('N/A')

    # df['cpi_sum_gene1'] = df['cpi_sum_gene1'].replace([''],'N/A')
    # df['cpi_sum_gene2'] = df['cpi_sum_gene2'].replace([''],'N/A')
    # df['cpi_sum_gene3'] = df['cpi_sum_gene3'].replace([''],'N/A')
    # df['cpi_sum_dip_name1'] = df['cpi_sum_dip_name1'].replace([''],'N/A')
    # df['cpi_sum_dip_name2'] = df['cpi_sum_dip_name2'].replace([''],'N/A')
    # df['cpi_sum_dip_name3'] = df['cpi_sum_dip_name3'].replace([''],'N/A')
    # df['cpi_sum_drug'] = df['cpi_sum_drug'].replace([''],'N/A')
    # df['cpi_sum_population'] = df['cpi_sum_population'].replace([''],'N/A')
    # df['cpi_sum_act_score1'] = df['cpi_sum_act_score1'].replace([''],'N/A')
    # df['cpi_sum_act_score2'] = df['cpi_sum_act_score2'].replace([''],'N/A')
    # df['cpi_sum_act_score3'] = df['cpi_sum_act_score3'].replace([''],'N/A')
    df['cpi_sum_strength'] = df['cpi_sum_strength'].replace([''],'N/A')
    df['cpi_sum_recommendations'] = df['cpi_sum_recommendations'].replace([''],'N/A')
    df['cpi_sum_recommendations_full'] = df['cpi_sum_recommendations_full'].replace([''],'N/A')
    # df['cpi_sum_recommendations_full_figure'] = df['cpi_sum_recommendations_full_figure'].replace([''],'N/A')
    # df['cpi_sum_comments'] = df['cpi_sum_comments'].replace([''],'N/A')
    df['cpi_sum_implications1'] = np.where((df.cpi_sum_gene1 != '') & (df.cpi_sum_implications1 == ''),'N/A',df.cpi_sum_implications1)
    df['cpi_sum_implications2'] = np.where((df.cpi_sum_gene2 != '') & (df.cpi_sum_implications2 == ''),'N/A',df.cpi_sum_implications2)
    df['cpi_sum_implications3'] = np.where((df.cpi_sum_gene3 != '') & (df.cpi_sum_implications3 == ''),'N/A',df.cpi_sum_implications3)
    df['cpi_sum_phenotype1'] = np.where((df.cpi_sum_gene1 != '') & (df.cpi_sum_phenotype1 == ''),'N/A',df.cpi_sum_phenotype1)
    df['cpi_sum_phenotype2'] = np.where((df.cpi_sum_gene2 != '') & (df.cpi_sum_phenotype2 == ''),'N/A',df.cpi_sum_phenotype2)
    df['cpi_sum_phenotype3'] = np.where((df.cpi_sum_gene3 != '') & (df.cpi_sum_phenotype3 == ''),'N/A',df.cpi_sum_phenotype3)
    df['cpi_sum_met_status_1'] = np.where((df.cpi_sum_gene1 != '') & (df.cpi_sum_met_status_1 == ''),'N/A',df.cpi_sum_met_status_1)
    df['cpi_sum_met_status_2'] = np.where((df.cpi_sum_gene2 != '') & (df.cpi_sum_met_status_2 == ''),'N/A',df.cpi_sum_met_status_2)
    df['cpi_sum_met_status_3'] = np.where((df.cpi_sum_gene3 != '') & (df.cpi_sum_met_status_3 == ''),'N/A',df.cpi_sum_met_status_3)
    # df['cpi_sum_gen_1_missing'] = df['cpi_sum_gen_1_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_1_total'] = df['cpi_sum_gen_1_total'].replace([''],'N/A')
    # df['cpi_sum_gen_2_missing'] = df['cpi_sum_gen_2_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_2_total'] = df['cpi_sum_gen_2_total'].replace([''],'N/A')
    # df['cpi_sum_gen_3_missing'] = df['cpi_sum_gen_3_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_3_total'] = df['cpi_sum_gen_3_total'].replace([''],'N/A')
    # df['cpi_sum_hla_tool_1_guide'] = df['cpi_sum_hla_tool_1_guide'].replace([''],'N/A')
    # df['cpi_sum_hla_tool_2_guide'] = df['cpi_sum_hla_tool_2_guide'].replace([''],'N/A')
    
    return df

def main():
    text1 = "This python script create to process on guideline module.It have 2 fuctional which is\n"
    text2 = "#1 allele_definitions_transform will convert those *.xlsx(those file store allele definitons) in to json file\n"
    text3 = "#2 allele_matcher will matching genome from VCF with suitable drug guideline"
    dest_text = text1+text2+text3
    parser = MyParser(description=dest_text)

    subparsers = parser.add_subparsers(dest="subparser_name")
 
    allele_definitions_transform_parser = subparsers.add_parser(name="allele_definitions_transform")
    allele_definitions_transform_parser.add_argument(
        "--allele_definitions",
        help="use this option follow with $path to specific where those *.xlsx which store allele definitions to be use in the process",
        required=False,
        default=defaults_allele_definitions_table
    )
    allele_definitions_transform_parser.add_argument(
        "--outputs",
        help="use this option follow with $path to specific where json of allele definition from the process should be",
        required=False,
        default=defaults_allele_definitions_transform
    )
    allele_definitions_transform_parser.add_argument(
        "--dbpmcgenomics",
        help="use this option follow with $path to specific where the tuple set of text should be write down",
        required=False,
        default=defaults_allele_definitions_dbpmcgenomics
    )

    allele_matcher_parser = subparsers.add_parser("allele_matcher")
    allele_matcher_parser.add_argument(
        "--allele_definitions",
        help="define path to directory which store allele json path",
        required=False,
        default=defaults_allele_definitions_transform
    )
    allele_matcher_parser.add_argument(
        "--function_mappings",
        help="define path to directory which store allele fucntion",
        required=False,
        default=defaults_function_mappings
    )
    allele_matcher_parser.add_argument(
        "--clinical_guideline_annotations",
        help="define path to directory which store guideline annotations",
        required=False,
        default=defaults_clinical_guideline_annotations
    )
    allele_matcher_parser.add_argument(
        "--ana_user_id",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--ana_id",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--ana_options_cpic",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--ana_best_candidate",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--ana_genes_cyp2d6",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--ana_options_hla",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--vcf_gz_file",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--diplotype_cyp2d6",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--diplotype_hla",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--outputs",
        help="file path to write json version of result files",
        required=True
    )
    allele_matcher_parser.add_argument(
        "--dbpmcgenomics",
        help="file path to write text version of result files",
        required=True
    )

    args = parser.parse_args()

    if args.subparser_name == "allele_definitions_transform":
        allele_definitions_transform_start_time = time.time()

        transform(
            args.allele_definitions,
            args.outputs
        )

        transform_dbpmcgenomics(
            args.outputs,
            args.dbpmcgenomics
        )

        print(f"run pharmvip_guideline allele_definitions_transform successfully in {time.time() - allele_definitions_transform_start_time:.2f} seconds")
        
    elif args.subparser_name == "allele_matcher":
        allele_matcher_start_time = time.time()

        if args.ana_options_cpic == "true":
            
            matcher(
                args.allele_definitions,
                args.ana_user_id,
                args.ana_id,
                args.ana_best_candidate,
                args.vcf_gz_file,
                args.outputs
            )

            diplotype_cpic = create_diplotype_cpic(args.outputs)
            
            diplotype_cyp2d6 = read_diplotype(args.diplotype_cyp2d6)
            diplotype_cyp2d6["sample_id"] = VCF(args.vcf_gz_file).samples[0]

            diplotype_dbpmcgenomics(args.ana_user_id, args.ana_id, diplotype_cpic, diplotype_cyp2d6, args.dbpmcgenomics)

            diplotype_hla = read_hla(args.diplotype_hla)

            diplotype_cpic = diplotype_cpic.append(diplotype_cyp2d6)
            diplotype_cpic = diplotype_cpic.append(diplotype_hla)
            diplotype_cpic = diplotype_cpic.sort_values(by=["gene"])
            diplotype_cpic = diplotype_cpic.reset_index(drop=True)
            
            summary_and_full_report = annotate(args.clinical_guideline_annotations, args.function_mappings, diplotype_cpic)
            summary_and_full_report = summary_and_full_report.sort_values(by=['cpi_sum_gene1', 'cpi_sum_gene2', 'cpi_sum_gene3', 'cpi_sum_drug'])
            summary_and_full_report = summary_and_full_report.reset_index(drop = True)
            summary_and_full_report = replace_blank(summary_and_full_report)
            export_guideline_report(summary_and_full_report, args.dbpmcgenomics, args.ana_user_id, args.ana_id)
       
        elif args.ana_options_cpic == "false" and args.ana_options_hla == "true":
            diplotype_hla = read_diplotype(args.diplotype_hla)

            summary_and_full_report = annotate(args.clinical_guideline_annotations, args.function_mappings, diplotype_hla, f"{args.clinical_guideline_annotations}/annotations_short/guideline_add_short.xlsx")
            summary_and_full_report = summary_and_full_report.sort_values(by=['cpi_sum_gene1', 'cpi_sum_gene2', 'cpi_sum_gene3', 'cpi_sum_drug'])
            summary_and_full_report = summary_and_full_report.reset_index(drop = True)
            summary_and_full_report = replace_blank(summary_and_full_report)
            export_guideline_report(summary_and_full_report, args.dbpmcgenomics, args.ana_user_id, args.ana_id)

        else:
            print(f"error with ana_options: {args.ana_options_cpic}, {args.ana_genes_cyp2d6}, {args.ana_options_hla}")
            exit()

        print(f"run pharmvip_guideline allele_matcher successfully in {time.time() - allele_matcher_start_time:.2f} seconds")
    else:
        parser.print_help()
