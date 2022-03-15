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
from pharmvip_guideline.allele_matcher.annotation import *

class  MyParser(argparse.ArgumentParser):
    def  error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    # begin to setting things up
    text1 = 'This python script create to process on guide_line module.It have 2 fuctional which is \n'
    text2 = '#1 allele_definitions_transform will convert those *.xls(those file store allele definitons) in to json file\n'
    text3 = '#2 allele_matcher will matching genome from VCF with suitable drug guideline '
    # dest_text = textwrap.dedent(text1+text2+text3)
    dest_text = text1+text2+text3
    parser = MyParser(description=dest_text)

    subparsers = parser.add_subparsers(dest="subparser_name")
 
    ## add parser for using transform fucntion (initiate table n stuff)
    allele_definitions_transform_parser = subparsers.add_parser(name="allele_definitions_transform")

    #add argrument of transform_parser for input path (directory where input stay eg. *.xls file)
    allele_definitions_transform_parser.add_argument(
        "--allele_definitions",
        # dest='$input_path',
        help='use this option follow with $path to specific where those *xls which store allele definitions to be use in the process ',
        required=False,
        default=defaults_allele_definitions_table # ~/resource/allele_definitions/$version_name/table
    )
    #add argrument of transform_parser for output path (directory where we write json file from the process)
    allele_definitions_transform_parser.add_argument(
        "--outputs",
        # dest='$output_path',
        help="use this option follow with $path to specific where json of allele definition from the process should be",
        required=False,
        default=defaults_allele_definitions_transform # ~/resource/allele_definitions/$version_name/transform/.....
    )
    allele_definitions_transform_parser.add_argument(
        "--dbpmcgenomics",
        help="use this option follow with path to specific where the tuple set of text should be write down",
        required=False,
        default=defaults_allele_definitions_dbpmcgenomics
    )

    #add parser to use matcher (which will find guidline)
    allele_matcher_parser = subparsers.add_parser("allele_matcher")
    allele_matcher_parser.add_argument(
        "--allele_definitions",
        help='define path to directory which store allele json path',
        required=False,
        default=defaults_allele_definitions_transform # ~/resource/allele_definitions/$version_name/transform/.....
    )
    allele_matcher_parser.add_argument(
        "--function_mappings",
        help='define path to directory which store allele fucntion',
        required=False,
        default=defaults_function_mappings # ~/resource/fuction_mapping/$version_name/.....
    )
    allele_matcher_parser.add_argument(
        "--clinical_guideline_annotations",
        help='define path to directory which store guideline annotations',
        required=False,
        default=defaults_clinical_guideline_annotations # ~/resource/clinical_guideline_anotations/$version_name/.....
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
    # print(vars(args))
    # print(args)

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

            summary_and_full_report = annotation(args.clinical_guideline_annotations, args.function_mappings, diplotype_cpic, f"{args.clinical_guideline_annotations}/annotations_short/guideline_add_short.xlsx")
            summary_and_full_report = handle_summary_and_full_report_layout(summary_and_full_report, diplotype_cpic)
            summary_and_full_report = drop_columns(summary_and_full_report, args.ana_options_cpic, args.ana_options_hla, args.ana_genes_cyp2d6)
            summary_and_full_report = summary_and_full_report.sort_values(by=['cpi_sum_gene1', 'cpi_sum_gene2', 'cpi_sum_gene3', 'cpi_sum_drug'])
            summary_and_full_report = summary_and_full_report.reset_index(drop = True)
            summary_and_full_report = replace_na(summary_and_full_report)
            summary_and_full_report = replace_matabolizer(summary_and_full_report)
            to_txt(summary_and_full_report, args.dbpmcgenomics, args.ana_user_id, args.ana_id)
       
        elif args.ana_options_cpic == "false" and args.ana_options_hla == "true":
            diplotype_hla = read_diplotype(args.diplotype_hla)

            summary_and_full_report = annotation(args.clinical_guideline_annotations, args.function_mappings, diplotype_hla, f"{args.clinical_guideline_annotations}/annotations_short/guideline_add_short.xlsx")
            summary_and_full_report = summary_and_full_report.sort_values(by=['cpi_sum_gene1', 'cpi_sum_gene2', 'cpi_sum_gene3', 'cpi_sum_drug'])
            summary_and_full_report = summary_and_full_report.reset_index(drop = True)
            summary_and_full_report = replace_na(summary_and_full_report)
            summary_and_full_report = replace_matabolizer(summary_and_full_report)
            to_txt(summary_and_full_report, args.dbpmcgenomics, args.ana_user_id, args.ana_id)
        else:
            print(f"error with ana_options: {args.ana_options_cpic}, {args.ana_genes_cyp2d6}, {args.ana_options_hla}")
            exit()

        print(f"run pharmvip_guideline allele_matcher successfully in {time.time() - allele_matcher_start_time:.2f} seconds")
    else :
        parser.print_help()