#!/home/ferrari/anaconda3/bin/python
# branch "select_type"

#todo: check why, if I run the script twice, the files don't get updated

import argparse
import sys
import os
import textwrap


def parse_args(defaults={"verbose":False,
                         "feature_to_extract":["genes","TSS"],
                         "from_where":"gene",
                         "protein_coding":False,
                         "out_dir":"./",
                         "before_tss":1000,
                         "after_tss":500,
                         "before_gene_start_site":0,
                         "after_gene_end_site":0,
                         "before_tes":500,
                         "after_tes":500,
                         "percentile_range":None}):

    """parse arguments from the command line"""

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(__description__),
        add_help= False
    )


    # positional - required
    parser.add_argument("gtf_file",
                       metavar="GTF",
                       help="gencode gtf file from which you want to extract genomic coordinates")
    # general arguments
    general = parser.add_argument_group('general arguments')
    general.add_argument("-o", "--output_directory",
                         dest="out_dir",
                         help="output directory",
                         default=defaults["out_dir"])
    general.add_argument("-h","--help",
                         action="help",
                         help="show this help message and exit")
    general.add_argument("-v","--verbose",
                         dest="verbose",
                         action="store_true",
                         help="verbose output (default: {})".format(defaults["verbose"]),
                         default=defaults["verbose"])
    general.add_argument("-f","--feature",
                         dest="FEATURE",
                         nargs="*",
                         choices=["genes","TSS","TES"],
                         default=defaults["feature_to_extract"])
    general.add_argument("-prot_cod","--protein_coding",
                         dest="PROTEIN_CODING",
                         action="store_true",
                         default=defaults["protein_coding"],
                         help="set this flag if you want to extract only protein coding genes")
    general.add_argument("-w","--from_what",
                         dest="from_what",
                         choices=["gene","transcript"],
                         default=defaults["from_where"],
                         help="Do you want to extract TSS/TES info from genes or transcripts?")
    general.add_argument("-b_tss","--before_tss",
                         dest="BEFORE_TSS",
                         type=int,
                         help="number of bp to include before the transcription start site",
                         default=defaults["before_tss"]
                         )
    general.add_argument("-a_tss","--after_tss",
                         dest="AFTER_TSS",
                         type=int,
                         help="number of bp to include after the transcription start site",
                         default=defaults["after_tss"]
                         )
    general.add_argument("-b_tes","--before_tes",
                         dest="BEFORE_TES",
                         type=int,
                         help="number of bp to include before the transcription end site",
                         default=defaults["before_tes"]
                         )
    general.add_argument("-a_tes","--after_tes",
                         dest="AFTER_TES",
                         type=int,
                         help="number of bp to include after the transcription end site",
                         default=defaults["after_tes"]
                         )
    general.add_argument("-b_gene","--before_gene",
                         dest="BEFORE_GENE",
                         type=int,
                         help="number of bp to include before the gene body start site",
                         default=defaults["before_gene_start_site"]
                         )
    general.add_argument("-a_gene","--after_gene",
                         dest="AFTER_GENE",
                         type=int,
                         help="number of bp to include after the gene body end site",
                         default=defaults["after_gene_end_site"]
                         )
    return parser





def info_to_dict(info):

    dict_info={}

    lista = info.split(";")
    lista = [i.strip() for i in lista]

    for j in lista:
        pair = j.split()
        if len(pair) == 2:
            key,value = pair[0],pair[1]
            if len(value.split('"')) > 1:
                value = value.split('"')[1]

            if not key in dict_info:
                dict_info[key] = value
            else:
                pass
                #print("Warning: duplicate keys have been found in the information section: {}".format(key))
        else:
            if len(pair) != 0:
                print("Warning: more than two values found: {}".format(pair))

    return dict_info





def line_to_dict(gtf_line_list):

    keys = ["chromosome_name",
            "annotation_source",
            "feature_type",
            "genomic_start_location",
            "genomic_end_location",
            "score",
            "genomic_strand",
            "genomic_phase_CDS",
            "additional_information"]

    line_dict = {}

    if len(gtf_line_list) == 9:
        for i in range(len(keys)):
            line_dict[keys[i]] = gtf_line_list[i]
        line_dict["additional_information"] = info_to_dict(line_dict["additional_information"])
    else:
        if not gtf_line_list[0].startswith("#"):
            print("Warning! The following line does not conform to the expected gencode gtf format:\n{}".format(" ".join(gtf_line_list)))

    return line_dict




def make_bed(list_dict, feature, from_where, arg):

    bed_out_list = []

    if feature == "genes" and list_dict["feature_type"] == from_where:
        if list_dict["genomic_strand"] == '+':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - arg.BEFORE_GENE),
                str(int(list_dict["genomic_end_location"]) + arg.AFTER_GENE),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])

        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - arg.AFTER_GENE),
                str(int(list_dict["genomic_end_location"]) + arg.BEFORE_GENE),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])



    elif feature == "TSS" and list_dict["feature_type"] == from_where:
        if list_dict["genomic_strand"] == '+':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - arg.BEFORE_TSS),
                str(int(list_dict["genomic_start_location"]) + arg.AFTER_TSS),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])



        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_end_location"]) - arg.AFTER_TSS),
                str(int(list_dict["genomic_end_location"]) + arg.BEFORE_TSS),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])



    elif feature == "TES" and list_dict["feature_type"] == from_where:
        if list_dict["genomic_strand"] == '+':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_end_location"]) - arg.BEFORE_TES),
                str(int(list_dict["genomic_end_location"]) + arg.AFTER_TES),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])


        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - arg.AFTER_TES),
                str(int(list_dict["genomic_start_location"]) + arg.BEFORE_TES),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"],
                list_dict["additional_information"]["gene_name"]
            ]
            if list_dict["feature_type"] == 'transcript':
                bed_out_list.append(list_dict["additional_information"]["transcript_id"])



    return bed_out_list


def main():
    l = 0
    parser = parse_args()
    arg = parser.parse_args()
    print(arg)
    with open(arg.gtf_file) as in_file:

        ref_dict = {}
        for feature_ in arg.FEATURE:
            ref_dict[feature_] = open("{}/{}.bed".format(arg.out_dir, feature_),'w')

        for line in in_file:
            lista = line.strip().split("\t")
            list_dict = line_to_dict(lista)
            if len(list_dict) > 0:
                if arg.PROTEIN_CODING:
                    if list_dict["additional_information"]["gene_type"] == "protein_coding":
                        for feature_ in arg.FEATURE:
                            printing_line = make_bed(list_dict, feature_, arg.from_what, arg)
                            if len(printing_line) > 0:
                                printing_line[-1] = printing_line[-1]+"\n"
                                ref_dict[feature_].write("\t".join(printing_line))

                elif not list_dict["chromosome_name"].endswith("MT") and not list_dict["chromosome_name"].endswith("M"):

                    for feature_ in arg.FEATURE:
                        printing_line = make_bed(list_dict, feature_, arg.from_what, arg)
                        if len(printing_line) > 0:
                            printing_line[-1] = printing_line[-1]+"\n"
                            ref_dict[feature_].write("\t".join(printing_line))




        for feature_ in arg.FEATURE:
            ref_dict[feature_].close()


if __name__ == "__main__" :
    main()
