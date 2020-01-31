#!/home/ferrari/anaconda3/bin/python
import argparse
import sys
import os
import math


def parse_args(defaults={"verbose":False,
                         "feature_to_extract":["gene_body","TSS"],
                         "from_where":"gene",
                         "out_dir":"./",
                         "before_tss":1000,
                         "after_tss":500,
                         "before_gene_start_site":0,
                         "after_gene_end_site":0,
                         "before_tes":500,
                         "after_tes":500,
                         "percentile_range":None,
                         "limit":False}):

    """parse arguments from the command line"""

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(__description__),
        add_help= False
    )


    # positional - required
    parser.add_argument("bed_file",
                       metavar="bed",
                       help="bed file from which you want to extract genomic coordinates")
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
                         choices=["gene_body","TSS","TES"],
                         default=defaults["feature_to_extract"],
                         help="the feature that you want to extract form the gtf file")
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
    general.add_argument("-l","--limit",
                         action="store_true",
                         help="set this flag if you want to limit the window to at most the feature length",
                         default=defaults["limit"]
                         )
    return parser

def is_there_header(bed_file):
    with open(bed_file) as in_file:
        line = in_file.readline().strip().split()
        try:
            int(line[1])
            return False
        except:
            return True

def make_bed(list_bed, feature, arg):

    if feature == "gene_body":
        if list_bed[5] == '+':
            bed_out_list = [
                list_bed[0],
                str(int(list_bed[1]) - arg.BEFORE_GENE),
                str(int(list_bed[2]) + arg.AFTER_GENE),
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs((int(list_bed[1]) - arg.BEFORE_GENE) - (int(list_bed[2]) + arg.AFTER_GENE)))
            ]

        elif list_bed[5] == '-':
            bed_out_list = [
                list_bed[0],
                str(int(list_bed[1]) - arg.AFTER_GENE),
                str(int(list_bed[2]) + arg.BEFORE_GENE),
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs((int(list_bed[1]) - arg.AFTER_GENE) - (int(list_bed[2]) + arg.BEFORE_GENE)))
            ]

    elif feature == "TSS":
        if list_bed[5] == '+':
            if arg.limit:
                max_length = math.fabs((int(list_bed[1]) - int(list_bed[2])))
                if max_length < arg.AFTER_TSS:
                    start = str(int(int(list_bed[1]) - arg.BEFORE_TSS))
                    end = str(int(int(list_bed[1]) + max_length))
                else:
                    start = str(int(int(list_bed[1]) - arg.BEFORE_TSS))
                    end = str(int(int(list_bed[1]) + arg.AFTER_TSS))
            else:
                start = str(int(int(list_bed[1]) - arg.BEFORE_TSS))
                end = str(int(int(list_bed[1]) + arg.AFTER_TSS))

            bed_out_list = [
                list_bed[0],
                start,
                end,
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs(int(end) - int(start)))
                ]

        elif list_bed[5] == '-':
            if arg.limit:
                max_length = math.fabs((int(list_bed[1]) - int(list_bed[2])))
                if max_length < arg.AFTER_TSS:
                    start = str(int(int(list_bed[2]) - max_length))
                    end = str(int(int(list_bed[2]) + arg.BEFORE_TSS))
                else:
                    start = str(int(int(list_bed[2]) - arg.AFTER_TSS))
                    end = str(int(int(list_bed[2]) + arg.BEFORE_TSS))
            else:
                start = str(int(int(list_bed[2]) - arg.AFTER_TSS))
                end = str(int(int(list_bed[2]) + arg.BEFORE_TSS))

            bed_out_list = [
                list_bed[0],
                start,
                end,
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs(int(start) - int(end)))
            ]

    elif feature == "TES" :
        if list_bed[5] == '+':
            if arg.limit:
                max_length = math.fabs((int(list_bed[1]) - int(list_bed[2])))
                if max_length < arg.BEFORE_TES:
                    start = str(int(int(list_bed[2]) - max_length))
                    end = str(int(int(list_bed[2]) + arg.AFTER_TES))
                else:
                    start = str(int(int(list_bed[2]) - arg.BEFORE_TES))
                    end = str(int(int(list_bed[2]) + arg.AFTER_TES))
            else:
                start = str(int(int(list_bed[2]) - arg.BEFORE_TES))
                end = str(int(int(list_bed[2]) + arg.AFTER_TES))

            bed_out_list = [
                list_bed[0],
                start,
                end,
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs(int(start) - int(end)))
            ]

        elif list_bed[5] == '-':
            if arg.limit:
                max_length = math.fabs((int(list_bed[1]) - int(list_bed[2])))
                if max_length < arg.BEFORE_TES:
                    start = str(int(int(list_bed[1]) - arg.AFTER_TES))
                    end = str(int(int(list_bed[1]) + max_length))
                else:
                    start = str(int(int(list_bed[1]) - arg.AFTER_TES))
                    end = str(int(int(list_bed[1]) + arg.BEFORE_TES))
            else:
                start = str(int(int(list_bed[1]) - arg.AFTER_TES))
                end = str(int(int(list_bed[1]) + arg.BEFORE_TES))

            bed_out_list = [
                list_bed[0],
                start,
                end,
                list_bed[3],
                list_bed[4],
                list_bed[5],
                list_bed[6],
                str(math.fabs(int(start) - int(end)))
            ]

    return bed_out_list


def main():
    parser = parse_args()
    arg = parser.parse_args()
    print(arg)
    kount = 0
    with open(arg.bed_file) as in_file:

        ref_dict = {}
        for feature_ in arg.FEATURE:
            ref_dict[feature_] = open("{}.bed".format(feature_),'w')

        if is_there_header(arg.bed_file):
            kount-=1

        for line in in_file:
            kount += 1
            if kount >= 1:
                lista = line.strip().split("\t")
                for feature_ in arg.FEATURE:
                    printing_line = make_bed(lista, feature_, arg)
                    if len(printing_line) > 0:
                        printing_line[-1] = printing_line[-1]+"\n"
                        ref_dict[feature_].write("\t".join(printing_line))
                    else:
                        print("error")
                        break
    for feature_ in arg.FEATURE:
        ref_dict[feature_].close()

if __name__ == "__main__":
    main()
