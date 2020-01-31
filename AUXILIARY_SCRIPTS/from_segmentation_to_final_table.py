#!/home/ferrari/anaconda3/bin/python

import argparse
import subprocess
import pandas as pd
import sys
import os

def compute_overlaps(args, temp_dir):

    """compute overlap between the regions of intereste provided and the segmentation"""

    for r in args.regions_of_interest:
        bed_name = ".".join(os.path.basename(r).split(".")[:-1])
        # print(bed_name)
        with open(os.path.join(temp_dir,"intersect_{}_with_segmentation.txt".format(bed_name)),"w") as out_file:
             subprocess.run(["bedtools","intersect","-a",r,"-b",args.segmentation_file,"-wao"], stdout=out_file)

def get_states(segmentation_f):

    p1 = subprocess.Popen("cut -f 4 {}".format(segmentation_f).split(), stdout=subprocess.PIPE)
    p2 = subprocess.Popen("sort", stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen("uniq", stdin=p2.stdout, stdout=subprocess.PIPE)

    return [i.decode() for i in p3.communicate()[0].split()]

def build_df_overlaps(intersect_file,args):
    id_file = ".".join(os.path.basename(intersect_file).split(".")[:-1])
    with open(intersect_file) as inter:
        dixio=dict()
        states = get_states(args.segmentation_file)
        i=0
        for line in inter:
            lista = line.strip().split()
            feat_length=int(lista[2])-int(lista[1])
            gene_ID = lista[3]
            if not gene_ID in dixio:
                dixio[gene_ID]=dict()
                dixio[gene_ID]['feature_length']=feat_length
                for s in [p+"_"+id_file for p in states]:
                    dixio[gene_ID][s]=[]
                if lista[-2].startswith("E"):
                    dixio[gene_ID][lista[-2]+"_"+id_file].append(int(lista[-1]))
            else:
                if lista[-2].startswith("E"):
                    dixio[gene_ID][lista[-2]+"_"+id_file].append(int(lista[-1]))

    for k in dixio:
        for s in dixio[k]:
            if s != 'feature_length':
                dixio[k][s]=round((sum(dixio[k][s])/dixio[k]['feature_length']),3)

    seg_df = pd.DataFrame.from_dict(dixio).T
    seg_df = seg_df.drop(columns = "feature_length")

    return seg_df

def main():

    #create input options
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--segmentation_file", help="segmentation file in bed format",
                            required=True)
    parser.add_argument("-r","--regions_of_interest",
                            nargs='*',
                            help="bed files with genomic features you want to use for the creation of the tSNE map \
                                  (output of 'extract_from_gtf.py')",
                            required=True)

    args = parser.parse_args()

    ################################################################################

    temp_dir = "./temp_dir"
    os.mkdir(temp_dir)

    compute_overlaps(args,temp_dir)
    list_dfs = []

    for overlap in os.listdir(temp_dir):
        list_dfs.append(build_df_overlaps(os.path.join(os.getcwd(),"temp_dir",overlap),args))

    while len(list_dfs) > 1:
        list_dfs[0] = list_dfs[0].join(list_dfs[1], how="outer")
        del list_dfs[1]


    list_dfs[0].to_csv("final_segmentation_matrix.tsv", sep="\t")


    for overlap in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir,overlap))

    os.rmdir(temp_dir)
    # print(list(list_dfs[0]))
    # print(len(list(list_dfs[0])))
    # print(list_dfs[1].head())
    # print(list_dfs[0].join(list_dfs[1], how="outer"))


if __name__ == "__main__":
    main()
#print(build_df_overlaps("intersect_TSS_with_segmentation.txt", args))
