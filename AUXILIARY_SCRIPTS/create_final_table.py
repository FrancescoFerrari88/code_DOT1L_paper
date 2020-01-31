#!/usr/bin/env python

import pandas as pd
import sys

print("creating final table ...")
out_folder, mark = sys.argv[1], sys.argv[2]

ann_uropa = pd.read_csv("{}/output_uropa/peaks_Annotation_allhits.bed".format(out_folder), sep="\t", header=None)

ann_homer = pd.read_csv("{}/output_homer/{}_Homer_PeaksAnnotation.annot".format(out_folder, mark), sep="\t")
ann_homer = ann_homer[[list(ann_homer)[0],list(ann_homer)[7]]]
ann_homer.columns = ["PeakID","Annotation_homer"]


deac = pd.read_csv("{}/resLFC_apeglm.tsv".format(out_folder), sep="\t")
deac["chr"] = [i.split("_")[0] for i in deac.index]
deac["start"] = [int(i.split("_")[1]) for i in deac.index]
deac["end"] = [int(i.split("_")[2]) for i in deac.index]
deac = deac[["chr","start","end","baseMean","log2FoldChange","lfcSE","pvalue","padj"]]

tr = pd.read_csv("{}/TSS.bed".format(out_folder),sep="\t", index_col=7, header=None)
tr_tr = tr[3].to_dict()
tr.index = tr[3]
tr_genes = tr[6].to_dict()


dixio_symbol=dict()

gr = ann_uropa.groupby(3)
for i in set(ann_uropa[3].values):
#     ann.groupby(3).get_group(i)[15].tolist()
    dixio_symbol[i]=list(set(gr.get_group(i)[17].tolist()))

dixio_ens=dict()

for i in set(ann_uropa[3].values):
#     ann.groupby(3).get_group(i)[15].tolist()
    dixio_ens[i]=list(set(gr.get_group(i)[15].tolist()))


inde = []
gggenes = []
for i,k in dixio_ens.items():
#     print(i,k)
    for j in k:
#         print(i,j)
        inde.append(i)
        gggenes.append(j)
        
final_annot = pd.DataFrame({"peakID":inde,"geneID":gggenes}, index=inde)

final_annot["symbol"] = [tr_genes[i] if i in tr_genes else i for i in final_annot["geneID"]]

final_annot = final_annot.merge(deac, how="right", left_index=True, right_index=True)
final_annot = final_annot.merge(ann_homer, how="right", left_index=True, right_on="PeakID")
final_annot["homer_genes"] = [str(i).split("(")[1].split(")")[0].split(",")[0] if len(str(i).split("("))>=2 else i for i in final_annot['Annotation_homer']]
final_annot["homer_genes"] = [tr_tr[i] if i in tr_tr else i for i in final_annot["homer_genes"]]
final_annot["Annotation_homer"] = [str(i).split("(")[0].strip() for i in final_annot["Annotation_homer"]]
final_annot.sort_values("padj", ascending=True, inplace=True)

final_annot = final_annot[["chr","start","end","peakID","baseMean","log2FoldChange","pvalue","padj","geneID","symbol","Annotation_homer","homer_genes"]]

final_annot.to_csv("{}/REFERENCE_TABLE_{}.bed".format(out_folder, mark), sep="\t", index=False)
#final_annot.to_excel("{}/REFERENCE_TABLE_{}.xls".format(out_folder, mark))

