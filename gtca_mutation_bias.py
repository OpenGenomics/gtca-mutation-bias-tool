#!/usr/bin/env python3
"""
Detects samples in a MAF with over-abundance of G>T mutations compared to C>A
mutations.

Usage: gtca_mutation_bias.py <MAF> <output file name> <filter stats file name>
"""

from __future__ import print_function

import sys


if len(sys.argv) != 4:
    sys.exit(__doc__)

maf = open(sys.argv[1], "r")
outputfile = open(sys.argv[2], "w")
filterstats = open(sys.argv[3], "w")

gt_dict = {}
ca_dict = {}
for line in maf:
    if line.startswith("#") or line.startswith("Hugo"):
        continue
    else:
        line = line.split("\t")
        tcga_id = line[15]
        ref = line[10]
        alt = line[12]
        t_depth = line[41]
        t_alt_count = line[39]
        try:
            vaf = float(t_depth) / float(t_alt_count)
        except:
            print(
                "Could not calculate VAF.",
                "Either t_depth or t_alt_count is not an integer.\n",
                "t_depth:\t{}\n".format(t_depth),
                "t_alt_count:\t{}".format(t_alt_count),
                file=sys.stderr
            )
            sys.exit(1)
        if (ref == "G" and alt == "T"):
            if tcga_id in gt_dict:
                gt_dict[tcga_id].append(vaf)
            else:
                gt_dict[tcga_id] = [vaf]
            if tcga_id not in ca_dict:
                ca_dict[tcga_id] = []

        if (ref == "C" and alt == "A"):
            if tcga_id in ca_dict:
                ca_dict[tcga_id].append(vaf)
            else:
                ca_dict[tcga_id] = [vaf]
            if tcga_id not in gt_dict:
                gt_dict[tcga_id] = []

maf.close()

vaf_dict = {}
for k in gt_dict:
    vaf_dict[k] = None
    min_gt_vaf = 0
    gt_count = len(gt_dict[k])
    ca_count = len(ca_dict[k])
    if (gt_count >= ca_count and ca_count > 0):
        min_gt_vaf = min(sorted(gt_dict[k], reverse=True)[0:int(ca_count)])
    vaf_dict[k] = min_gt_vaf

maf = open(sys.argv[1], "r")
filter_dict = {}
ca_dict2 = {}
for line in maf:
    if line.startswith("#") or line.startswith("Hugo"):
        continue
    else:
        # Tumor_Sample_Barcode|Chromosome|Start_Position|End_Position|Reference_Allele|Tumor_Seq_Allele2\tStrandBias
        line = line.split("\t")
        Tumor_Sample_Barcode = line[15]
        Chromosome = line[4]
        Start_Position = line[5]
        End_Position = line[6]
        Reference_Allele = line[10]
        Tumor_Seq_Allele2 = line[12]
        this_key = '|'.join([Tumor_Sample_Barcode,
                             Chromosome,
                             Start_Position,
                             End_Position,
                             Reference_Allele,
                             Tumor_Seq_Allele2])

        if (Reference_Allele == "C" and Tumor_Seq_Allele2 == "A"):
            if Tumor_Sample_Barcode in ca_dict2:
                ca_dict2[Tumor_Sample_Barcode] += 1
            else:
                ca_dict2[Tumor_Sample_Barcode] = 1

        if (Reference_Allele == "G" and Tumor_Seq_Allele2 == "T"):
            if Tumor_Sample_Barcode in filter_dict:
                filter_dict[Tumor_Sample_Barcode][0] += 1
            else:
                # number of variants, number filtered, min_vaf
                filter_dict[Tumor_Sample_Barcode] = [
                    1,
                    0,
                    vaf_dict[Tumor_Sample_Barcode]
                ]

            this_vaf = float(line[41]) / float(line[39])
            if this_vaf < vaf_dict[Tumor_Sample_Barcode]:
                filter_dict[Tumor_Sample_Barcode][1] += 1
                outputfile.write(this_key + "\tStrandBias\n")
            else:
                outputfile.write(this_key + "\n")
        else:
            outputfile.write(this_key + "\n")

maf.close()
outputfile.close()

filterstats.write('\t'.join(["TCGA_ID",
                             "Number_of_G>T_variants",
                             "Number_of_G>T_variants_filtered",
                             "Minimum_VAF_for_inclusion",
                             "Number_of_C>A_variants"]) + "\n")
for k, v in filter_dict.items():
    filterstats.write(
        k +
        "\t" +
        '\t'.join([str(x) for x in v]) +
        "\t" +
        str(ca_dict2[k]) +
        "\n"
    )

filterstats.close()
