#!/usr/bin/env python
#
# Usage: collapse_tissue_beds.py METADATA_FILE [hg19|hg38|all]
#


import os
import sys


def is_header(line):
    return line.find('File accession\t') == 0


def collapse_beds(name, accessions, assembly):
    bednames = [os.path.join(assembly, 'raw', acc + '.bed.gz') for acc in accessions]
    output_name = os.path.join(assembly, name.replace(" ", "_") + '.bed.gz')
    if not "NTR:0004647" in output_name:
        print(" ".join(["zcat",
                        " ".join(bednames),
                        "| bedtools sort -i -",
                        "| bedtools merge -i -",
                        "| gzip -c >",
                        "\"{}\"".format(output_name)
                        ]))
        print(r'printf "File {} generated\n"'.format(output_name))


metadatafile = sys.argv[1]

hg19_bed_groups = {}
hg38_bed_groups = {}
tissue_ids = {}

with open(metadatafile) as md:
    cols = {}

    for line in md.readlines():
        ls = line.split('\t')
        if is_header(line):
            cols['accession'] = ls.index('File accession')
            cols['tissue'] = ls.index('Biosample term name')
            cols['tissue_id'] = ls.index('Biosample term id')
            cols['lifestage'] = ls.index('Biosample life stage')
            cols['assembly'] = ls.index('Assembly')
        else:
            tissue = ls[cols['tissue']]
            lifestage = ls[cols['lifestage']]

            bed_groups = hg19_bed_groups if ls[cols['assembly']] == 'hg19' else hg38_bed_groups

            if not tissue in bed_groups:
                bed_groups[tissue] = {}
            if not lifestage in bed_groups[tissue]:
                bed_groups[tissue][lifestage] = []

            bed_groups[tissue][lifestage].append(ls[cols['accession']])

            tissue_ids[tissue] = ls[cols['tissue_id']]

# bed_groups = hg19_bed_groups
# for t in bed_groups.keys():
#    counts = ["\n\t"+ls+"("+str(len(bed_groups[t][ls]))+")" for ls in bed_groups[t].keys()]
#    print t, "".join(counts)

# print bed_groups

#
# Analysis of gene regulation in embryonic life stages can be critical to finding causes for developmental disorders.
# Hence, we extracted data for embryonic lifestage (where available) and make it available separately. 
# Data from all lifestages (e.g. adult, child, newborn, unknown and embryonic) is collapsed into one bed file without suffix, e.g. 
#  - fibroblast of lung (embryonic) - contains only embryonic peaks
#  - fibroblast of lung             - contains merged peaks for all life stages
#

EMBRYO = 'embryonic'

assemblies = ['hg19', 'hg38'] if sys.argv[2] == 'all' else [sys.argv[2]]

for assembly in assemblies:

    bed_groups = hg19_bed_groups if assembly == 'hg19' else hg38_bed_groups

    for t, tissue in bed_groups.items():

        accessions = []
        for l in tissue.values():
            accessions.extend(l)

        collapse_beds(" ".join([tissue_ids[t], t]), accessions, raw_bed_dir, collapsed_bed_dir)
        
        if EMBRYO in tissue:
            collapse_beds(" ".join([tissue_ids[t], t, "embryonic"]), tissue[EMBRYO], raw_bed_dir, collapsed_bed_dir)


