#!/usr/bin/env python
#
# Usage: collapse_tissue_beds.py METADATA_FILE [hg19|hg38|all]
#

import sys, os

def is_header(line):
     return line.find('File accession\t')==0

def collapse_beds(name, accessions, raw_bed_dir, collapsed_bed_dir):
    bednames = [os.path.join(raw_bed_dir, acc+'.bed.gz') for acc in accessions]
    output_name = os.path.join(collapsed_bed_dir, name.replace(" ","_")+'.bed.gz')
    print " ".join(["zcat",
                    " ".join(bednames),
                    "| bedtools sort -i -",
                    "| bedtools merge -i -",
                    "| gzip -c >",
                    output_name
                   ])


metadatafile = sys.argv[1]
raw_bed_dir = sys.argv[2]
collapsed_bed_dir = sys.argv[3]

hg19_bed_groups = {}
hg38_bed_groups = {}
tissue_ids = {}

with open(metadatafile) as md:

    cols={}
    
    for line in md.xreadlines():
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
            
            bed_groups = hg19_bed_groups if ls[cols['assembly']] =='hg19' else hg38_bed_groups
            
            if not tissue in bed_groups:
               bed_groups[tissue] = {}
            if not lifestage in bed_groups[tissue]:
               bed_groups[tissue][lifestage] = []

            bed_groups[tissue][lifestage].append(ls[cols['accession']])
            
            tissue_ids[ls[cols['tissue']]] = ls[cols['tissue_id']]

#bed_groups = hg19_bed_groups
#for t in bed_groups.keys():
#    counts = ["\n\t"+ls+"("+str(len(bed_groups[t][ls]))+")" for ls in bed_groups[t].keys()]
#    print t, "".join(counts)

#print bed_groups

#
# Analysis of gene regulation in embryonic life stages can be critical to finding causes for developmental disorders.
# Hence, we extracted data for embrionic lifestage (where available) and make it available separately. 
# When data from other lifestages (e.g. adult, child, newborn, unknown) is available for a given tissue/celltype, it is collapsed into "other", e.g. fibroblast of lung (embryonic), fibroblast of lung (other)
# If no embryonic data is available, we make no distinction, e.g. pancreas.
#

EMBRYO = 'embryonic'

#assemblies = ['hg19','hg38'] if sys.argv[2]=='all' else [sys.argv[2]]
assemblies = ['hg19']

for assembly in assemblies:
    
    bed_groups = hg19_bed_groups if assembly =='hg19' else hg38_bed_groups

    for t in bed_groups.keys():
        
        tissue = bed_groups[t]
        accessions=[]
        
        for s in tissue.keys():
            if s is not EMBRYO: 
                accessions.extend(tissue[s])
            
        if EMBRYO in tissue:
            collapse_beds(" ".join([tissue_ids[t],t,"embryonic"]), tissue[EMBRYO], raw_bed_dir, collapsed_bed_dir)
            collapse_beds(" ".join([tissue_ids[t],t,"other"]), accessions, raw_bed_dir, collapsed_bed_dir)
        else:
            collapse_beds(" ".join([tissue_ids[t],t]), accessions, raw_bed_dir, collapsed_bed_dir)
            
     
        
        
