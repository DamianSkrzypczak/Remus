#!/usr/bin/env python
#
# Usage: collapse_tissue_beds.py METADATA_FILE RAW_BED_DIR COLLAPSED_BED_DIR [hg19|hg38]
#

import os, sys
sys.path.append(os.getcwd())

from remus.data_import.encode_data_import import map_raw_bed_files_to_tissues, get_collapse_beds_script, SUPPORTED_GENOMES

metadatafile = sys.argv[1]
raw_bed_dir = sys.argv[2]
collapsed_bed_dir = sys.argv[3]


genome_build = 'both'
if len(sys.argv) > 4:
    
    genome_build = sys.argv[4]
    
    if genome_build not in SUPPORTED_GENOMES:
        exit('Last (4th) argument must be empty (include both available genome builds) or one of '+ str(SUPPORTED_GENOMES))
    

# map names of raw BED files to tissue term_id, life_stage, and genome build
tissue_ids, hg19_bed_groups, hg38_bed_groups = map_raw_bed_files_to_tissues(metadatafile)

EMBRYO = 'embryonic'
EXCLUDES = ["NTR:0004647"]

genome_build_groups = [hg19_bed_groups, hg38_bed_groups] 
if genome_build=='hg19': genome_build_groups = [hg19_bed_groups]
if genome_build=='hg38': genome_build_groups = [hg38_bed_groups] 
script = get_collapse_beds_script(raw_bed_dir, collapsed_bed_dir, 
                                    genome_build_groups, tissue_ids, 
                                    special_life_stages=[EMBRYO], exclude_terms=EXCLUDES)

print(script)
    


