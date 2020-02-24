#!/usr/bin/env python
#
# Usage: collapse_tissue_beds.py METADATA_FILE RAW_BED_DIR COLLAPSED_BED_DIR [hg19|hg38]
#
# Not specifiying genome build produces collapsed date for hg19 and hg38, with liftovers bothways
#

import sys

from encode_data_import import map_raw_bed_files_to_tissues, get_collapse_beds_script, get_collapsed_bed_dir_map
from encode_data_import import SUPPORTED_GENOME_BUILDS, LIFTOVER_BOTH

metadatafile = sys.argv[1]
raw_bed_dir = sys.argv[2]
collapsed_bed_dir = sys.argv[3]

#SPECIAL_LIFE_STAGES = ['embryonic']
SPECIAL_LIFE_STAGES = []
liftover_to = LIFTOVER_BOTH

# map names of raw BED files to tissue term_id, life_stage, and genome build
includes = {"File Status": ["released"],
            "Output type": ["optimal IDR thresholded peaks", "pseudoreplicated IDR thresholded peaks"]}
excludes = {"Biosample term id": ["NTR:0001226"]} #parathyroid adenoma
tissue_ids, bed_groups = map_raw_bed_files_to_tissues(metadatafile, includes, excludes)

if len(sys.argv) > 4:
    genome_build = sys.argv[4]
    if genome_build not in SUPPORTED_GENOME_BUILDS:
        exit('Last (4th) argument must be empty (include both available genome builds) or one of '+ str(SUPPORTED_GENOME_BUILDS.keys()))
    
    # limit to selected genome build and drop liftover
    bed_groups = {genome_build: bed_groups[genome_build]}
    liftover_to=None
 
script = get_collapse_beds_script(raw_bed_dir, collapsed_bed_dir, 
                                  bed_groups, tissue_ids, 
                                  special_life_stages=SPECIAL_LIFE_STAGES,
                                  liftover_to=liftover_to)

print(script)
    


