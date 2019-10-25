#!/usr/bin/env python
#
# Usage: split_screen_beds.py METADATA_FILE RAW_BED_DIR SPLIT_BED_DIR
#
#

import os, sys, glob
import gzip

from encode_data_import import passes_includes_and_excludes, get_collapse_beds_script, get_collapsed_bed_dir_map
from encode_data_import import SUPPORTED_GENOME_BUILDS


def is_header(line):
    return line.find('ID\tAccession\t') == 0


def choose_existing(id1, id2, dir):
    p1 = os.path.join(dir, id1+".bed.gz")
    p2 = os.path.join(dir, id2+".bed.gz")
    e1, e2 = os.path.exists(p1), os.path.exists(p2)
    if e1 and not e2:
        return id1
    elif not e1 and e2:
        return id2
    else:
        raise Exception(f"One and only one of files [{p1}, {p2}] should exist.")


def map_raw_bed_files_to_tissues(metadatafile, raw_bed_dir, include_dict=None, exclude_dict=None):
    """
    Note that the format for ENCODE SCREEN metadata is different that for ENCODE Data Matrix metadata.
    metadatafile: a tsv matrix of sample and file metadata. Can be downloaded from ENCODE website.
    include_dict: a dictionary of colname:[value,value,..] specifing which records should be included
    return: a tuple of:
            (1) dictionary mapping term_ids to term_name. both are used in file name later
            (2) a 3 level dictionary of BED files mapped to life_stage, tissue, and genome build, e.g.:
                dict['hg19']['lung']['embryonic'] gives a list of raw BED files for embryonic lungs samples in hg19 coordinates
    """

    bed_groups = {}
    tissue_ids = {}

    with open(metadatafile) as md:

        past_header = False
        cols = {}

        for line in md.readlines():
            ls = line.split('\t')
            if not past_header and is_header(line):
                past_header = True
                # map column names to indices
                for i, name in enumerate(ls):
                    cols[name] = i

            elif past_header and \
                 passes_includes_and_excludes(ls, cols, include_dict, exclude_dict) and \
                 ls[cols["Description"]].startswith("5-group"):

                assembly = ls[cols['Genome assembly']]
                tissue = ls[cols['Biosample term name']]
                lifestage = ls[cols['Life stage']]

                if tissue not in bed_groups:
                    bed_groups[tissue] = {}
                if lifestage not in bed_groups[tissue]:
                    bed_groups[tissue][lifestage] = {}
                if assembly not in bed_groups[tissue][lifestage]:
                    bed_groups[tissue][lifestage][assembly] = []

                # order of bed and bigBed files is not static, and we only download BED files,
                # so the right accession is chosen based on existence of the downloaded BED
                bed_accession1 = ls[cols['Files']].split("/")[2]
                bed_accession2 = ls[cols['Files']].split("/")[5]
                bed_accession = choose_existing(bed_accession1, bed_accession2, raw_bed_dir)

                bed_groups[tissue][lifestage][assembly].append(bed_accession)

                biosample_type = ls[cols['Biosample ontology']]
                ontology_term  = '_'.join((biosample_type.split("/")[2]).split("_")[-2:])
                tissue_ids[tissue] = ontology_term

    return tissue_ids, bed_groups


def split_beds(raw_bed_dir, skip_splitting=False,
               groups=dict(dnase_only="6,218,147", promoter="255,0,0",
                           enhancer="255,205,0", ctcf="0,176,240",
                           dnase_only_dnase_NA="140,140,140", inactive="225,225,225")):
    """ Split records in the BED file according to grouping in the last column """

    beds = glob.glob(os.path.join(raw_bed_dir,"*.bed.gz"))
    bed_types = ['enhancers','promoters','insulators','chromatin']
    split_bed_dirs = [os.path.join(raw_bed_dir, type) for type in bed_types]
    for dir in split_bed_dirs:
        if not os.path.exists(dir) and not os.path.isdir(dir):
            os.mkdir(dir)

    if not skip_splitting:

        for bed in beds:
            id = os.path.basename(bed)
            split_beds = [os.path.join(dir, id) for dir in split_bed_dirs]

            with gzip.open(bed, "rt") as b, \
                    gzip.open(split_beds[0], "wt") as enhb, \
                    gzip.open(split_beds[1], "wt") as promb, \
                    gzip.open(split_beds[2], "wt") as insb, \
                    gzip.open(split_beds[3], "wt") as chromb:

                for l in b:
                    group = l.strip().split("\t")[8]
                    if group == groups['enhancer']:
                        enhb.write(l)
                    elif group == groups['promoter']:
                        promb.write(l)
                    elif group == groups['ctcf']:
                        insb.write(l)
                    elif group == groups['dnase_only']:
                        chromb.write(l)
                    elif group not in groups.values():
                        raise Exception(f"Unknown group {group} in the input BED file: {bed}")

    return split_bed_dirs


def get_remove_empty_beds(dirs):
    paths = ' '.join([f"{dir}/*/*.bed.gz" for dir in dirs])
    return "\n\necho Removing empty BED files\n" + \
            f"for f in {paths}; do \n" + \
            "  if [ `zcat $f | wc -l` == \"0\" ]; then \n" + \
            "      rm $f ${f}.tbi \n" + \
            "  fi \n" + \
            "done \n"


def get_remove_empty_dirs(dirs):
    echo = "\necho Removing empty dirs\n"
    return echo + '\n'.join([f"rmdir {dir}/GRCh38 {dir}/hg19_with_liftover" for dir in dirs])


def move_beds_to_data(bed_dirs, raw_bed_dir):
    names = [os.path.basename(d) for d in bed_dirs]
    root = os.path.dirname(raw_bed_dir)
    script = "\n\necho Moving collapsed BEDs to data directory\n"
    for i in range(0, len(bed_dirs)):
        path = bed_dirs[i]
        new_path = os.path.join(root, names[i], "screen")
        link_path = os.path.join("..", names[i], "screen")
        script += f"mv {path} {new_path}; ln -s {link_path} {path}\n"

    return script


if __name__ == "__main__":

    metadatafile = sys.argv[1]
    raw_bed_dir = sys.argv[2]
    output_bed_dir = sys.argv[3]

    split_raw_bed_dirs = split_beds(raw_bed_dir)

    # map names of raw BED files to tissue term_id, life_stage, and genome build
    excludes = {"Biosample term name": []}
    tissue_ids, bed_groups = map_raw_bed_files_to_tissues(metadatafile, raw_bed_dir, exclude_dict=excludes)

    script = ""
    collapsed_bed_dirs = [os.path.join(output_bed_dir, os.path.basename(dir)) for dir in split_raw_bed_dirs]
    for i in range(0,4):
        script += get_collapse_beds_script(split_raw_bed_dirs[i], collapsed_bed_dirs[i],
                                      bed_groups, tissue_ids,
                                      special_life_stages=['embryonic'],
                                      liftover_to={SUPPORTED_GENOME_BUILDS['hg19']: [SUPPORTED_GENOME_BUILDS['hg38']]})

    script += get_remove_empty_beds(collapsed_bed_dirs)
    script += get_remove_empty_dirs(collapsed_bed_dirs)
    script += move_beds_to_data(collapsed_bed_dirs, output_bed_dir)
    print(script)



