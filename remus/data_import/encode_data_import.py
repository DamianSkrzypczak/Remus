
import os

SUPPORTED_GENOMES = ['hg19', 'hg38']

#
# Analysis of gene regulation in some (e.g. embryonic) life stages can be critical to finding causes for developmental disorders.
# Hence, we allow separate extraction of *special* lifestage (where available) and put it in BED file with a lifestage suffix
# Data from all lifestages (e.g. adult, child, newborn, unknown and embryonic) is collapsed into one bed file without suffix, e.g. 
#  - fibroblast_of_lung_embryonic - contains only embryonic peaks
#  - fibroblast_of_lung           - contains merged peaks for all life stages
#


def get_collapse_beds_script(raw_bed_dir, collapsed_bed_dir, bed_groups_list, tissue_ids, special_life_stages=[], exclude_terms=[], print_echo = True):
    """
    bed_groups_list: list of bed_groups
    tissue_ids: dictionary mapping term_ids and term_names for tissues/celltypes
    special_life_stages: life stages that should be collapsed and stored separately
    exclude_terms: tissues/celltypes with the term_id from this list will be skipped
    print_echo: should status messages be included in the script
    return: string with the script content
    """
    script = ""
    
    #for bed_groups in [hg19_bed_groups, hg38_bed_groups]:
    for bed_groups in bed_groups_list:
        
        for t, tissue in bed_groups.items():
            
            group_name = " ".join([tissue_ids[t], t])
            
            # skip exluded tissues
            if any([t in group_name for t in exclude_terms]):
                script += "echo Based on exclusion list file for tissue/celltype [%s] was not created.\n\n" % group_name 
                continue
    
            accessions = []
            for l in tissue.values():
                accessions.extend(l)
    
            cmd, echo = get_collapse_beds_command(group_name, accessions, raw_bed_dir, collapsed_bed_dir)
            script += '\n'.join([cmd, echo if print_echo else ""])
            
            for sls in special_life_stages:
                if sls in tissue:
                    cmd, echo = get_collapse_beds_command(" ".join([group_name, sls]), tissue[sls], raw_bed_dir, collapsed_bed_dir)
                    script += '\n'.join([cmd, echo if print_echo else ""])

    return script

def get_collapse_beds_command(name, accessions, raw_bed_dir, collapsed_bed_dir):
    """
    name: basename of file where collapsed data should be written
    accessions: list of accession ids for raw BED files to collapse
    raw_bed_dir: path to directory with raw BED files
    collapsed_bed_dir: path to directory where collapsed BED should be placed
    return: a pair of strings, where first is the shell command collapsing the BEDs
            and the second is echo message stating what file was generated
    """
    
    bednames = [os.path.join(raw_bed_dir, acc + '.bed.gz') for acc in accessions]
    output_name = os.path.join(collapsed_bed_dir, name.replace(" ", "_") + '.bed.gz')
    cmd = " ".join(["zcat",
                    " ".join(bednames),
                    "| bedtools sort -i -",
                    "| bedtools merge -i -",
                    "| gzip -c > {}".format(output_name)])
    echo = "echo Generated %s\n\n" % output_name
    
    return cmd, echo



def is_header(line):
    return line.find('File accession\t') == 0


def map_raw_bed_files_to_tissues(metadatafile):
    """
    metadatafile: a tsv matrix of sample and file metadata. Can be downloaded from ENCODE website
    return: a tuple of:
            (1) dictionary mapping term_ids to term_name. both are used in file name later
            (2) dictionary of hg19 raw BED groups (by tissue and life stage)
            (3) dictionary of hg38 raw BED groups (by tissue and life stage)
    """
    
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
    
    return tissue_ids, hg19_bed_groups, hg38_bed_groups
