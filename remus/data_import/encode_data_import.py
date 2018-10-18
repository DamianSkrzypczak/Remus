
import os

SUPPORTED_GENOMES = ['hg19', 'GRCh38']

#
# Analysis of gene regulation in some (e.g. embryonic) life stages can be critical to finding causes for developmental disorders.
# Hence, we allow separate extraction of *special* lifestage (where available) and put it in BED file with a lifestage suffix
# Data from all lifestages (e.g. adult, child, newborn, unknown and embryonic) is collapsed into one bed file without suffix, e.g. 
#  - fibroblast_of_lung_embryonic - contains only embryonic peaks
#  - fibroblast_of_lung           - contains merged peaks for all life stages
#

def format_cmd(cmd_echo_tuple, print_echo):
    cmd, echo = cmd_echo_tuple
    return '\n'.join([cmd, echo if print_echo else ""])+"\n\n"

def make_path(directory, name, extension='.bed.gz'):
    return os.path.join(directory, name.replace(" ", "_") + extension)

def get_collapse_beds_script(raw_bed_dir, collapsed_bed_dirs, 
                             bed_groups, tissue_ids, special_life_stages=[], exclude_terms=[], 
                             liftover_chain_files=None, liftoverTo=None, print_echo = True):
    """
    raw_bed_dir:
    collapsed_bed_dirs:
    genome_bed_groups: dictionary of bed_groups. top level is genome build, next tissue, and lifestage
    tissue_ids: dictionary mapping term_ids and term_names for tissues/celltypes
    special_life_stages: life stages that should be collapsed and stored separately
    exclude_terms: tissues/celltypes with the term_id from this list will be skipped
    liftover_chain_files: a dict mapping genome build with liftover chain file to use, if None no liftover is performed
    liftoverTo: dict mapping genome build to lift over
    print_echo: should status messages be included in the script
    return: string with the script content
    """

    script = ""
    
    for genome_build, build_bed_groups in bed_groups.items():
        
        #for build_bed_groups in bed_groups[genome_build].values():
        
            for t, tissue in build_bed_groups.items():
                
                group_name = " ".join([tissue_ids[t], t])
                
                # skip excluded tissues
                if any([t in group_name for t in exclude_terms]):
                    script += "echo Based on exclusion list file for tissue/celltype [%s] was not created.\n\n" % group_name 
                    continue
        
                # aggregate accession ids from all life stages in this tissue
                accessions = []
                for l in tissue.values():
                    accessions.extend(l)
        
                # create file paths
                raw_bed_paths = [make_path(raw_bed_dir, acc) for acc in accessions]
                collapsed_bed_path = make_path(collapsed_bed_dirs[genome_build], group_name)
                
                # aggregate all BEDs in the tissue
                script += format_cmd(get_collapse_beds_command(raw_bed_paths, collapsed_bed_path), print_echo)
                
                # and optionally liftover
                if liftover_chain_files:
                    liftover_bed_path = make_path(collapsed_bed_dirs[liftoverTo[genome_build]], group_name)
                    script += format_cmd(get_liftover_command(liftover_chain_files[genome_build], collapsed_bed_path , liftover_bed_path ), print_echo)
                   
                # process special lifestage groups in similar fashion
                for sls in special_life_stages:
                    if sls in tissue:
                        collapsed_bed_path = make_path(collapsed_bed_dirs[genome_build], " ".join([group_name, sls]))
                        raw_bed_paths = [make_path(raw_bed_dir, acc) for acc in tissue[sls]]
                        script +=  format_cmd(get_collapse_beds_command(raw_bed_paths, collapsed_bed_path), print_echo)

                        if liftover_chain_files:
                            liftover_bed_path = make_path(collapsed_bed_dirs[liftoverTo[genome_build]], " ".join([group_name, sls]))
                            script += format_cmd(get_liftover_command(liftover_chain_files[genome_build], collapsed_bed_path, liftover_bed_path ), print_echo)
    
    
    return script


def get_liftover_command(chain_file, input_bed, output_bed, unmapped_file='/dev/null', liftover_exec='external_resources/liftOver'):
    """
    chain_file: liftover chainfile to use for liftover
    input_bed: input BED file (gzipped)
    output_bed: output BED file (gzipped)
    unmapped_file
    liftover_exec:
    """
    not_compressed_bed_name = output_bed[:-len('.gz')]
    liftover = " ".join([liftover_exec, input_bed, chain_file, not_compressed_bed_name, unmapped_file])
    gzip     = "gzip %s" % not_compressed_bed_name
    cmd = "\n".join([liftover, gzip])
    
    echo = "echo Lifted over \"%s\" to \"%s\" using chain %s" % (input_bed, output_bed, chain_file)
    
    return cmd, echo


def get_collapse_beds_command(raw_beds, collapsed_bed):
    """
    raw_beds: list of paths to raw BEDs to collapse
    collapsed_bed: destination path for collapsed BED
    return: a pair of strings, where first is the shell command collapsing the BEDs
            and the second is echo message stating what file was generated
    """
    
    cmd = " ".join(["zcat",
                    " ".join(raw_beds),
                    "| bedtools sort -i -",
                    "| bedtools merge -i -",
                    "| gzip -c > \"{}\"".format(collapsed_bed)])
    echo = "echo \"Generated %s\"" % collapsed_bed
        
    return cmd, echo



def is_header(line):
    return line.find('File accession\t') == 0


def map_raw_bed_files_to_tissues(metadatafile):
    """
    metadatafile: a tsv matrix of sample and file metadata. Can be downloaded from ENCODE website
    return: a tuple of:
            (1) dictionary mapping term_ids to term_name. both are used in file name later
            (2) a 3 level dictionary of BED files mapped to life_stage, tissue, and genome build, e.g.:
                dict['hg19']['lung']['embryonic'] gives a list of raw BED files for embryonic lungs samples in hg19 coordinates
    """
    
    bed_groups = {}
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
                assembly = ls[cols['assembly']]
                tissue = ls[cols['tissue']]
                lifestage = ls[cols['lifestage']]
                
                if assembly not in bed_groups:
                    bed_groups[assembly] = {}
    
                if tissue not in bed_groups[assembly]:
                    bed_groups[assembly][tissue] = {}
                if lifestage not in bed_groups[assembly][tissue]:
                    bed_groups[assembly][tissue][lifestage] = []
    
                bed_groups[assembly][tissue][lifestage].append(ls[cols['accession']])
    
                tissue_ids[tissue] = ls[cols['tissue_id']]
    
    return tissue_ids, bed_groups
