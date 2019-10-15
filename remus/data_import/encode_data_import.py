
import os, sys



SUPPORTED_GENOME_BUILDS = {'hg19':'hg19', 'hg38':'GRCh38'}

LIFTOVER_EXECUTABLE = 'external_resources/liftOver'

LIFTOVER_BOTH = {SUPPORTED_GENOME_BUILDS['hg19']: [SUPPORTED_GENOME_BUILDS['hg38']],
                 SUPPORTED_GENOME_BUILDS['hg38']: [SUPPORTED_GENOME_BUILDS['hg19']]}

LIFTOVER_CHAINS = {(SUPPORTED_GENOME_BUILDS['hg19'], SUPPORTED_GENOME_BUILDS['hg38']) : "external_resources/hg19ToHg38.over.chain.gz",
                   (SUPPORTED_GENOME_BUILDS['hg38'], SUPPORTED_GENOME_BUILDS['hg19']) : "external_resources/hg38ToHg19.over.chain.gz"}


def get_collapsed_bed_dir_map(collapsed_bed_dir_prefix, liftover_to):
    """
    liftover_to: dict mapping genome build to lift over genome build list
    collapsed_bed_dir_prefix: path prefix to use
    """
    collapsed_bed_dirs = {b : os.path.join(collapsed_bed_dir_prefix, b) for b in SUPPORTED_GENOME_BUILDS.values()}
    
    if liftover_to:
        for b, lo_builds in liftover_to.items():
            lo_names = [get_liftover_name(b, lo_build) for lo_build in lo_builds]
            for lo_name in lo_names:
                collapsed_bed_dirs[lo_name] = os.path.join(collapsed_bed_dir_prefix, lo_name)

    return(collapsed_bed_dirs)


def get_liftover_name(from_build, to_build):
    return to_build+"_loFrom_"+from_build


def format_cmd(cmd_echo_tuple, print_echo):
    cmd, echo = cmd_echo_tuple
    return '\n'.join([cmd, echo if print_echo else ""])+"\n\n"


def make_path(directory, name, extension='.bed.gz'):
    return os.path.join(directory, name.replace(" ", "_") + extension)


#
# Analysis of gene regulation in some (e.g. embryonic) life stages can be critical to finding causes for developmental disorders.
# Hence, we allow separate extraction of *special* lifestage (where available) and put it in BED file with a lifestage suffix
# Data from all lifestages (e.g. adult, child, newborn, unknown and embryonic) is collapsed into one bed file without suffix, e.g. 
#  - fibroblast_of_lung_embryonic - contains only embryonic peaks
#  - fibroblast_of_lung           - contains merged peaks for all life stages
#
def get_collapse_beds_script(raw_bed_dir, collapsed_bed_dir, 
                             bed_groups, tissue_ids, special_life_stages=[],
                             liftover_to=LIFTOVER_BOTH, liftover_chain_files=LIFTOVER_CHAINS, print_echo = True):
    """
    raw_bed_dir:
    collapsed_bed_dir:
    genome_bed_groups: dictionary of bed_groups. top level is tissue, next lifestage, and genome build
    tissue_ids: dictionary mapping term_ids and term_names for tissues/celltypes
    special_life_stages: life stages that should be collapsed and stored separately
    liftover_to: dict mapping genome build to lift over genome build list
    liftover_chains: a dict mapping pairs of genome builds with liftover chain file to use
    print_echo: should status messages be included in the script
    return: string with the script content
    """

    collapsed_bed_dirs = get_collapsed_bed_dir_map(collapsed_bed_dir, liftover_to)

    separator = "##########################\n\n"
    script = "\n".join(["mkdir -v -p %s" % d for d in collapsed_bed_dirs.values()]) + "\n\n"
    script += "\n".join(["mkdir -v -p %s" % os.path.join(collapsed_bed_dir, build+"_with_liftover") for build in SUPPORTED_GENOME_BUILDS.values()]) + "\n\n"
    script += separator
                
    for t, tissue in bed_groups.items():
            
        group_name = " ".join([tissue_ids[t], t])

        # iterate over genome builds that this tissue is available in
        tissue_gbs = set()
        for ls in tissue:
            tissue_gbs.update(tissue[ls].keys())
        
        for genome_build in tissue_gbs:
    
            # aggregate accession ids from all life stages in this tissue
            accessions = []
            for ls in tissue:
                if genome_build in tissue[ls]:
                    accessions.extend(tissue[ls][genome_build])
    
            # create file paths
            raw_bed_paths = [make_path(raw_bed_dir, acc) for acc in accessions]
            collapsed_bed_path = make_path(collapsed_bed_dirs[genome_build], group_name)
            
            # aggregate all BEDs in the tissue
            script += format_cmd(get_collapse_beds_command(raw_bed_paths, collapsed_bed_path), print_echo)
            
            # and optionally liftover
            if liftover_to and genome_build in liftover_to:
                
                for lo_build in liftover_to[genome_build]:
                    liftover_bed_path = make_path(collapsed_bed_dirs[get_liftover_name(genome_build, lo_build)], group_name, '.liftover.bed.gz')
                    cmd = get_liftover_command(liftover_chain_files[(genome_build, lo_build)], 
                                                collapsed_bed_path, liftover_bed_path, liftover_bed_path+'.unmapped')
                    script += format_cmd(cmd, print_echo)
                
               
            # process special lifestage groups in similar fashion
            for sls in special_life_stages:
                if sls in tissue:
                    collapsed_bed_path = make_path(collapsed_bed_dirs[genome_build], " ".join([group_name, sls]))
                    raw_bed_paths = [make_path(raw_bed_dir, acc) for acc in tissue[sls][genome_build]]
                    script += format_cmd(get_collapse_beds_command(raw_bed_paths, collapsed_bed_path), print_echo)

                    if liftover_to and genome_build in liftover_to:                        
                        for lo_build in liftover_to[genome_build]:
                            liftover_bed_path = make_path(collapsed_bed_dirs[get_liftover_name(genome_build, lo_build)], "_".join([group_name, sls]), '.liftover.bed.gz')
                            cmd = get_liftover_command(liftover_chain_files[(genome_build, lo_build)], 
                                                        collapsed_bed_path, liftover_bed_path)
                            script += format_cmd(cmd, print_echo)
                    

        # merge lifted over data with original genome build data, i.e. 
        # - hg19 and GRCh38 lifted over to hg19
        # - GRCh38 with hg19 lifted over to GRCh38

        if liftover_to:
            script += get_merge_with_liftover_commands(tissue, group_name, special_life_stages, 
                                                        liftover_to, collapsed_bed_dirs, print_echo)

        script += separator
    
    return script


def get_merge_with_liftover_commands(tissue, group_name, special_life_stages, liftover_to, collapsed_bed_dirs, print_echo):

    tissue_gbs = set()
    for ls in tissue:
        tissue_gbs.update(tissue[ls].keys())    

    script = ""

    for lo_from_gb in liftover_to:

        for lo_to_gb in liftover_to[lo_from_gb]:

            if lo_from_gb in tissue_gbs and lo_to_gb not in tissue_gbs:
                # simply copy the lifted-over file, because there is no data for lo_to_gb
                script += f"\necho \"Copying lifted over data for \'{group_name}\'\""
                lo_bed = make_path(collapsed_bed_dirs[get_liftover_name(lo_from_gb, lo_to_gb)], group_name,
                               '.liftover.bed.gz')
                bed_with_lo = make_path(collapsed_bed_dirs[lo_to_gb] + "_with_liftover", group_name)
                script += f"\ncp \"{lo_bed}\" \"{bed_with_lo}\"\n\n"

            elif lo_from_gb in tissue_gbs and lo_to_gb in tissue_gbs:
                # merge the original genome data with liftover data

                script += f"\necho \"Merging original and lifted over data for \'{group_name}\' ({lo_to_gb})\""

                beds_to_merge = [make_path(collapsed_bed_dirs[lo_to_gb], group_name),
                                 make_path(collapsed_bed_dirs[get_liftover_name(lo_from_gb, lo_to_gb)], group_name,
                                           '.liftover.bed.gz')]
                collapsed_with_liftover_bed_path = make_path(collapsed_bed_dirs[lo_to_gb] + "_with_liftover",
                                                             group_name)
                script += format_cmd(get_collapse_beds_command(beds_to_merge, collapsed_with_liftover_bed_path),
                                     print_echo)

                for sls in [s for s in special_life_stages if s in tissue]:
                    beds_to_merge = [make_path(collapsed_bed_dirs[lo_to_gb], " ".join([group_name, sls])),
                                     make_path(collapsed_bed_dirs[get_liftover_name(lo_from_gb, lo_to_gb)],
                                               "_".join([group_name, sls]), '.liftover.bed.gz')]
                    collapsed_with_liftover_bed_path = make_path(collapsed_bed_dirs[lo_to_gb] + "_with_liftover",
                                                                 "_".join([group_name, sls]))
                    script += format_cmd(get_collapse_beds_command(beds_to_merge, collapsed_with_liftover_bed_path),
                                         print_echo)

            else:
                # if lo_from_gb is not among tissue gb, there was no liftover made - nothing to merge
                # if none liftover gbs are among tissue gbs - there is nothing to do
                pass

    return script



def get_liftover_command(chain_file, input_bed, output_bed, unmapped_file='/dev/null', liftover_exec=LIFTOVER_EXECUTABLE):
    """
    chain_file: liftover chainfile to use for liftover
    input_bed: input BED file (gzipped)
    output_bed: output BED file (gzipped)
    unmapped_file
    liftover_exec:
    """
    not_compressed_bed_name = output_bed[:-len('.gz')]
    liftover = "%s \"%s\" %s \"%s\" \"%s\"" % (liftover_exec, input_bed, chain_file, not_compressed_bed_name, unmapped_file)
    gzip     = "bgzip \"%s\"" % not_compressed_bed_name
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
    not_compressed_bed_name = collapsed_bed[:-len('.gz')]
    tmp_bed_name = not_compressed_bed_name+".chunks"

    cmd = ""

    # long lists of input BED files must be processed in chunks
    chunk_size = 20
    for i in range(0, len(raw_beds), chunk_size):
        beds_chunk = raw_beds[i:(i+chunk_size)]

        zcat = " ".join(["zcat",
                         " ".join(["\""+b+"\"" for b in beds_chunk]),
                         "| bedtools sort -i -",
                         "| bedtools merge -i -",
                         " {}> \"{}\"".format("" if i == 0 else ">", tmp_bed_name)])
        cmd += '\n'+zcat

    if len(raw_beds) > chunk_size:  # if the bed list was chunked, sort and merge the tmp_bed
        cmd += "\nbedtools sort -i \"{}\"| bedtools merge -i - > \"{}\"".format(tmp_bed_name, not_compressed_bed_name)
        cmd += "\nrm \"{}\"\n".format(tmp_bed_name)
    else:
        cmd += "\nmv \"{}\" \"{}\"\n".format(tmp_bed_name, not_compressed_bed_name)

    bgzip = "bgzip -f \"%s\"" % not_compressed_bed_name
    tabix = "tabix -p bed \"%s\"" % collapsed_bed
    cmd = '\n'.join([cmd, bgzip, tabix])

    echo = "echo \"Generated %s\"" % collapsed_bed

    return cmd, echo



def is_header(line):
    return line.find('File accession\t') == 0


def passes_includes_and_excludes(ls, cols, include_dict, exclude_dict):
    """
    Checks whether this record (BED) should be included
    ls: list of line elements (line.split())
    cols: map of column names to indices
    include_dict: dictionaty of colname:[values] of accepted column values
    return: True if include and exlude dicts are empty or
            column values specified in the dict match accepted values and don't match unaccepted values
    """
    if include_dict:
        for k in include_dict:
            if ls[cols[k]] not in include_dict[k]:
                sys.stderr.write(f"Dropping dataset because {k} is {ls[cols[k]]}\n")
                return False

    if exclude_dict:
        for k in exclude_dict:
            if ls[cols[k]] in exclude_dict[k]:
                sys.stderr.write(f"Dropping dataset because {k} is {ls[cols[k]]}\n")
                return False

    return True


def map_raw_bed_files_to_tissues(metadatafile, include_dict=None, exclude_dict=None):
    """
    metadatafile: a tsv matrix of sample and file metadata. Can be downloaded from ENCODE website
    include_dict: a dictionary of colname:[value,value,..] specifiing which records should be included
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

                # map column names to indices
                for i,name in enumerate(ls):
                    cols[name] = i

            elif passes_includes_and_excludes(ls, cols, include_dict, exclude_dict):

                assembly = ls[cols['Assembly']]
                tissue = ls[cols['Biosample term name']]
                lifestage = ls[cols['Biosample life stage']] if 'Biosample life stage' in cols else 'all'
                
                if tissue not in bed_groups:
                    bed_groups[tissue] = {}
                if lifestage not in bed_groups[tissue]:
                    bed_groups[tissue][lifestage] = {}
                if assembly not in bed_groups[tissue][lifestage]:
                    bed_groups[tissue][lifestage][assembly] = []
    
                bed_groups[tissue][lifestage][assembly].append(ls[cols['File accession']])
    
                tissue_ids[tissue] = (ls[cols['Biosample term id']]).replace(':', '_')
    
    return tissue_ids, bed_groups
