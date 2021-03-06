import os
import re
import logging
from collections import defaultdict

from remus.bio.bed.beds_loading import BedLoader
from remus.bio.bed.beds_operations import BedOperations


def convert_genome_build(genome, hg19_expected="hg19", hg38_expected="GRCh38"):
    if re.match("(hg37|hg19|b37)", genome, re.IGNORECASE):
        return hg19_expected
    elif re.match("(hg38|grch38|b38)", genome, re.IGNORECASE):
        return hg38_expected
    raise InvalidGenomeBuildException(genome)


def get_dirname_for_merged_with_liftover(genome_build):
    return genome_build + "_with_liftover";


class RegulatoryRegionsFilesRegistry:
    
    instances = None      # dictionary of singleton objects
    
    FANTOM5_PROMOTERS_KEY = "PR_F"
    FANTOM5_ENHANCERS_KEY = "ENH_F"

    SCREEN_PROMOTERS_KEY = "PR_S"
    SCREEN_ENHANCERS_KEY = "ENH_S"
    SCREEN_CHROMATIN_KEY = "CHR_S"

    ENCODE_ENHANCERS_KEY = "ENH_E"
    ENCODE_CHROMATIN_KEY = "CHR_E"

    DATA_DIRECTORIES_MAP = {
        "enhancers/encode": ENCODE_ENHANCERS_KEY,
        "enhancers/fantom5": FANTOM5_ENHANCERS_KEY,
        "enhancers/screen": SCREEN_ENHANCERS_KEY,
        "chromatin/encode": ENCODE_CHROMATIN_KEY,
        "chromatin/screen": SCREEN_CHROMATIN_KEY,
        "promoters/fantom5": FANTOM5_PROMOTERS_KEY,
        "promoters/screen": SCREEN_PROMOTERS_KEY
    }
    
    @staticmethod
    def get_registry(genome_build='hg19'):

        genome_build = convert_genome_build(genome_build)

        if RegulatoryRegionsFilesRegistry.instances is None:
            RegulatoryRegionsFilesRegistry.instances = {}
        if genome_build not in RegulatoryRegionsFilesRegistry.instances:
            RegulatoryRegionsFilesRegistry.instances[genome_build] = RegulatoryRegionsFilesRegistry(genome_build)
        return RegulatoryRegionsFilesRegistry.instances[genome_build]

    def __init__(self, genome_build, merge_lifted_over=True, root="data",
                 directories_and_symbols=DATA_DIRECTORIES_MAP,
                 extensions=(".bed", ".bed.gz")):
        self.logger = logging.getLogger(self.__class__.__name__)
        sources_map = self._make_sources_map(directories_and_symbols, genome_build, merge_lifted_over, extensions, root)
        self._available_tissues = self._create_available_tissues_map(sources_map)

    def _make_sources_map(self, directories_and_symbols, genome_build, merge_lifted_over, extensions, root):
        
        merge_with_liftover_genome_build = get_dirname_for_merged_with_liftover(genome_build)
        
        self.logger.info("Making sources map for root: %s ; paths: %s ; genome: %s (%smerging with lifted over coordinates); and extensions: %s" \
                         % (root, str(directories_and_symbols),
                            genome_build, "" if merge_lifted_over else "not ",
                            str(extensions)))
        
        pattern = re.compile(r"(\w+_\d+)_(.+?)(|_embryonic)\.")
        
        sources = defaultdict(dict)
        for path in directories_and_symbols:
                       
            genome_build_dir = genome_build
            if merge_lifted_over and not os.path.isdir(os.path.join(root, path, merge_with_liftover_genome_build)):
                self.logger.warning("BED dir for path [%s] and genome [%s] (merged with lifted over) does not exist. Using defult genome build dir: [%s]" \
                                    % (path, merge_with_liftover_genome_build, genome_build))
            else:
                genome_build_dir = merge_with_liftover_genome_build
        
                
            if not os.path.isdir(os.path.join(root, path, genome_build_dir)):
                self.logger.warning("BED dir for path [%s] and genome [%s] does not exist. Skipping." % (path, genome_build))
                continue 
                
            for bed in os.listdir(os.path.join(root, path, genome_build_dir)):
                
                termid_and_name = pattern.match(bed)
                
                if not termid_and_name and any([bed.endswith(ext) for ext in extensions]):
                    self.logger.warn("File %s not matching expected file name format" % bed)
                    continue
                    
                if any([bed.endswith(ext) for ext in extensions]):
                    termid = termid_and_name.group(1)
                    name = termid_and_name.group(2)
                    life_stage = termid_and_name.group(3)
                    self.logger.debug("%s was split into [%s],[%s],[%s]" % (bed, termid, name, life_stage))
                    
                    symbol = directories_and_symbols[path]
                    
                    sources[termid, life_stage][symbol] = os.path.join(root, path, genome_build_dir, bed)
                    sources[termid, life_stage]["name"] = name.replace("_", " ")
        
        self.logger.debug("Sources map:\n%s" % str(sources))
        return sources

    def _create_available_tissues_map(self, sources_map):
        tissues_map = {}
        for (_, life_stage), files_data in sources_map.items():
            source_name = files_data.pop("name")
            name = "{}{} ({})".format(source_name, life_stage.replace("_"," "), ", ".join(sorted(files_data.keys())))
            tissues_map[name] = files_data
        
        self.logger.debug("Tissues map:\n%s" % str(tissues_map))
        return tissues_map

    @property
    def available_tissues(self):
        self.logger.debug("Available tissues:\n%s" % str(self._available_tissues.keys()))
        return list(self._available_tissues.keys())

    # not used any more. Substituted by get_bed_fragment
    def get_bed(self, tissue, source_symbol):
        return self.get_bed_fragment(tissue, source_symbol, None)

    def get_bed_fragment(self, tissue, source_symbol, regions):
        """ 
        Get slice of a BED. Filtering on non-tabixed BED files is not supported.
        If regions is None, entire BED is returned
        """
        
        self.logger.info('Requested {}tissue [{}] from source [{}]'
                        .format("fragment [%s] from " % regions if regions else "", 
                        tissue, 
                        source_symbol))

        if tissue not in self._available_tissues:
            raise InvalidTissueNameException("Querried tissue [%s] was not among available tissue keys:\n%s" % 
                    (tissue, str(self._available_tissues.keys())))
       
        try:
            bed_path = self._available_tissues[tissue][source_symbol]
            track_name = source_symbol + "(" + tissue.split('(')[0].strip().replace(" ", "_") + ")"
            self.logger.info('Found %s. Adding name %s' % (bed_path, track_name))

            full_bed = BedLoader(bed_path)

            if regions is None:
                return BedOperations.add_name(full_bed.bed, track_name)
            else:
                beds = [ full_bed.filter_by(i) for i in regions ]
                if any(beds):
                    filtered_bed = BedOperations.union([e for e in beds if e], merge=False).result
                    return BedOperations.add_name(filtered_bed, track_name)

        except KeyError:
            self.logger.info('No tissue [%s] in source [%s]' % (tissue, source_symbol))

        return None

    def get_matching_tissues(self, pattern, limit):
        pattern = pattern if pattern else " "
        try:
            pattern = re.sub("\s+", ".*", pattern, re.IGNORECASE)
            matches = sorted([i for i in self.available_tissues if re.search(pattern, i, re.IGNORECASE)])
            return matches[:limit] if limit else matches
        except re.error:
            return []
            


class InvalidTissueNameException(Exception):
    pass


class InvalidGenomeBuildException(Exception):
    pass

