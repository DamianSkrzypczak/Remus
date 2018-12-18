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


class RegulatoryRegionsFilesRegistry:
    
    instances = None      # dictionary of singleton objects
    
    FANTOM5_TSS_KEY = "TSS_F5"
    FANTOM5_ENHANCERS_KEY = "ENH_F5"
    ENCODE_ENHANCERS_KEY = "ENH_ENC"
    ENCODE_CHROMATIN_KEY = "CHRM"
    
    DATA_DIRECTORIES_MAP = {
        "enhancers/encode": ENCODE_ENHANCERS_KEY,
        "enhancers/fantom5": FANTOM5_ENHANCERS_KEY,
        "chromatin": ENCODE_CHROMATIN_KEY,
        "tss/fantom5": FANTOM5_TSS_KEY
    }
    
    @staticmethod
    def get_registry(genome_build='hg19'):

        genome_build = convert_genome_build(genome_build)

        if RegulatoryRegionsFilesRegistry.instances == None:
            RegulatoryRegionsFilesRegistry.instances = {}
        if genome_build not in RegulatoryRegionsFilesRegistry.instances:
            RegulatoryRegionsFilesRegistry.instances[genome_build] = RegulatoryRegionsFilesRegistry(genome_build)
        return RegulatoryRegionsFilesRegistry.instances[genome_build]
    
    
    def __init__(self, genome_build, root="data",
                 directories_and_symbols=DATA_DIRECTORIES_MAP,
                 extensions=(".bed", ".bed.gz")):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        sources_map = self._make_sources_map(directories_and_symbols, genome_build, extensions, root)
        
        self._available_tissues = self._create_available_tissues_map(sources_map)
        

    def _make_sources_map(self, directories_and_symbols, genome_build, extensions, root):
        
        self.logger.info("Making sources map for root: %s ; paths: %s ; genome: %s ; and extensions: %s" \
                            % (root, str(directories_and_symbols), genome_build, str(extensions)))
        
        pattern = re.compile(r"(\w+:\d+)_(.+?)(|_embryonic)\.")
        
        sources = defaultdict(dict)
        for path in directories_and_symbols:
            
            if not os.path.isdir(os.path.join(root, path, genome_build)):
                self.logger.warn("BED dir for path [%s] and genome [%s] does not exist. Skipping." % (path, genome_build))
                continue 
                
            for bed in os.listdir(os.path.join(root, path, genome_build)):
                
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
                    
                    sources[termid, life_stage][symbol] = os.path.join(root, path, genome_build, bed)
                    sources[termid, life_stage]["name"] = name.replace("_expressed_enhancers", "").replace("_promoters", "").replace("_", " ")
        
        self.logger.debug("Sources map:\n%s" % str(sources))
        return sources

    def _create_available_tissues_map(self, sources_map):
        tissues_map = {}
        for (_, life_stage), files_data in sources_map.items():
            source_name = files_data.pop("name")
            name = "{}{} ({})".format(source_name, life_stage, ", ".join(sorted(files_data.keys())))
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
            self.logger.info('Found %s' % bed_path)
            full_bed = BedLoader(bed_path)
            
            if regions is None:
                return full_bed.bed
            else:
                beds = [ full_bed.filter_by(i) for i in regions ]
                if any(beds):
                    return BedOperations.union([e for e in beds if e], merge=False).result
                    
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

