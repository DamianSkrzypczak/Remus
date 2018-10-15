import os
import re
from collections import defaultdict

from remus.bio.bed.beds_loading import BedLoader


DATA_DIRECTORIES_MAP = {
    "enhancers/encode": "ENH_EN",
    "enhancers/fantom5": "ENH_F5",
    "chromatin": "CHRM",
    "tss/fantom5": "TSS_F5"
}


class RegulatoryRegionsFilesRegistry:
    def __init__(self, root="data",
                 directories_and_symbols=DATA_DIRECTORIES_MAP,
                 extensions=(".bed", ".bed.gz")):
        sources_map = self.make_sources_map(directories_and_symbols, extensions, root)
        self._available_tissues = self._create_available_tissues_map(sources_map)

    @staticmethod
    def _create_available_tissues_map(sources_map):
        tissues_map = {}
        for onto_num, files_data in sources_map.items():
            source_name = files_data.pop("name")
            name = "{} ({})".format(source_name, ", ".join(files_data.keys()))
            tissues_map[name] = files_data
        return tissues_map

    @staticmethod
    def make_sources_map(directories_and_symbols, extensions, root):
        sources = defaultdict(dict)
        for path in directories_and_symbols:
            for bed in os.listdir(os.path.join(root, path)):
                ontology_and_name = re.match(r"(\w+:\d+)_([^.]+)", bed)
                if ontology_and_name and any([bed.endswith(ext) for ext in extensions]):
                    ontology = ontology_and_name.group(1)
                    name = ontology_and_name.group(2)
                    symbol = directories_and_symbols[path]
                    sources[ontology][symbol] = os.path.join(root, path, bed)
                    sources[ontology]["name"] = name.replace("_expressed_enhancers", "").replace("_promoters", "").replace("_", " ")
        return sources

    @property
    def available_tissues(self):
        return list(self._available_tissues.keys())

    def get_tissue_path(self, name, source_symbol):
        return self._available_tissues[name][source_symbol]

    def get_bed(self, tissue, source_symbol):
        # log.DEBUG('Requested tissue [%s] for source [%s]' % (tissue, source_symbol))
        try:
            tissue_path = self.get_tissue_path(tissue, source_symbol)
        except KeyError:
            # log.DEBUG('No tissue [%s] for source [%s]' % (tissue, source_symbol))
            return None
                
        return BedLoader(tissue_path).bed

    def get_matching_tissues(self, pattern, limit):
        limit = limit if limit else -1
        pattern = pattern if pattern else " "
        try:
            pattern = re.sub("\s+", ".*", pattern, re.IGNORECASE)
            return sorted([i for i in self.available_tissues if re.search(pattern, i, re.IGNORECASE)])[:limit]
        except LookupError:
            return []
