import os
import pprint
import re
from collections import defaultdict

from remus.bio.bed.beds_loading import BedLoader


class RegulatoryRegionsFilesRegistry2:
    def __init__(self, tissues_root=os.path.join("data", "enhancers"), extension=".bed"):
        self._root = tissues_root
        self._tissues_sub_dirs = next(os.walk(self._root))[1]
        self._available_tissues = {}
        for subdir in self._tissues_sub_dirs:
            for t in os.listdir(os.path.join(self._root, subdir)):
                if t.endswith(extension):
                    t_name = "{1} ({0})".format(subdir, self.strip_filename(t))
                    self._available_tissues[t_name] = os.path.join(self._root, subdir, t)
        self._all_tissues_string = " ".join(self.available_tissues)

    @staticmethod
    def strip_filename(filename):
        filename = filename.replace("_enhancers", "")
        filename = filename.replace("_expressed", "")
        return re.match("(\w+[\d:_]+)(.*)(\.bed)", filename).groups()[1].replace("_", " ")

    @property
    def available_tissues(self):
        return list(self._available_tissues.keys())

    def get_tissue_path(self, name):
        return self._available_tissues[name]

    def get_bed(self, tissue):
        tissue_path = self.get_tissue_path(tissue)
        return BedLoader(tissue_path).bed

    def split_tissues_by_sub_dirs(self, tissues_names):
        results = {}
        for sub_dir in self._tissues_sub_dirs:
            results[sub_dir] = [tn for tn in tissues_names if tn.endswith(sub_dir)]

    def get_matching_tissues(self, pattern, limit):
        limit = limit if limit else -1
        pattern = pattern if pattern else " "
        try:
            pattern = re.sub("\s+", ".*", pattern, re.IGNORECASE)
            return sorted([i for i in self.available_tissues if re.search(pattern, i, re.IGNORECASE)])[:limit]
        except LookupError:
            return []


DATA_DIRECTORIES_MAP = {
    "enhancers/encode": "ENH_EN",
    "enhancers/fantom5": "ENH_F5",
    "chromatin": "CHRM",
    "transcription_start_sites": "TSS"
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
                    sources[ontology]["name"] = name.replace("_expressed_enhancers", "").replace("_", " ")
        return sources

    @property
    def available_tissues(self):
        return list(self._available_tissues.keys())

    def get_tissue_path(self, name, source_symbol):
        return self._available_tissues[name][source_symbol]

    def get_bed(self, tissue, source_symbol):
        try:
            tissue_path = self.get_tissue_path(tissue, source_symbol)
        except KeyError:
            return None
        print(tissue_path)
        return BedLoader(tissue_path).bed

    def get_matching_tissues(self, pattern, limit):
        limit = limit if limit else -1
        pattern = pattern if pattern else " "
        try:
            pattern = re.sub("\s+", ".*", pattern, re.IGNORECASE)
            return sorted([i for i in self.available_tissues if re.search(pattern, i, re.IGNORECASE)])[:limit]
        except LookupError:
            return []


if __name__ == '__main__':
    x = RegulatoryRegionsFilesRegistry2()
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(x.available_tissues)
