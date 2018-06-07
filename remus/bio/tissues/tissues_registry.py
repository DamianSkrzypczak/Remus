import os
import re

from remus.bio.bed.beds_loading import BedLoader


class TissuesFilesRegistry:
    def __init__(self, tissues_root=os.path.join("db", "tissues"), extension=".bed"):
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
        pattern = re.sub("\s+", ".*", pattern, re.IGNORECASE)
        return sorted([i for i in self.available_tissues if re.search(pattern, i, re.IGNORECASE)])[:limit]
