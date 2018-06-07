import os

from remus.bio.bed.beds_loading import BedLoader


class TranscriptionStartSitesRegistry:

    def __init__(self, root=os.path.join("db", "transcription_start_sites")):
        self._root = root
        self.filename = "promoter_data.bed"

    def get_bed(self):
        return BedLoader(os.path.join(self._root, self.filename)).bed
