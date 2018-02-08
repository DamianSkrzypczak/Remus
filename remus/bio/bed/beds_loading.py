from pybedtools import BedTool

from remus.bio.time_measurement import time_it


class BedLoader:
    # TODO: docstring
    def __init__(self, src, from_string=False, bed_tool=BedTool):
        # TODO: docstring
        self._bed_tool = bed_tool
        self._src = src
        self._from_string = from_string
        self._bed, self._time_elapsed = self._load_bed()

    @time_it
    def _load_bed(self):
        return self._bed_tool(self._src, from_string=self._from_string)

    @property
    def bed(self):
        # TODO: docstring
        return self._bed

    @property
    def time_elapsed(self, decimal_precision=6):
        # TODO: docstring
        return round(self._time_elapsed, decimal_precision)
