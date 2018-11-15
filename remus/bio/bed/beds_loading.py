import os

from pybedtools import BedTool

from remus.bio.time_measurement import time_it


class BedLoader:

    def __init__(self, src, from_string=False, bed_tool=BedTool):
        self._bed_tool = bed_tool
        self._src = src
        self._from_string = from_string
        self._tabixed = not from_string and src[-3:]==".gz" and os.path.exists(src+'.tbi')
        self._bed, self._time_elapsed = self._load_bed()

    @time_it
    def _load_bed(self):
        return self._bed_tool(self._src, from_string=self._from_string)

    @property
    def bed(self):
        return self._bed

    @property
    def time_elapsed(self, decimal_precision=6):
        return round(self._time_elapsed, decimal_precision)

    @property
    def tabixed(self):
        return self._tabixed

    def filter_by(self, interval):
        if not self._tabixed:
            raise NotImplementedError("Attemted filtering BED without index")
            
        try:
            return self._bed.tabix_intervals(interval)
        except ValueError:
            return None
    
