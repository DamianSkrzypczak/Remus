
from remus.bio.time_measurement import time_it

class BedOperationResult:
    
    def __init__(self, result, time):
        self._result = result
        self._time = time
    
    @property
    def result(self):
        return self._result

    @property
    def time_elapsed(self, decimal_precision=3):
        return round(self._time, decimal_precision)



class BedOperations:
    
    
    @staticmethod
    def union(beds, merge=False, **kwargs):
        result, time = BedOperations._union(beds, merge, **kwargs)
        return BedOperationResult(result, time)
    
    @staticmethod
    def intersect(beds, merge=False, **kwargs):
        result, time = BedOperations._intersect(beds, merge, **kwargs)
        return BedOperationResult(result, time)
        
    @staticmethod
    def get_promoter_region(bed, upstream, downstream):
        result, time  = BedOperations._get_promoter_region(bed, upstream, downstream)
        return BedOperationResult(result, time)
        
    
    @staticmethod
    @time_it
    def _intersect(beds, merge=False, **kwargs):
        """ produces a BED with intersection of features from input BEDs. 
        Intervals are sorted and optionally merged"""
        
        if len(beds) == 0: 
            raise MissingBedsException('Empty BED list for intersect operation')
                        
        accumulation = beds[0]
        for bed in beds[1:]:
            accumulation = accumulation.intersect(bed, **kwargs)
            
        accumulation = accumulation.sort()
        if merge:
            accumulation = accumulation.merge()
            
        return accumulation


    @staticmethod
    @time_it
    def _union(beds, merge, **kwargs):
        """ Produces a BED with union of features in all input BEDs. 
        Intervls are sorted and optionally merged """
        
        if len(beds) == 0: 
            raise MissingBedsException('Empty BED list for union operation')
        
        accumulation = beds[0].merge() if merge else beds[0]
        for bed in beds[1:]:
            accumulation = accumulation.cat(bed, postmerge=merge, **kwargs)
            
        accumulation = accumulation.sort()
            
        return accumulation

    @staticmethod
    @time_it
    def _get_promoter_region(bed, upstream, downstream):
        """ gets promoter region for each feature in a BED file """
        from pybedtools.featurefuncs import TSS
        return bed.each(TSS, upstream=int(upstream), downstream=int(downstream)).saveas()



class MissingBedsException(Exception):
    pass

