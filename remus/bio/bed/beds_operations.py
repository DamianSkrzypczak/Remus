from pybedtools import BedTool

from remus.bio.time_measurement import time_it

class BedOperations:
    
    
    @staticmethod
    #@time_it
    def intersect(beds, merge=False, **kwargs):
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
    #@time_it
    def union(beds, merge=False, **kwargs):
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
    #@time_it
    def get_promoter_region(bed, upstream, downstream, genome):
        """ gets promoter region for each feature in a BED file """
        from pybedtools.featurefuncs import TSS
        return bed.each(TSS, upstream=int(upstream), downstream=int(downstream)).saveas()



class MissingBedsException(Exception):
    pass

