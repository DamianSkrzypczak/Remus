from collections import OrderedDict

from pybedtools import BedTool

from remus.bio.time_measurement import time_it


class BedsMutualOperation:
    operations = OrderedDict(
        [
            ("intersection", "intersect"),
            ("subtraction", "subtract"),
            ("union", "cat")
        ]
    )

    def __init__(self, beds, operation, **operation_kwargs):
        self._operation_kwargs = operation_kwargs
        if len(beds) == 1:
            self._result_bed = beds[0]
            self._time_elapsed = 0
        else:
            self._beds = beds
            chosen_operation = self._choose_operation(operation)
            self._result_bed, self._time_elapsed = self._execute_with_accumulation(chosen_operation)

    def _choose_operation(self, operation):
        if operation in self.operations:
            return getattr(BedTool, self.operations[operation])
        else:
            raise OperationError('{} has no operation named "{}"'.format(self.__class__, operation))

    @time_it
    def _execute_with_accumulation(self, operation):
        accumulation = self._beds[0]
        for bed in self._beds[1:]:
            accumulation = operation(accumulation, bed, **self._operation_kwargs).sort().merge()
        return accumulation

    @property
    def result(self):

        return self._result_bed

    @property
    def time_elapsed(self, decimal_precision=6):

        return round(self._time_elapsed, decimal_precision)


class BedsFlanker:
    def __init__(self, beds, upstream, downstream, genome):
        self._beds = beds
        self.results_with_times = [self._get_promoter(b, upstream, downstream, genome) for b in self._beds]
        self._times_elapsed = [r[1] for r in self.results_with_times]
        self._result_beds = [r[0] for r in self.results_with_times]

    @time_it
    def _get_promoter(self, bed, upstream, downstream, genome):
        bed_reduced_to_1_size = bed.flank(**{"l": 1, "r": 0, "genome": genome})
        bed_flanked = bed_reduced_to_1_size.flank(**{"l": upstream, "r": downstream, "s": True, "genome": genome})
        bed_sorted = bed_flanked.sort()
        bed_merged_strand_wise_without_distance_between_features_merged = bed_sorted.merge(**{"s": True, "d": 1})
        return bed_merged_strand_wise_without_distance_between_features_merged

    @property
    def results(self):
        return self._result_beds

    def times_elapsed(self, decimal_precision=6):
        return [round(t, decimal_precision) for t in self._times_elapsed]


class MissingBedsException(Exception):
    pass


class OperationError(Exception):
    pass
