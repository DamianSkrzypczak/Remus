from collections import OrderedDict

from pybedtools import BedTool

from remus.bio.time_measurement import time_it


class BedsOperation:
    # TODO: docstring

    operations = OrderedDict(
        [
            ("intersection", "intersect"),
            ("subtraction", "subtract"),
            ("union", "cat")
        ]
    )

    def __init__(self, beds, operation):
        # TODO: docstring
        if len(beds) < 2:
            raise MissingBedsException("{} requires a collection of at least two Beds".format(self.__class__))
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
            accumulation = operation(accumulation, bed).sort().merge()
        return accumulation

    @property
    def result(self):
        # TODO: docstring
        return self._result_bed

    @property
    def time_elapsed(self, decimal_precision=6):
        # TODO: docstring
        return round(self._time_elapsed, decimal_precision)


class MissingBedsException(Exception):
    # TODO: docstring
    pass


class OperationError(Exception):
    # TODO: docstring
    pass
