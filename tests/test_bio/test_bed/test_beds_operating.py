import unittest
from unittest import mock

from ddt import ddt, data

from remus.bio.bed.beds_operations import BedsMutualOperation, OperationError


@ddt
class TestBedsOperator(unittest.TestCase):
    operations = BedsMutualOperation.operations

    def setUp(self):
        self.bed1 = mock.Mock()
        self.bed2 = mock.Mock()

    @data(*BedsMutualOperation.operations)
    def test_tissue_operations_dispatching(self, operation):
        with mock.patch("remus.bio.bed.beds_operations.BedTool") as bed_tool_mock:
            test_obj = BedsMutualOperation(beds=[self.bed1, self.bed2], operation=operation)
            result = test_obj.result
            expected_operation = self._get_operation(operation)
            operation = getattr(bed_tool_mock, expected_operation)
            expected_mock_objects = [self.bed1, self.bed2]
            operation.assert_called_once_with(*expected_mock_objects)
            self.assertEqual(operation(*expected_mock_objects).sort(), result)

    def _get_operation(self, operation):
        expected_operation = self.operations[operation]
        return expected_operation

    @data(*BedsMutualOperation.operations)
    def test_operation_time_measurement(self, operation):
        with mock.patch("remus.bio.bed.beds_operations.BedTool"):
            test_obj = BedsMutualOperation(beds=[self.bed1, self.bed2], operation=operation)
            test_obj.result(operation)
            self.assertIsNotNone(test_obj.time_elapsed)

    def test_unavailable_operation_error(self):
        with self.assertRaises(OperationError):
            test_obj = BedsMutualOperation(beds=[self.bed1, self.bed2], operation="Not Valid Operation")
