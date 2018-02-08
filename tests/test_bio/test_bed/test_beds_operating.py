import unittest
from unittest import mock

from ddt import ddt, data

from remus.bio.bed.beds_operations import BedsOperation, MissingBedsException, OperationError


@ddt
class TestBedsOperator(unittest.TestCase):
    # TODO: docstring
    operations = BedsOperation.operations

    def setUp(self):
        # TODO: docstring
        self.bed1 = mock.Mock()
        self.bed2 = mock.Mock()

    @data(*BedsOperation.operations)
    def test_tissue_operations_dispatching(self, operation):
        # TODO: docstring
        with mock.patch("remus.bio.bed.beds_operations.BedTool") as bed_tool_mock:
            test_obj = BedsOperation(beds=[self.bed1, self.bed2], operation=operation)
            result = test_obj.result
            expected_operation = self._get_operation(operation)
            operation = getattr(bed_tool_mock, expected_operation)
            expected_mock_objects = [self.bed1, self.bed2]
            operation.assert_called_once_with(*expected_mock_objects)
            self.assertEqual(result, operation(*expected_mock_objects).sort().merge())

    def _get_operation(self, operation):
        # TODO: docstring
        expected_operation = self.operations[operation]
        return expected_operation

    @data(*BedsOperation.operations)
    def test_operation_time_measurement(self, operation):
        # TODO: docstring
        with mock.patch("remus.bio.bed.beds_operations.BedTool"):
            test_obj = BedsOperation(beds=[self.bed1, self.bed2], operation=operation)
            test_obj.result(operation)
            self.assertIsNotNone(test_obj.time_elapsed)

    def test_unavailable_operation_error(self):
        # TODO: docstring
        with self.assertRaises(OperationError):
            test_obj = BedsOperation(beds=[self.bed1, self.bed2], operation="Not Valid Operation")
            test_obj.result("Dummy Operation that not exists")

    def test_wrong_number_of_beds(self):
        # TODO: docstring
        with self.assertRaises(MissingBedsException):
            BedsOperation(beds=[self.bed1], operation=list(self.operations.keys())[0])
