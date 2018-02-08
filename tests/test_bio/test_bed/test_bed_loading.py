import unittest
from unittest import mock

from pybedtools import BedTool

from remus.bio.bed.beds_loading import BedLoader


class TestBedLoader(unittest.TestCase):
    # TODO: docstring
    def setUp(self):
        # TODO: docstring
        self.bed_tool_mock = mock.create_autospec(BedTool)
        self.test_bed_path = "TestPath"

    def test_bed_loading(self):
        # TODO: docstring
        self.test_obj = BedLoader(self.test_bed_path, bed_tool=self.bed_tool_mock)
        beds = self.test_obj.bed
        self.bed_tool_mock.assert_called_once_with(self.test_bed_path, from_string=False)
        self.assertEqual(beds, self.bed_tool_mock(self.test_bed_path, from_string=False))

    def test_loading_time_measurement(self):
        # TODO: docstring
        self.test_obj = BedLoader(self.test_bed_path, bed_tool=self.bed_tool_mock)
        self.assertIsNotNone(self.test_obj.time_elapsed)
