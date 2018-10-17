
import unittest

from remus.bio.regulatory_regions.registry import RegulatoryRegionsFilesRegistry

class RegulatoryRegionsFilesRegistryTest(unittest.TestCase):


    def setUp(self):
       self.reg = RegulatoryRegionsFilesRegistry()
       
       self.dummy_sources_map = \
                        {'ONTO1': {'src1.1' : 'path1.1', 
                                   'src1.2' : 'path1.2', 
                                   'name'   : 'ONTO1_name',
                                   'src1.3' : 'path1.3'},
                         'aONTO2': {'src2.3': 'path2.3', 
                                   'src2.1' : 'path2.1', 
                                   'name'   : 'aONTO2_name', 
                                   'src2.2' : 'path2.2'},
                         'ONTO0': {'src3.3' : 'path3.3', 
                                   'name'   : 'ONTO0_name'}
                         }

       self.dummy_tissue_map = {'ONTO1_name (src1.1, src1.2, src1.3)': 
                                        {'src1.1': 'path1.1', 
                                         'src1.2': 'path1.2', 
                                         'src1.3': 'path1.3'},
                                 'aONTO2_name (src2.1, src2.2, src2.3)' : 
                                        {'src2.3': 'path2.3', 
                                         'src2.1' : 'path2.1', 
                                         'src2.2' : 'path2.2'},
                                 'ONTO0_name (src3.3)' : 
                                        {'src3.3' : 'path3.3'}
                                }


    def test_reproducibility_of_tissue_map(self):
        new_reg = RegulatoryRegionsFilesRegistry()        
        self.assertEqual(self.reg.available_tissues, new_reg.available_tissues)
        
    def test__create_available_tissues_map(self):
        tm = self.reg._create_available_tissues_map(self.dummy_sources_map)
        self.assertEqual(sorted(self.dummy_tissue_map.keys()), sorted(tm.keys()))
        for k in self.dummy_tissue_map:
            self.assertEqual(self.dummy_tissue_map[k], tm[k])
            
    def tearDown(self):
        pass
