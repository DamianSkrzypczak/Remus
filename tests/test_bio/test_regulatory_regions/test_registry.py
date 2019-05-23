
import re
import unittest

from remus.bio.regulatory_regions.registry import RegulatoryRegionsFilesRegistry

class RegulatoryRegionsFilesRegistryTest(unittest.TestCase):


    def setUp(self):
       self.reg = RegulatoryRegionsFilesRegistry.get_registry('hg19')
       
       self.dummy_sources_map = \
                        {('ONTO1',''): 
                                  {'src1.1' : 'path1.1', 
                                   'src1.2' : 'path1.2', 
                                   'name'   : 'ONTO1name',
                                   'src1.3' : 'path1.3'},
                         ('aONTO2','_embryonic'): 
                                  {'src2.3': 'path2.3', 
                                   'src2.1' : 'path2.1', 
                                   'name'   : 'aONTO2name', 
                                   'src2.2' : 'path2.2'},
                         ('ONTO0','_other'): 
                                  {'src3.3' : 'path3.3', 
                                   'name'   : 'ONTO0name'}
                         }

       self.dummy_tissue_map = {'ONTO1name (src1.1, src1.2, src1.3)': 
                                        {'src1.1': 'path1.1', 
                                         'src1.2': 'path1.2', 
                                         'src1.3': 'path1.3'},
                                 'aONTO2name embryonic (src2.1, src2.2, src2.3)' : 
                                        {'src2.3': 'path2.3', 
                                         'src2.1' : 'path2.1', 
                                         'src2.2' : 'path2.2'},
                                 'ONTO0name other (src3.3)' : 
                                        {'src3.3' : 'path3.3'}
                                }


    ## not really needed any more since the registry is a singleton
    def test_reproducibility_of_tissue_map(self):
        new_reg = RegulatoryRegionsFilesRegistry.get_registry('hg19')        
        self.assertEqual(self.reg.available_tissues, new_reg.available_tissues)
        
    def test_create_available_tissues_map(self):
        tm = self.reg._create_available_tissues_map(self.dummy_sources_map)
        self.assertEqual(sorted(self.dummy_tissue_map.keys()), sorted(tm.keys()))
        for k in self.dummy_tissue_map:
            self.assertEqual(self.dummy_tissue_map[k], tm[k])
            
    def test_get_matching_tissues(self):
        """ TODO to use files in data/ instead of reg._available_tissues """
        limit=100
        text = 'eMbryonic'
        matching = self.reg.get_matching_tissues(text, limit)
    
        pattern = re.compile(text, re.IGNORECASE)
        matching_test = sorted([k for k in self.reg._available_tissues.keys() if pattern.search(k) ])
        
        self.assertEqual(matching, matching_test[:limit])
        
    def tearDown(self):
        pass
