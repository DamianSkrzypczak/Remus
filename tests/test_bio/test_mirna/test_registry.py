

import unittest

from remus.bio.mirna.registry import MirTarBaseRegistry

class MiRNATargetRegistryTest(unittest.TestCase):


    def setUp(self):
        self.mir_symbols = ['hsa-let-7a-2-3p', 'hsa-let-7c*', 'hsa-miR-105-5p', \
                            'hsa-miR-1179','hsa-miR-103b', 'hsa-miR-1255b-2-3p', \
                            'hsa-miR-548aq-3p']
                            
        self.mirna_genes = ['MIRLET7A2', 'MIRLET7C', 'MIR105', 'MIR1179', \
                            'MIR103B', 'MIR1255B2', 'MIR548AQ']
                            
        self.reg = MirTarBaseRegistry()

    def test_get_mirna_gene_symbols(self):
        mir_genes = self.reg.get_mirna_gene_symbols(self.mir_symbols)
        self.assertEqual(self.mirna_genes, mir_genes)

    def test_get_mirna_gene_symbol(self):
        mir_gene = self.reg.get_mirna_gene_symbols([self.mir_symbols[0]])
        self.assertEqual([self.mirna_genes[0]], mir_gene)

    def test_get_mirna_gene_symbol_empty(self):
        mir_gene = self.reg.get_mirna_gene_symbols([])
        self.assertEqual([], mir_gene)

    def test_HNF1B_gene_query(self):
        mirnas = self.reg.get_mirnas_targetting_gene('HNF1B')
        self.assertEqual(['hsa-miR-23a-3p'], mirnas)
        
        hnf1b_all = set(['hsa-miR-215-5p', 'hsa-miR-192-5p', 'hsa-miR-23a-3p'])
        mirnas = self.reg.get_mirnas_targetting_gene('HNF1B', True)
        self.assertEqual(len(hnf1b_all), len(mirnas))
        for m in mirnas:
            self.assertTrue(m in hnf1b_all)

    def tearDown(self):
        self.reg.teardown_registry()
