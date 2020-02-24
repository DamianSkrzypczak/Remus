

import unittest

from remus.bio.mirna.registry import MirTarBaseRegistry, MirWalkRegistry, MiRNATargetRegistryFactory

class MiRNATargetRegistryTest(unittest.TestCase):


    def setUp(self):
        self.mir_symbols = ['hsa-let-7a-2-3p', 'hsa-let-7c*', 'hsa-miR-105-5p', \
                            'hsa-miR-1179','hsa-miR-103b', 'hsa-miR-1255b-2-3p', \
                            'hsa-miR-548aq-3p']
                            
        self.mirna_genes = ['MIRLET7A2', 'MIRLET7C', 'MIR105', 'MIR1179', \
                            'MIR103B', 'MIR1255B2', 'MIR548AQ']
                
        
        self.reg = MiRNATargetRegistryFactory.get_instance(MiRNATargetRegistryFactory.MIRTARBASE_KEY)
        self.reg2 = MiRNATargetRegistryFactory.get_instance(MiRNATargetRegistryFactory.MIRWALK_KEY)

    def test_get_mirna_gene_symbols(self):
        mir_genes = self.reg.get_mirna_gene_symbols(self.mir_symbols)
        self.assertEqual(self.mirna_genes, mir_genes)

    def test_get_mirna_gene_symbol(self):
        mir_gene = self.reg.get_mirna_gene_symbols([self.mir_symbols[0]])
        self.assertEqual([self.mirna_genes[0]], mir_gene)

    def test_get_mirna_gene_symbol_empty(self):
        mir_gene = self.reg.get_mirna_gene_symbols([])
        self.assertEqual([], mir_gene)

    def test_HNF1B_gene_query_mirtarbase(self):
        mirnas = self.reg.get_mirnas_targetting_gene('HNF1B')
        self.assertEqual(['hsa-miR-23a-3p'], mirnas)
        
        hnf1b_all = set(['hsa-miR-215-5p', 'hsa-miR-192-5p', 'hsa-miR-23a-3p'])
        mirnas = self.reg.get_mirnas_targetting_gene('HNF1B', include_weak_support=True)
        self.assertEqual(len(hnf1b_all), len(mirnas))
        for m in mirnas:
            self.assertTrue(m in hnf1b_all)

    def test_3gene_query_mirtarbase(self):
        genes = ['HNF1B', 'HNF1A', 'PKD1']
        mirnas = self.reg.get_mirnas_targetting_genes(genes)
        three_genes_mirnas = ['hsa-miR-15b-5p', 'hsa-miR-20a-5p', 'hsa-miR-23a-3p']
        self.assertEqual(len(three_genes_mirnas), len(mirnas))
        for m in mirnas:
            self.assertTrue(m in three_genes_mirnas)
            
    def test_HNF1B_gene_query_mirwalk(self):
        mirnas = self.reg2.get_mirnas_targetting_gene('HNF1B', min_confidence=0.5)
        self.assertEqual(729, len(mirnas))

        mirnas = self.reg2.get_mirnas_targetting_gene('HNF1B', min_confidence=1.0)
        self.assertEqual(205, len(mirnas))
        

    def test_3genes_query_mirwalk(self):
        genes = ['HNF1B', 'HNF1A', 'PKD1']
        mirnas = self.reg2.get_mirnas_targetting_genes(genes, min_confidence=0.5)
        self.assertEqual(1420, len(mirnas))

        mirnas = self.reg2.get_mirnas_targetting_genes(genes, min_confidence=1.0)
        self.assertEqual(486, len(mirnas))

    def tearDown(self):
        pass
