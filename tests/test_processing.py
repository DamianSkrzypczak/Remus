

import unittest
import warnings
from unittest.mock import patch, patch

from pybedtools import BedTool
from remus.processing import BedsProcessor
from remus.bio.genes.registry import GenesDBRegistry

#
# TODO
#
# Would be awesome to separate real genes and their coordinates from the tests
# and use only dummy gene coordinates. Could probably be done by mocking GemesDBRegistry.
# 
class TestRemusProcessing(unittest.TestCase):


    def setUp(self):
        
        self.genes_reg = GenesDBRegistry()
        
        self.genome = 'hg19'
        
        #gene1coord = 'chr1\t100000\t101000\t+\tgene1.1\t0\n'+\
        #             'chr1\t100100\t101000\t+\tgene1.2\t0\n'
                     
        #gene2coord = 'chr2\t100000\t101000\t-\tgene2.1\t0\n'+\
        #             'chr2\t100100\t101000\t-\tgene2.2\t0\n'
                     
        #tss_coord  = ''
        
        #self.gene1bed = BedTool(gene1coord, from_string=True)
        #self.gene1bed = BedTool(gene1coord, from_string=True)
        #self.gene2bed = BedTool(gene2coord, from_string=True)
        #self.tss_bed  = BedTool(tss_coord,  from_string=True)
        
        self.hnf1b = 'HNF1B'
        self.hnf4a = 'HNF4A'
        self.unknown_gene = 'unknown gene symbol'
        
        self.hnf1b_bed = "chr17\t36046433\t36093736\tuc021tvu.1\t6\t-\n"+\
                         "chr17\t36046433\t36105096\tuc010wdi.2\t9\t-\n"+\
                         "chr17\t36046433\t36105096\tuc002hok.4\t9\t-\n"+\
                         "chr17\t36046433\t36105096\tuc021tvv.1\t8\t-\n"+\
                         "chr17\t36046433\t36105096\tuc021tvw.1\t7\t-\n"+\
                         "chr17\t36090480\t36093814\tuc010cve.1\t2\t-\n"
                         
        self.hnf4a_bed = "chr20\t42984440\t43036115\tuc010zwo.1\t4\t+\n"+\
                         "chr20\t42984440\t43053276\tuc002xlt.3\t8\t+\n"+\
                         "chr20\t42984440\t43061485\tuc002xlu.4\t10\t+\n"+\
                         "chr20\t42984440\t43061485\tuc002xlv.4\t10\t+\n"+\
                         "chr20\t43029895\t43053276\tuc002xly.4\t8\t+\n"+\
                         "chr20\t43029895\t43061485\tuc010ggq.4\t11\t+\n"+\
                         "chr20\t43029895\t43061485\tuc002xlz.4\t10\t+\n"+\
                         "chr20\t43029895\t43061485\tuc002xma.4\t10\t+\n"
       
        
      
     
    #def test_get_genes_bed_empty_genelist(self):
    #    self.assertEqual(BedsProcessor.get_genes_bed([], self.genome), [])
         
         
    def test_get_genes_bed_output_content(self):
        
        with patch("remus.processing.g") as g:
            
            g.genes_registry = self.genes_reg
            
            hnf1b = BedsProcessor.get_genes_bed([self.hnf1b], self.genome)           
            self.assertEqual(len(hnf1b), 1)
            self.assertEqual(''.join([str(i) for i in hnf1b]), self.hnf1b_bed)
            
            two_genes = BedsProcessor.get_genes_bed([self.hnf1b, self.hnf4a], self.genome)
            self.assertEqual(len(two_genes), 1)
            self.assertEqual(''.join([str(i) for i in two_genes]), self.hnf1b_bed + self.hnf4a_bed)


                
           
    def test_get_genes_bed_unknown_gene(self):
        
        with patch("remus.processing.g") as g:
            
            g.genes_registry = self.genes_reg 

            bed = BedsProcessor.get_genes_bed([self.unknown_gene], self.genome)

            self.assertEqual(len(bed),1)
            self.assertEqual(bed[0].count(), 0)
        
           
    
    def test_extraction_of_genes_tss(self):
       
        # pybedtools does not close file_handles immediately, which yields a stack of warnings
        # ignore them!
        warnings.simplefilter("ignore", ResourceWarning)
       
        upstreams, downstreams = [0,1,10,200], [0,1,10,200]
       
        with patch("remus.processing.g") as g:
            
            g.genes_registry = self.genes_reg 
            
            genes = ['HNF1B']
            for upstream in upstreams:
                for downstream in downstreams:
            
                    bed = BedsProcessor._get_joined_flanked_genes(genes, self.genome, upstream, downstream)
            
                    self.assertEqual(len(bed), len((self.hnf1b_bed.strip()).split('\n')))
                    
                    intervals = [str(i).split("\t") for i in bed]
                    for i in intervals:
                        self.assertEqual(int(i[2])-int(i[1]), upstream + downstream)
                        self.assertEqual(len(i), 6)
            
            
            genes = ['HNF1B', 'HNF4A']
            for upstream in upstreams:
                for downstream in downstreams:
            
                    bed = BedsProcessor._get_joined_flanked_genes(genes, self.genome, upstream, downstream)
            
                    self.assertEqual(len(bed), len((self.hnf1b_bed + self.hnf4a_bed).strip().split('\n')))
            
                    intervals = [str(i).split("\t") for i in bed]
                    for i in intervals:
                        self.assertEqual(int(i[2])-int(i[1]), upstream + downstream)
                        self.assertEqual(len(i), 6)

                
            


