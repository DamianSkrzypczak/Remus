

import unittest
import warnings
#from unittest.mock import patch, patch

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
        
        self.genes_reg = GenesDBRegistry.get_instance()
        
        self.genome = 'hg19'
               
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
                         
        gene1coord = 'chr1\t100000\t101000\tgene1.1\t0\t+\n'+\
                     'chr1\t100100\t101000\tgene1.2\t0\t+'
                     
        gene2coord = 'chr2\t100000\t101000\tgene2.1\t0\t-\n'+\
                     'chr2\t100100\t101000\tgene2.2\t0\t-'
                             
        self.gene1bed = BedTool(gene1coord, from_string=True)
        self.gene2bed = BedTool(gene2coord, from_string=True)

        
    def _bed2str(self, bed):
        return "\n".join([str(i).strip() for i in bed])
        
    #################################
    #
    #   Gene BED extraction tests
    #
    ################
      
      
    #
    # throws IndexError. Handle it as exception?
    #     
    #def test_get_genes_bed_empty_genelist(self):
    #    self.assertEqual([], BedsProcessor.get_genes_bed([], self.genome))
         
         
    def test_get_genes_bed_output_content(self):
        
        # pybedtools does not close file_handles immediately, which yields a stack of warnings
        # ignore them!        
        warnings.simplefilter("ignore", ResourceWarning)

        hnf1b = BedsProcessor.get_genes_bed([self.hnf1b], self.genome)[0]           
        self.assertEqual(len(hnf1b), 6)
        self.assertEqual(''.join([str(i) for i in hnf1b]), self.hnf1b_bed)

        
        two_genes = BedsProcessor.get_genes_bed([self.hnf1b, self.hnf4a], self.genome)[0]
        self.assertEqual(len(two_genes), 14)
        self.assertEqual(''.join([str(i) for i in two_genes]), self.hnf1b_bed + self.hnf4a_bed)



                
           
    def test_get_genes_bed_unknown_gene(self):
        bed = BedsProcessor.get_genes_bed([self.unknown_gene], self.genome)[0]
        self.assertEqual(len(bed), 0)
            
        
        
        
        
    ###########################
    #
    #   Gene flanking tests
    #
    ####################
    
    def _test_get_gene_promoter_sites(self, genes, expected_bed):

        # pybedtools does not close file_handles immediately, which yields a stack of warnings
        # ignore them!        
        warnings.simplefilter("ignore", ResourceWarning)
        
        upstreams, downstreams = [0,1,10,200], [0,1,10,200]
    
        for upstream in upstreams:
            for downstream in downstreams:
        
                bed = BedsProcessor._get_gene_promoter_sites(genes, self.genome, upstream, downstream)
        
                self.assertEqual(len(bed), len((expected_bed.strip()).split('\n')))
                
                intervals = [str(i).split("\t") for i in bed]
                for i in intervals:
                    self.assertEqual(upstream + downstream, int(i[2])-int(i[1]))
                    self.assertEqual(6, len(i))
            
    #
    # throws IndexError. Handle it as exception?
    #
    #def test_get_gene_promoter_sites_empty_list(self):
    #    self.assertEqual([], self._test_get_gene_promoter_sites([], ""))
    
    
    def test_get_gene_promoter_sites_single_gene(self):
        self._test_get_gene_promoter_sites([self.hnf1b], self.hnf1b_bed)
       
       
    def test_get_gene_promoter_sites_list_of_genes(self):
        self._test_get_gene_promoter_sites([self.hnf1b, self.hnf4a], self.hnf1b_bed + self.hnf4a_bed)
        
    
    
    
    ###########################
    #
    #   Interval overlapping tests
    #
    ####################



    def test_process_with_overlapping_combine_mode_all(self):
        
        mode = "all"
        
        beds = [self.gene1bed]
        overlap = BedsProcessor._combine_beds(beds, mode)
        # nothing happens if the list has a single bed
        self.assertEqual(2, len(overlap))
        self.assertEqual(self._bed2str(self.gene1bed), self._bed2str(overlap))

        # nothing in common between the two beds
        beds = [self.gene1bed, self.gene2bed]
        overlap = BedsProcessor._combine_beds(beds, mode)
        self.assertEqual(0, len(overlap))
        self.assertEqual("", self._bed2str(overlap))

        # shared small interval
        gene3=BedTool('chr1\t100300\t100310\tgene3.1\t0\t+', from_string=True)
        
        beds = [self.gene1bed, gene3]
        overlap = BedsProcessor._combine_beds(beds, mode)
        self.assertEqual(2, len(overlap))
        self.assertEqual("chr1\t100300\t100310\tgene1.1\t0\t+\n" +\
                         "chr1\t100300\t100310\tgene1.2\t0\t+", self._bed2str(overlap))
        
        # the operation is not symmetric, unless columns 4+ are dropped!!
        beds = [gene3, self.gene1bed]
        overlap = BedsProcessor._combine_beds(beds, mode)
        self.assertEqual(2, len(overlap))
        self.assertEqual("chr1\t100300\t100310\tgene3.1\t0\t+\n"+\
                         "chr1\t100300\t100310\tgene3.1\t0\t+", self._bed2str(overlap))
        
        
    def test_process_with_overlapping_combine_mode_any(self):
        mode = "any"

        beds = [self.gene1bed]
        overlap = BedsProcessor._combine_beds(beds, mode)
        # nothing happens if the list has a single bed
        self.assertEqual(2, len(overlap))
        self.assertEqual(self._bed2str(self.gene1bed), self._bed2str(overlap))

        overlap = BedsProcessor._combine_beds(beds, mode, merge=True)
        # merges overlapping intervals
        self.assertEqual(1, len(overlap))
        self.assertEqual('chr1\t100000\t101000', self._bed2str(overlap))


        beds = [self.gene1bed, self.gene2bed]
        overlap = BedsProcessor._combine_beds(beds, mode)
        self.assertEqual(4, len(overlap))
        self.assertEqual('\n'.join([self._bed2str(self.gene1bed), self._bed2str(self.gene2bed)]), self._bed2str(overlap))
        
        overlap = BedsProcessor._combine_beds(beds, mode, merge=True)
        self.assertEqual(2, len(overlap))
        expected_output='chr1\t100000\t101000\n'+\
                        'chr2\t100000\t101000'
        self.assertEqual(expected_output, self._bed2str(overlap))
        
        
        
    #
    # TODO
    # Process_with_overlapping returns an BedMutualOperation object, or an ampty list
    # Make the result type coherent
    #
    def test_process_with_overlapping_unknown_combine_mode(self):
        mode = 'unknown'
        self.assertEqual([], BedsProcessor._combine_beds(None, mode))
    

