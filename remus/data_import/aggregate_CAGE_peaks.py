##
# Based on FANTOM5 (F5) ontology (arg1) the script groups individual F5 samples into preselected organs, tissues and celltypes (facets).
# Only primary cells and tissue samples from human are used. 
# Robust CAGE peaks (TPM>10) for samples are extracted from expression matrix (arg2) and aggregated over the facets. 
# Output is saved in user-specified path (arg3) in facet BED files containing genomic location of CAGE peaks and aggregated score for the facet.
#
# Usage: 
# aggregate_CAGE_peaks.py  OBO_FILE  TSS_EXPRESSION_MATRIX_FILE  OUTPUT_PATH
#
# Pawel Sztromwasser, 9/2018
#


import obonet, networkx
import gzip
import sys, os

class F5Ontology():
    
    # organs used in SlideBase
    SLIDEBASE_ORGANS = set(['UBERON:0000029','UBERON:0000059','UBERON:0000178','UBERON:0000341','UBERON:0000473','UBERON:0000945',
            'UBERON:0000948','UBERON:0000955','UBERON:0000970','UBERON:0000989','UBERON:0000992','UBERON:0000995','UBERON:0000996',
            'UBERON:0001013','UBERON:0001043','UBERON:0001044','UBERON:0001134','UBERON:0001135','UBERON:0001255','UBERON:0001264',
            'UBERON:0001723','UBERON:0001736','UBERON:0001831','UBERON:0001981','UBERON:0001987','UBERON:0002046','UBERON:0002048',
            'UBERON:0002097','UBERON:0002106','UBERON:0002107','UBERON:0002108','UBERON:0002110','UBERON:0002113','UBERON:0002240',
            'UBERON:0002331','UBERON:0002360','UBERON:0002367','UBERON:0002370','UBERON:0002372','UBERON:0003112','UBERON:0004054'])
               
    ENCODE_ORGANS = set(['UBERON:0000006','UBERON:0000059','UBERON:0000473','UBERON:0000945','UBERON:0000948','UBERON:0000955',
            'UBERON:0000966','UBERON:0000970','UBERON:0000992','UBERON:0000995','UBERON:0000996','UBERON:0001114','UBERON:0001150',
            'UBERON:0001157','UBERON:0001159','UBERON:0001211','UBERON:0001224','UBERON:0001264','UBERON:0001323','UBERON:0001383',
            'UBERON:0001496','UBERON:0001499','UBERON:0001515','UBERON:0001621','UBERON:0001723','UBERON:0001774','UBERON:0001870',
            'UBERON:0001875','UBERON:0001987','UBERON:0002037','UBERON:0002046','UBERON:0002048','UBERON:0002080','UBERON:0002084',
            'UBERON:0002106','UBERON:0002107','UBERON:0002108','UBERON:0002113','UBERON:0002129','UBERON:0002167','UBERON:0002168',
            'UBERON:0002190','UBERON:0002240','UBERON:0002324','UBERON:0002331','UBERON:0002367','UBERON:0002369','UBERON:0002370',
            'UBERON:0003662','UBERON:0003663','UBERON:0004264','UBERON:0004538','UBERON:0004539','UBERON:0004550','UBERON:0004648',
            'UBERON:0005270','UBERON:0006631','UBERON:0006920','UBERON:0007610','UBERON:0008367','UBERON:0008450','UBERON:0008952',
            'UBERON:0010414','UBERON:0011907','UBERON:0018115','UBERON:0018116','UBERON:0018117','UBERON:0018118','UBERON:0036149'])

    OTHER_ORGANS = set()
    
    SLIDEBASE_CELLTYPES = set(['CL:0000047','CL:0000056','CL:0000062','CL:0000067','CL:0000071','CL:0000077','CL:0000080','CL:1000487',
            'CL:0000084','CL:0000094','CL:0000097','CL:0000098','CL:0000127','CL:0000134','CL:0000136','CL:0000138','CL:0000148',
            'CL:0000182','CL:0000188','CL:0000235','CL:0000312','CL:0000359','CL:0000388','CL:0000451','CL:0000499','CL:0000540',
            'CL:0000558','CL:0000575','CL:0000576','CL:0000622','CL:0000623','CL:0000632','CL:0000669','CL:0000731','CL:0000746',
            'CL:0000767','CL:0000775','CL:0000945','CL:0002138','CL:0002166','CL:0002224','CL:0002231','CL:0002252','CL:0002327',
            'CL:0002334','CL:0002363','CL:0002367','CL:0002368','CL:0002504','CL:0002518','CL:0002536','CL:0002548','CL:0002549',
            'CL:0002550','CL:0002552','CL:0002554','CL:0002556','CL:0002557','CL:0002559','CL:0002563','CL:0002565','CL:0002577',
            'CL:0002586','CL:0002598','CL:0002599','CL:0002600','CL:0002601','CL:0002620','CL:0002621','CL:1000306','CL:1000398'])

    ENCODE_CELLTYPES = set(['CL:0000056','CL:0000084','CL:0000127','CL:0000188','CL:0000236','CL:0000307','CL:0000312',
            'CL:0000351','CL:0000515','CL:0000545','CL:0000546','CL:0000623','CL:0000624','CL:0000625','CL:0000706','CL:0000746',
            'CL:0000765','CL:0000775','CL:0000815','CL:0000895','CL:0000899','CL:0001054','CL:0001059','CL:0002188','CL:0002231',
            'CL:0002252','CL:0002304','CL:0002306','CL:0002327','CL:0002328','CL:0002399','CL:0002518','CL:0002536','CL:0002539',
            'CL:0002547','CL:0002548','CL:0002550','CL:0002551','CL:0002552','CL:0002553','CL:0002555','CL:0002557','CL:0002558',
            'CL:0002565','CL:0002584','CL:0002586','CL:0002590','CL:0002603','CL:0002604','CL:0002606','CL:0002618','CL:0002620',
            'CL:0010001','CL:1001568','CL:1001606','CL:1001608','CL:2000010','CL:2000012','CL:2000013','CL:2000014','CL:2000016',
            'CL:2000017','CL:2000041','CL:2000043','CL:2000044','CL:2000045','NTR:0004646','NTR:0004647'])

    OTHER_CELLTYPES = set()
    
    ORGAN_IDS = SLIDEBASE_ORGANS | ENCODE_ORGANS | OTHER_ORGANS
    CELLTYPE_IDS = SLIDEBASE_CELLTYPES | ENCODE_CELLTYPES | OTHER_CELLTYPES
    
    HOMO_SAMPIENS_ID = 'NCBITaxon:9606'
    CELL_LINE_SAMPLE_ID='FF:0000003'
    
    def __init__(self, obo_file_path):
        self.g = obonet.read_obo(obo_file_path)
        self._id2name = {id_: data['name'] for id_, data in self.g.nodes(data=True)}
        self._name2id = {self._id2name[k] : k for k in self._id2name.keys()}
        
        # slim the ontology for the purpose of classifying samples to organs (removes unneeded relations)
        edge_labels_to_remove=['develops_from', 'is_model_for', 'treated_with']
        e2rm = [(a,b) for (a,b,c) in self.g.edges(keys=True) if c in edge_labels_to_remove]
        import copy
        self.slim_g = copy.copy(self.g)
        self.slim_g.remove_edges_from(e2rm)
        self.slim_g.remove_edges_from(e2rm) # to be certain that we remove them
        #print('remaining:', set([c for (a,b,c) in slim_g.edges(keys=True)]))
        
    def id2name(self, ids):
        if isinstance(ids, list) or isinstance(ids, set):
            return [self.id2name(identifier) for identifier in ids]
        return self._id2name[ids] if ids in self._id2name else None
    
    def name2id(self, names):
        if isinstance(names, list) or isinstance(names, set):
            return [self.name2id(name) for name in names]
        return self._name2id[names] if names in self._name2id else None
    
    def get_samples_for_terms(self, term_ids):
        """ Returns a dictionary organ_id:[samples] build of organ_ids given as argument,
            and list of samples belonging to that organ (according to ontology) """
        
        human_samples = networkx.ancestors(self.slim_g, F5Ontology.HOMO_SAMPIENS_ID)
        cell_line_samples = networkx.ancestors(self.slim_g, F5Ontology.CELL_LINE_SAMPLE_ID)
        
        samples={}
        for term_id in term_ids:
            if term_id not in self.slim_g:
                samples[term_id] = []
            else:
                ancestor_terms = networkx.ancestors(self.slim_g, term_id)          # take all celltypes composing an organ
                samples[term_id] = list((set(human_samples) - set(cell_line_samples)) & set(ancestor_terms)) # exclude cell_line samples and intersect
                # remove non-root terms
                samples[term_id] = [s for s in samples[term_id] if len(networkx.ancestors(self.slim_g, s))==0]
 
        return samples

    def get_organ_for_sample(self, sample_id, allow_missing = False):
        return self._get_term_for_sample(sample_id, F5Ontology.ORGAN_IDS)
    
    def get_celltypes_for_samples(self, sample_id, allow_missing = False):
        return self._get_term_for_sample(sample_id, F5Ontology.CELLTYPE_IDS)

    def _get_terms_for_sample(self, sample_id, term_ids, allow_missing = False):
        """ Returns a list of term_ids (organs / celltypes). One sample can be part of several terms. """    
        if sample_id not in self.slim_g:
            if allow_missing:
                descendants = set()
            else:
                raise Exception("Node %s is missing from the ontology." % sample_id)
        else:
            descendants = networkx.descendants(self.slim_g, sample_id)
         
        return list(descendants & term_ids)


def get_sample_ids_from_expression_table(file_handle):  

    sample_ids = []
    for l in file_handle:
        if l.find('##')==0: continue
        elif l.find('00Annotation')==0:
            header = l.split('\t')
            sample_ids = [ (e.split('.')[3]).strip() for e in header[1:] ]
            break
        else: raise Exception('No header found') # (should not happen)
    
    return sample_ids


def decode_genomic_location(field):
    s = field.split(":")
    chrom = s[0]
    s = s[1].split("..")
    start = s[0]
    s = s[1].split(",")
    end = s[0]
    strand = s[1]
    return chrom, start, end, strand


def mean(l):
    return sum(l)/len(l)


#### MAIN ####

if __name__ == '__main__':

    OBO_FILE = sys.argv[1]
    EXPRESSION_TABLE_FILE = sys.argv[2]
    OUTPUT_DIR = sys.argv[3]

    EXPRESSION_CUTOFF = 10
    EXPRESSION_AGGREGATE_FUNCTION = mean
    
    print("Reading ontology file...")
    f5o = F5Ontology(OBO_FILE)
    tsd = f5o.get_samples_for_terms(F5Ontology.ORGAN_IDS | F5Ontology.CELLTYPE_IDS)  ## term-samples dict
    
    # Prepare list of output files, one per organ/celltype. 
    # Organ/celltypes missing from the ontology are skipped
    print("Initiating output BED files...")
    bed_names = {t: os.path.join(OUTPUT_DIR,
                                 t.replace(":", "_") + "_" + f5o.id2name(t).replace(" ", "_") + ".bed"
                                 ) for t in tsd if len(tsd[t]) > 0}
    bed_files = {t: open(bed_names[t], 'wt') for t in bed_names}
    
    # iterate over expression table, and append single records to organ BED files.
    print("Parsing CAGE expression matrix...")
    with gzip.open(EXPRESSION_TABLE_FILE, 'rt') as f: 
    
        # jump to header and read sample IDs
        sample_ids = get_sample_ids_from_expression_table(f)
        
        # Create a dict of indices. 
        # Some sample IDs have replicates, so the dict values are lists
        #
        s_col={}
        for i, s_id in enumerate(sample_ids):
            if s_id in s_col:
                s_col[s_id] = s_col[s_id] + [i]
            else:
                s_col[s_id] = [i]
        
        #s_col = {sample_ids[i] : (i+1) for i in range(1,len(sample_ids))} 
    
        # skip two lines of normalization stats
        assert f.readline().find("01STAT:MAPPED") == 0
        assert f.readline().find("02STAT:NORM_FACTOR") == 0
        
        row_cnt=0
        for l in f:
            lsplit = l.split('\t')
            chrom, start, end, strand = decode_genomic_location(lsplit[0])
            
            new_record = '\t'.join([chrom, start, end, strand])
            
            if row_cnt%100==0: print(new_record)
            row_cnt+=1
            
            for facet, samples in tsd.items():
                
                # skip if there is no samples for term/facet
                if len(samples) == 0: continue
                
                s_ids = [s[len('FF:'):] for s in samples]
                
                # remove ontology sample_ids missing from expression matrix (isolated cases)
                missing = [s for s in s_ids if s not in s_col]
                if len(missing)>0:
                    #print("Following sample IDs for facet [%s] (%s) are missing from expression matrix: %s" % (f5o.id2name(facet), facet, str(["FF:"+s for s in missing])))
                    s_ids = [s for s in s_ids if s not in missing]
                    
                # extract promoter expression values for samples in this facet
                # values for technical replicates (list of columns with the same extract ID) are aggregated separately, to not bias the facet score           
                expr_list = []
                for s in s_ids:
                    cols = s_col[s]
                    if len(cols)>1:    # technical replicates 
                        exprs = [float(lsplit[c]) for c in cols]
                        expr_list.append(EXPRESSION_AGGREGATE_FUNCTION(exprs))
                    else:
                        expr_list.append(float(lsplit[cols[0]]))           
                
                # calculate aggregated expression/activity of the promoter in the facet
                # provided that at least one of samples meets the cutoff criteria
                if max(expr_list) >= EXPRESSION_CUTOFF:
                    #expression = sum(expression_list)/len(expression_list)
                    score = EXPRESSION_AGGREGATE_FUNCTION(expr_list)
                    bed_files[facet].write(new_record + ('\t%.2f\n' % score))
        
    for _,f in bed_files.items():
        f.close()
    
    print("Done.")
