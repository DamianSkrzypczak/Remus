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
    
    OTHER_ORGANS = set()
    
    ORGAN_IDS = SLIDEBASE_ORGANS | OTHER_ORGANS
    
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
            return [self._id2name[identifier] for identifier in ids]
        return self._id2name[ids]
    
    def name2id(self, names):
        if isinstance(names, list) or isinstance(names, set):
            return [self._name2id[name] for name in names]
        return self._name2id[name]
    
    def get_samples_for_organs(self, organs):
        """ Returns a dictionary organ_id:[samples] build of organ_ids given as argument,
            and list of samples belonging to that organ (according to ontology) """
        organ_samples={}
        for organ_id in organs:
            human_samples = networkx.ancestors(self.slim_g, F5Ontology.HOMO_SAMPIENS_ID)
            cell_line_samples = networkx.ancestors(self.slim_g, F5Ontology.CELL_LINE_SAMPLE_ID)
            celltypes = networkx.ancestors(self.slim_g, organ_id)          # take all celltypes composing an organ
            organ_samples[organ_id] = list((set(human_samples) - set(cell_line_samples)) & set(celltypes)) # exclude cell_line samples and intersect
            # remove non-root terms
            organ_samples[organ_id] = [s for s in organ_samples[organ_id] if len(networkx.ancestors(self.slim_g, s))==0]
 
        return organ_samples
    
    
    def get_organ_for_sample(self, sample_id, allow_missing = False):
        """ Returns a list of organ_ids. One sample can be part of several organs. """    
        try:
            descendants = networkx.descendants(self.slim_g, sample_id)
        except networkx.exception.NetworkXError as e:
            if allow_missing: 
                return []
            else:
                raise Error("Node %s is missing from the ontology.\n Original exception: %s" % (sample_id, str(e)))
 
        return list(descendants & F5Ontology.ORGAN_IDS)

        
        
        

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
    osd = f5o.get_samples_for_organs(F5Ontology.ORGAN_IDS)  ## organ-samples dict
    
    
    print("Initiating output BED files...")
    # prepare list of output files, one per organ
    organ_bed_names = {o : os.path.join(OUTPUT_DIR, "_".join([o,f5o.id2name(o).replace(" ","_"),"promoters.bed"])) for o in osd.keys()}
    organ_bed_files = {o : open(organ_bed_names[o], 'wt') for o in osd.keys()}
    
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
            
            for organ, samples in osd.items():
                
                s_ids = [s[len('FF:'):] for s in samples]
                
                # remove ontology sample_ids missing from expression matrix (isolated cases)
                missing = [s for s in s_ids if s not in s_col]
                if len(missing)>0:
                    #print("Following sample IDs for organ [%s] (%s) are missing from expression matrix: %s" % (f5o.id2name(organ), organ, str(["FF:"+s for s in missing])))
                    s_ids = [s for s in s_ids if s not in missing]
                    
                # extract promoter expression values for samples in this organ
                # values for technical replicates (list of columns with the same extract ID) are aggregated separately, to not bias the organ score           
                expr_list = []
                for s in s_ids:
                    cols = s_col[s]
                    if len(cols)>1:    # technical replicates 
                        exprs = [float(lsplit[c]) for c in cols]
                        expr_list.append(EXPRESSION_AGGREGATE_FUNCTION(exprs))
                    else:
                        expr_list.append(float(lsplit[cols[0]]))           
                
                # calculate aggregated expression/activity of the promoter in the organ
                # provided that at least one of samples meets the cutoff criteria
                if max(expr_list) >= EXPRESSION_CUTOFF:
                    #expression = sum(expression_list)/len(expression_list)
                    score = EXPRESSION_AGGREGATE_FUNCTION(expr_list)
                    organ_bed_files[organ].write(new_record + ('\t%.2f\n' % score))
        
    for _,f in organ_bed_files.items():
        f.close()
    
    print("Done.")
    
