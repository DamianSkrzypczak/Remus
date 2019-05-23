#!/usr/bin/env bash


# Tree Configuration
DATA_ROOT=data

GENES=${DATA_ROOT}/genes
GENES_RAW=${GENES}/raw

MIRNA=${DATA_ROOT}/mirna
MIRNA_RAW=${MIRNA}/raw

F5_TSS_HG19=${DATA_ROOT}/tss/fantom5/hg19
F5_TSS_HG38=${DATA_ROOT}/tss/fantom5/GRCh38

F5_ENH=${DATA_ROOT}/enhancers/fantom5
F5_ENH_HG19=${F5_ENH}/hg19
F5_ENH_HG38=${F5_ENH}/GRCh38
F5_ENH_RAW=${F5_ENH}/raw

ENC_ENH=${DATA_ROOT}/enhancers/encode
ENC_ENH_RAW=${ENC_ENH}/raw

CHROMATIN=${DATA_ROOT}/chromatin
CHROMATIN_RAW=${CHROMATIN}/raw

LIFTOVER_EXEC=external_resources/liftOver
LIFTOVER_HG19_HG38_CHAIN=external_resources/hg19ToHg38.over.chain.gz

# Make tree
printf "Making directories tree under %s\n" ${DATA_ROOT}
mkdir -p ${GENES_RAW} ${GENES} -v
mkdir -p ${MIRNA_RAW} ${MIRNA} -v
mkdir -p ${F5_ENH_RAW} ${F5_ENH_HG19} ${F5_ENH_HG38} -v
mkdir -p ${F5_TSS_HG19} ${F5_TSS_HG38} -v
mkdir -p ${ENC_ENH_RAW} ${ENC_ENH} -v
mkdir -p ${CHROMATIN_RAW} ${CHROMATIN} -v


# Create genes.db
printf "Creating genes db\n"
PREDEFINED_GENES_DB_SOURCES=predefined_genes_db_sources
cp ${PREDEFINED_GENES_DB_SOURCES}/*.tsv ${GENES_RAW}/
python3 remus/data_import/create_genes_db.py -i ${GENES_RAW} -o ${GENES}/genes.db

# Create miRNA targets.db
printf "Creating miRNA targets db\n"
PREDEFINED_MIRNA_DB_SOURCES=predefined_mirna_db_sources
cp ${PREDEFINED_MIRNA_DB_SOURCES}/*.tsv.gz ${MIRNA_RAW}/
wget -O ${PREDEFINED_MIRNA_DB_SOURCES}/hsa_miRWalk_3UTR.7z http://mirwalk.umm.uni-heidelberg.de/download/hsa_miRWalk_3UTR.7z
7zr x -so ${PREDEFINED_MIRNA_DB_SOURCES}/hsa_miRWalk_3UTR.7z | gzip -c > ${MIRNA_RAW}/mirwalk_3UTR.tsv.gz
python3 remus/data_import/create_mirna_target_db.py -i ${MIRNA_RAW} -o ${MIRNA}/targets.db


# Download FANTOM5 CAGE expression matrix and ontology (to find transcription start sites)
PREDEFINED_F5_TSS_SOURCES=predefined_fantom5_tss_data_sources
printf "Acquiring transcription start sites FANTOM5 data\n"
wget -O ${PREDEFINED_F5_TSS_SOURCES}/hg19.cage_peak_phase1and2combined_tpm.osc.txt.gz 'http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm.osc.txt.gz'
wget -O ${PREDEFINED_F5_TSS_SOURCES}/ff-phase2-170801.obo.txt http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Ontology/ff-phase2-170801.obo.txt
# aggregate samples by organs, tissues and cell-types and store location of TSSs in BED files
python3 remus/data_import/aggregate_CAGE_peaks.py ${PREDEFINED_F5_TSS_SOURCES}/ff-phase2-170801.obo.txt ${PREDEFINED_F5_TSS_SOURCES}/hg19.cage_peak_phase1and2combined_tpm.osc.txt.gz ${F5_TSS_HG19}
# compress and index BED files
for b in ${F5_TSS_HG19}/*.bed; do
    bedtools sort -i ${b} > ${b}.sbed && mv ${b}.sbed ${b} && bgzip ${b} && tabix -p bed ${b}.gz
done
# liftover to hg38
for b in ${F5_TSS_HG19}/*.bed.gz; do
    hg38_bed=${F5_TSS_HG38}/`basename ${b%.gz}`
    ${LIFTOVER_EXEC} ${b} ${LIFTOVER_HG19_HG38_CHAIN} ${hg38_bed} ${hg38_bed}.unmapped
    bedtools sort -i ${hg38_bed} > ${hg38_bed}.sorted && mv ${hg38_bed}.sorted ${hg38_bed}
    bgzip ${hg38_bed} && tabix -p bed ${hg38_bed}.gz
    echo "Lifted over ${b} to ${hg38_bed}.gz"
done




# Download FANTOM5 enhancers 
printf "Acquiring enhancers fantom5 data\n"
wget http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz -P ${F5_ENH_RAW}
printf "... extracting celltype data"
tar -xzf ${F5_ENH_RAW}/facet_expressed_enhancers.tgz -C ${F5_ENH_HG19} --wildcards CL:*
printf "... extracting organ data"
tar -xzf ${F5_ENH_RAW}/facet_expressed_enhancers.tgz -C ${F5_ENH_HG19} --wildcards UBERON*
# compress and index BED files
for b in ${F5_ENH_HG19}/*.bed; do
    bgzip ${b} && tabix -p bed ${b}.gz
done
# liftover to hg38
for b in ${F5_ENH_HG19}/*.bed.gz; do
    hg38_bed=${F5_ENH_HG38}/`basename ${b%.gz}`
    ${LIFTOVER_EXEC} ${b} ${LIFTOVER_HG19_HG38_CHAIN} ${hg38_bed} ${hg38_bed}.unmapped
    bedtools sort -i ${hg38_bed} > ${hg38_bed}.sorted && mv ${hg38_bed}.sorted ${hg38_bed}
    bgzip ${hg38_bed} && tabix -p bed ${hg38_bed}.gz
    echo "Lifted over ${b} to ${hg38_bed}.gz"
done



# Download ENCODE enhancers data (ChIP-seq)
PREDEFINED_ENC_ENH_SOURCES=predefined_encode_enhancer_data_sources
printf "Acquiring ENCODE enhancers data\n"
# download raw BED files
awk -F '\t' '$43 ~ hg19 {print $42}' ${PREDEFINED_ENC_ENH_SOURCES}/ENCODE_enhancers_ChipSeq.metadata.tsv | wget -i - -P ${ENC_ENH_RAW}
# generate collapsing script & run it
python3 remus/data_import/collapse_encode_enhancer_beds.py ${PREDEFINED_ENC_ENH_SOURCES}/ENCODE_enhancers_ChipSeq.metadata.tsv ${ENC_ENH_RAW} ${ENC_ENH} > ${ENC_ENH}/collapse_and_liftover.sh
chmod u+x ${ENC_ENH}/collapse_and_liftover.sh && ${ENC_ENH}/collapse_and_liftover.sh
# delete raw BEDs to save space
# rm -r ${ENC_ENH}


# Download ENCODE accessible chromatin data (~3GB)
PREDEFINED_DNASEQ_SOURCES=predefined_dnaseq_data_sources
printf "Acquiring ENCODE accessible chromatin data\n"
# download raw BED files
awk -F '\t' '$43 ~ hg19 {print $42}' ${PREDEFINED_DNASEQ_SOURCES}/ENCODE_DNase_seq.metadata.tsv | wget -i - -P ${CHROMATIN_RAW}
# generate collapsing script & run it
python3 remus/data_import/collapse_encode_chromatin_beds.py ${PREDEFINED_DNASEQ_SOURCES}/ENCODE_DNase_seq.metadata.tsv ${CHROMATIN_RAW} ${CHROMATIN} > ${CHROMATIN}/collapse_and_liftover.sh
chmod u+x ${CHROMATIN}/collapse_and_liftover.sh && ${CHROMATIN}/collapse_and_liftover.sh
# delete raw BEDs to save space
# rm -r ${CHROMATIN_RAW}





