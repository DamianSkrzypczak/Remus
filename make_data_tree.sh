#!/usr/bin/env bash

DATA_ROOT=data
GENES_DATA=${DATA_ROOT}/genes
GENES_DB_SOURCES=${GENES_DATA}/genes_db_sources
TISSUES=${DATA_ROOT}/tissues
TISSUES_CELLTYPE=${DATA_ROOT}/tissues/celltype
TISSUES_ORGAN=${DATA_ROOT}/tissues/organ
TSS=${DATA_ROOT}/transcription_start_sites
CHROMATIN=${DATA_ROOT}/chromatin
CHROMATIN_RAW=${DATA_ROOT}/chromatin/raw


#Make dirs
printf "Making directories tree under %s\n" ${DATA_ROOT}
mkdir -p ${GENES_DB_SOURCES} -v
mkdir -p ${TISSUES_CELLTYPE} -v
mkdir -p ${TISSUES_ORGAN} -v
mkdir -p ${TSS} -v
mkdir -p ${CHROMATIN} -v
mkdir -p ${CHROMATIN_RAW} -v


# Download genes
# TODO: add file collecting for genes_db_sources instead of copying predefined files
rm -Rf ${GENES_DB_SOURCES}
cp -R predefined_genes_db_sources ${GENES_DB_SOURCES}

# create genes_db
printf "Creating genes db\n"
python create_genes_db.py -i ${GENES_DB_SOURCES} -o ${GENES_DATA}/genes.db


# Download enhancers fantom5
printf "Acquiring enhancers fantom5 data\n"
wget http://enhancer.binf.ku.dk/presets/facet_expressed_enhancers.tgz -P ${TISSUES}
printf "... extracting celltype data"
tar -xzf ${TISSUES}/facet_expressed_enhancers.tgz -C ${TISSUES_CELLTYPE} --wildcards CL:*
printf "... extracting organ data"
tar -xzf ${TISSUES}/facet_expressed_enhancers.tgz -C ${TISSUES_ORGAN} --wildcards UBERON*


# Download transcription start sites
printf "Acquiring transcription start sites fantom5 data\n"
wget -O ${TSS}/promoter_data.bed 'http://promoter.binf.ku.dk/viewer.php?match=and&sort-by=donotsort&end-site=249250621&start-site=1&chr-number=ALL&toggle=basic&return=download'

#
# Download ENCODE accessible chromatin data
PREDEFINED_DNASEQ_SOURCES=predefined_dnaseq_data_sources
printf "Acquiring ENCODE accessible chromatin data\n"
# download raw BED files
awk -F '\t' '$43 ~ hg19 {print $42}' ${PREDEFINED_DNASEQ_SOURCES}/ENCODE_DNase_seq.metadata.tsv | wget -i - -P ${CHROMATIN_RAW}
# generate collapsing script & run it
python ${PREDEFINED_DNASEQ_SOURCES}/collapse_tissue_beds.py ${PREDEFINED_DNASEQ_SOURCES}/ENCODE_DNase_seq.metadata.tsv hg19 > ${CHROMATIN}/collapse_hg19.sh
chmod u+x ${CHROMATIN}/collapse_hg19.sh && ${CHROMATIN}/collapse_hg19.sh
# delete raw BEDs to save space
#rm -r ${CHROMATIN_RAW}





