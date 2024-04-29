#!/bin/bash

: <<'END_COMMENT'

This script builds V(D)J reference in igblast format

END_COMMENT

VALID_ARGS=$(getopt -o a --long allow-minor-alleles -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --allow-minor-alleles)
        echo "Enabled --allow-minor-alleles"
        save_all_alleles=true
        shift
        ;;
    --) shift;
        break
        ;;
  esac
done

OUTPUT_DIR=/tmp

IGBLAST_DIR=${OUTPUT_DIR}/ncbi-igblast
IGBLAST_VERSION=1.22.0
wget -q https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ${OUTPUT_DIR}/ncbi-igblast.tar.gz && \
tar -xvzf ${OUTPUT_DIR}/ncbi-igblast.tar.gz --one-top-level=${IGBLAST_DIR} --strip-component 1 && \
rm ${OUTPUT_DIR}/ncbi-igblast.tar.gz

SEQKIT_PATH=${OUTPUT_DIR}/seqkit
SEQKIT_VERSION=2.8.1
wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O ${OUTPUT_DIR}/seqkit.tar.gz && \
tar -xvzf ${OUTPUT_DIR}/seqkit.tar.gz --one-top-level=${SEQKIT_PATH} --strip-component 1 && \
rm ${OUTPUT_DIR}/seqkit.tar.gz

get_organism_specie() {
  if [[ $1 == "human" ]]; then
    echo "Homo_sapiens"
  else
    echo "Mus_musculus"
  fi
}

get_receptor_name() {
  if [[ ${1:0:1} == "I" ]]; then
    echo "Ig"
  else
    echo "TCR"
  fi
}

REF_DIR=${OUTPUT_DIR}/igblast.reference

for ORGANISM in human mouse
do
  for LOCUS in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ TRBV TRBD TRBJ TRAV TRAJ TRDV TRDD TRDJ TRGV TRGJ
  do
    SPECIE=$(get_organism_specie $ORGANISM)
    RECEPTOR=$(get_receptor_name $LOCUS)
    OUTPUT_NAME=${ORGANISM}.${RECEPTOR}.${LOCUS:3}

    # download reference from https://www.imgt.org/ database
    wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${SPECIE}/${LOCUS:0:2}/${LOCUS}.fasta -O ${OUTPUT_DIR}/${OUTPUT_NAME}.imgt

    # convert imgt fasta -> igblast fasta
    ${IGBLAST_DIR}/bin/edit_imgt_file.pl ${OUTPUT_DIR}/${OUTPUT_NAME}.imgt > ${OUTPUT_DIR}/${OUTPUT_NAME}.all.fasta

    if [ $save_all_alleles ] ; then
      mv ${OUTPUT_DIR}/${OUTPUT_NAME}.all.fasta ${REF_DIR}/database/${OUTPUT_NAME}
    else
      # select only minor alleles *01
      $SEQKIT_PATH grep ${OUTPUT_DIR}/${OUTPUT_NAME}.all.fasta -r -p "\*01" -o ${REF_DIR}/database/${OUTPUT_NAME}
      rm ${OUTPUT_DIR}/${OUTPUT_NAME}.all.fasta
    fi

    # make blast db index
    ${IGBLAST_DIR}/bin/makeblastdb -parse_seqids -dbtype nucl -in ${REF_DIR}/database/${OUTPUT_NAME}

    # remove unnecessary fasta files
    rm ${OUTPUT_DIR}/${OUTPUT_NAME}.imgt ${REF_DIR}/database/${OUTPUT_NAME}
  done
done

# add another necessary reference files
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar -q -O ${OUTPUT_DIR}/ncbi_human_c_genes.tar && \
tar -xvf ${OUTPUT_DIR}/ncbi_human_c_genes.tar -C ${REF_DIR}/database && \
cp -r ${IGBLAST_DIR}/internal_data ${REF_DIR} && \
cp -r ${IGBLAST_DIR}/optional_file ${REF_DIR} && \
tar czf ${REF_DIR}.tar.gz -C ${REF_DIR} database internal_data optional_file  # create archive with V(D)J ref

echo "Archive with IgBLAST reference here ${REF_DIR}.tar.gz"
