#!/bin/bash

: <<'END_COMMENT'

This script builds archive with V(D)J reference and vidjil presets

END_COMMENT

REPOSITORY_DIR=/tmp/vidjil
VIDJIL_VERSION=release-2024.02
git clone -b $VIDJIL_VERSION https://gitlab.inria.fr/vidjil/vidjil.git $REPOSITORY_DIR &&
cd $REPOSITORY_DIR/germline && \
make germline  # generate vidjil germline

GERMLINE_DIR=${REPOSITORY_DIR}/germline

# add C-genes presets to the main preset
jq -s '.[0] * .[1]' ${GERMLINE_DIR}/homo-sapiens.g ${GERMLINE_DIR}/homo-sapiens-isotypes.g > ${GERMLINE_DIR}/homo-sapiens.full.g && \
mv ${GERMLINE_DIR}/homo-sapiens.full.g ${GERMLINE_DIR}/homo-sapiens.g

# create archive with vidjil ref
REF_PATH=/tmp/vidjil.germline.tar.gz
tar czf $REF_PATH \
  -C ${GERMLINE_DIR} gallus-gallus.g  homo-sapiens-cd.g  homo-sapiens.g \
    homo-sapiens-isotypes.g homo-sapiens-isoforms.g mus-musculus.g  rattus-norvegicus.g  sus-scrofa.g \
    homo-sapiens gallus-gallus mus-musculus rattus-norvegicus sus-scrofa

echo "Archive with vidjil reference here ${REF_PATH}"
