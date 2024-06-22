#!/bin/bash

: <<'END_COMMENT'

This script downloads OLGA model files and make archive from them

END_COMMENT

REPOSITORY_DIR=/tmp/OLGA

git clone --quiet https://github.com/statbiophys/OLGA.git $REPOSITORY_DIR

OUTPUT_DIR=/tmp/olga-models
MODELS_DIR=$REPOSITORY_DIR/olga/default_models/

tar czf ${OUTPUT_DIR}.tar.gz -C $MODELS_DIR .
echo "Archive with OLGA models here ${OUTPUT_DIR}.tar.gz"
