#!/bin/bash

set -e
set -u

ACTIVATE_DIR=$PREFIX/etc/conda/activate.d
DEACTIVATE_DIR=$PREFIX/etc/conda/deactivate.d

cp -v $RECIPE_DIR/activate.sh $ACTIVATE_DIR/jvarkit-activate.sh
cp -v $RECIPE_DIR/deactivate.sh $DEACTIVATE_DIR/jvarkit-deactivate.sh

mkdir -p $PREFIX/dist
./gradlew  -Djvarkit.disable.test=true vcfhead vcftail vcffilterjdk vcf2table samjdk bioalcidaejdk vcffilterso sam2tsv mkminibam
mv -v dist/*.jar $PREFIX/dist/


