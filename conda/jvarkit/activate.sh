export CONDA_BACKUP_JVARKIT_HOME="${JVARKIT_HOME}"
export CONDA_BACKUP_JVARKIT_DIST="${JVARKIT_DIST}"


export JVARKIT_HOME=${CONDA_PREFIX}
export JVARKIT_DIST="${JVARKIT_HOME}/dist"

echo "[INFO] \${JVARKIT_DIST} is set to ${JVARKIT_DIST}. Example 'java -jar \${JVARKIT_DIST}/vcf2table.jar -h'." 1>&2
