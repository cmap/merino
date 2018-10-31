#!/usr/bin/env bash

while [[ $# > 1 ]]
do

key="$1"

case $key in
    -config_root)
    CONFIG_ROOT="$2"
    shift # past argument
    ;;
    -project_dir)
    PROJECT_DIR="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

echo CONFIG_ROOT = "${CONFIG_ROOT}"
echo PROJECT_DIR = "${PROJECT_DIR}"

IFS=',' read -r -a plates <<< "${PLATES}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
PLATE_NAME="${plates[${batch_index}]}"

echo PLATE_NAME = "${PLATE_NAME}"

# Activate conda environment

source activate merino

cd /cmap/

python merino/setup_merino.py develop

python merino/normalize/card.py -proj_dir ${PROJECT_DIR} -plate_name ${PLATE_NAME}

source deactivate