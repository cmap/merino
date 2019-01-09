#!/usr/bin/env bash

while [[ $# > 1 ]]
do

key="$1"

case $key in
    -config_root)
    CONFIG_ROOT="$2"
    shift # past argument
    ;;
    -project_code)
    PROJECT_CODE="$2"
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

IFS=',' read -r -a plates <<< "${PLATES}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
PLATE_NAME="${plates[${batch_index}]}"

IFS='_' read -r -a plate_token <<< "${PLATE_NAME}";

PROJECT_DIR="${CONFIG_ROOT}${PROJECT_CODE}/${plate_token[0]}_${plate_token[1]}_${plate_token[2]}"

echo CONFIG_ROOT = "${CONFIG_ROOT}"
echo PROJECT_DIR = "${PROJECT_DIR}"


echo PLATE_NAME = "${PLATE_NAME}"

# Activate conda environment

source activate merino

cd /cmap/

python merino/setup_merino.py develop

python /cmap/merino/normalization/card.py -proj_dir ${PROJECT_DIR} -plate_name ${PLATE_NAME}

source deactivate
