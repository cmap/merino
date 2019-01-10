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

PROJECT_DIR="${CONFIG_ROOT}${PROJECT_CODE}"

echo CONFIG_ROOT = "${CONFIG_ROOT}"
echo PROJECT_DIR = "${PROJECT_DIR}"

IFS=',' read -r -a replicate_sets <<< "${REPLICATE_SETS}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
REPLICATE_SET_NAME="${replicate_sets[${batch_index}]}"

echo REPLICATE_SET_NAME = "${REPLICATE_SET_NAME}"

# Activate conda environment

source activate merino

cd /cmap/

python merino/setup_merino.py develop

python /cmap/merino/normalization/weave.py -proj_dir ${PROJECT_DIR} -replicate_set_name ${REPLICATE_SET_NAME} -all_inputs -aggregate_output
exit_code=$?

source deactivate
exit $exit_code