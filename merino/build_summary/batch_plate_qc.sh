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

echo CONFIG_ROOT = "${CONFIG_ROOT}"
echo PROJECT_CODE = "${PROJECT_CODE}"

BUILD_FOLDER="${CONFIG_ROOT}${PROJECT_CODE}/build"
echo BUILD_FOLDER = "${BUILD_FOLDER}"


IFS=',' read -r -a plates <<< "${PLATES}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
PLATE="${plates[${batch_index}]}"
echo PLATE = "${PLATE}"

IFS='_' read -r -a plate_token <<< "${PLATE}";

PROJECT_DIR="${CONFIG_ROOT}${PROJECT_CODE}/${plate_token[0]}_${plate_token[1]}_${plate_token[2]}"
echo PROJECT_DIR = "${PROJECT_DIR}"

QC_FOLDER="${CONFIG_ROOT}${PROJECT_CODE}/qc"
echo QC_FOLDER = "${QC_FOLDER}"

mkdir -p "${QC_FOLDER}"
# Activate conda environment

source activate merino

cd /cmap/merino/

python setup.py develop

python /cmap/merino/merino/build_summary/plate_summary.py -project_folder ${PROJECT_DIR} -qc_folder ${QC_FOLDER} -plate_name ${PLATE}
exit_code=$?

source deactivate
exit $exit_code

