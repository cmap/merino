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

PROJECT_DIR="${CONFIG_ROOT}${PROJECT_CODE}"
echo PROJECT_DIR = "${PROJECT_DIR}"

BUILD_FOLDER="${PROJECT_DIR}/build"
echo BUILD_FOLDER = "${BUILD_FOLDER}"

QC_FOLDER="${PROJECT_DIR}/qc"
echo QC_FOLDER = "${QC_FOLDER}"

mkdir -p "${QC_FOLDER}"
# Activate conda environment

source activate merino

cd /cmap/

python merino/setup.py develop

python /cmap/merino/build_summary/build_summary.py -build_folder ${BUILD_FOLDER} -qc_folder ${QC_FOLDER} -project_name ${PROJECT_CODE}
exit_code=$?

source deactivate
exit $exit_code

