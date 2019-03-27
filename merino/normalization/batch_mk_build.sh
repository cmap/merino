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
    -cohort_name)
    COHORT_NAME="$2"
    shift
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
echo COHORT_NAME = "${COHORT_NAME}"

# Activate conda environment
BUILD_FOLDER="${PROJECT_DIR}/build"
mkdir -p "${BUILD_FOLDER}"

source activate merino

cd /cmap/merino/

python setup.py develop

python /cmap/merino/merino/normalization/mk_build_file.py -proj_dir ${PROJECT_DIR} -cohort_name ${COHORT_NAME} -build_folder ${BUILD_FOLDER} -aggregate_out
exit_code=$?

source deactivate
exit $exit_code