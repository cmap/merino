#!/usr/bin/env bash

while [[ $# > 1 ]]
do

key="$1"

case $key in
    -config_root)
    CONFIG_ROOT="$2"
    shift # past argument
    ;;
    -config_filepath)
    CONFIG_FILEPATH="$2"
    shift # past argument
    ;;
    -assay_type)
    ASSAY_TYPE="$2"
    shift # past argument
    ;;
    -pert_time)
    PERT_TIME="$2"
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
echo CONFIG_FILEPATH  = "${CONFIG_FILEPATH}"
echo ASSAY_TYPE     = "${ASSAY_TYPE}"
echo PERT_TIME   = "${PERT_TIME}"
echo PROJECT_CODE = "${PROJECT_CODE}"

IFS=',' read -r -a plates <<< "${PLATES}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
PLATE="${plates[${batch_index}]}"
echo "PLATE IS: ${PLATE}"

IFS='_' read -r -a plate_token <<< "${PLATE}";

PLATE_MAP_PATH="${CONFIG_ROOT}${PROJECT_CODE}/map_src/${plate_token[0]}.src"

echo PLATE_MAP_PATH = "${PLATE_MAP_PATH}"
OUTFILE="${CONFIG_ROOT}${PROJECT_CODE}/${plate_token[0]}_${plate_token[1]}_${plate_token[2]}"

echo OUTFILE = "${OUTFILE}"
# Activate conda environment

source activate merino

cd /cmap/

python merino/setup_merino.py develop

if [ "${ASSAY_TYPE}" = "DP78" ];
then
    DP7_PLATE="${plate_token[0]}_DP7_${PERT_TIME}_${plate_token[3]}_${plate_token[4]}"
    DP8_PLATE="${plate_token[0]}_DP8_${PERT_TIME}_${plate_token[3]}_${plate_token[4]}"

    DP7_CSV_PATH="${CONFIG_ROOT}${PROJECT_CODE}/lxb/${DP7_PLATE}/${DP7_PLATE}.csv"
    DP8_CSV_PATH="${CONFIG_ROOT}${PROJECT_CODE}/lxb/${DP8_PLATE}/${DP8_PLATE}.csv"
    DAVEPOOL_ID_CSV_FILEPATH_PAIRS="DP7 ${DP7_CSV_PATH} DP8 ${DP8_CSV_PATH}"
    echo "DAVEPOOL_ID_CSV_FILEPATH_PAIRS ${DAVEPOOL_ID_CSV_FILEPATH_PAIRS}"

    python /cmap/merino/assemble/assemble.py -config_filepath ${CONFIG_FILEPATH} -assay_type ${ASSAY_TYPE} -pert_time ${PERT_TIME} -pmp ${PLATE_MAP_PATH} -dp_csv ${DAVEPOOL_ID_CSV_FILEPATH_PAIRS} -out ${OUTFILE}
    exit_code=$?
else
    CSV_FILEPATH="${CONFIG_ROOT}${PROJECT_CODE}/lxb/${PLATE}/${PLATE}.jcsv"
    echo CSV_FILEPATH = "${CSV_FILEPATH}"
    python /cmap/merino/assemble/assemble.py -config_filepath ${CONFIG_FILEPATH} -assay_type ${ASSAY_TYPE} -pert_time ${PERT_TIME} -pmp ${PLATE_MAP_PATH} -csv ${CSV_FILEPATH} -out ${OUTFILE}
    exit_code=$?
fi

# Deactivate conda environment
source deactivate
exit $exit_code
