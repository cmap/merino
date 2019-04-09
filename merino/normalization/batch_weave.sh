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
    -group_by)
    GROUP_BY="$3"
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
IFS=',' read -r -a replicate_sets <<< "${REPLICATE_SETS}"

batch_index=${AWS_BATCH_JOB_ARRAY_INDEX}
REPLICATE_SET_NAME="${replicate_sets[${batch_index}]}"

PROJECT_DIR="${CONFIG_ROOT}${PROJECT_CODE}/${REPLICATE_SET_NAME}"

GROUP_BY="${GROUP_BY}"

echo CONFIG_ROOT = "${CONFIG_ROOT}"
echo PROJECT_DIR = "${PROJECT_DIR}"
echo REPLICATE_SET_NAME = "${REPLICATE_SET_NAME}"
echo GROUP_BY = "${GROUP_BY}"


# Activate conda environment

source activate merino

cd /cmap/merino/

python setup.py develop

#"-group_by", "-gb", help="Field(s) to group by for MODZ. If you are using more than one field separate columns names with a comma"
#                        ,type=str,default='pert_well'

python /cmap/merino/merino/normalization/weave.py -proj_dir ${PROJECT_DIR} -replicate_set_name ${REPLICATE_SET_NAME} -group_by ${GROUP_BY} -all_inputs -aggregate_output
exit_code=$?

source deactivate
exit $exit_code