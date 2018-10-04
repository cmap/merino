#!/bin/bash
# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to > 0 the /etc/hosts part is not recognized ( may be a bug )


while [[ $# > 1 ]]
do
key="$1"

case $key in
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
    -plate_map_path)
    PLATE_MAP_PATH="$2"
    shift # past argument
    ;;
    -csv_filepath)
    CSV_FILEPATH="$2"
    shift # past argument
    ;;
    -davepool_id_csv_filepath_pairs)
    DAVEPOOL_ID_CSV_FILEPATH_PAIRS="$2 $3 $4 $5"
    shift
    shift
    shift
    shift
    ;;
    -outfile)
    OUTFILE="$2"
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

echo CONFIG_FILEPATH  = "${CONFIG_FILEPATH}"
echo ASSAY_TYPE     = "${ASSAY_TYPE}"
echo PERT_TIME   = "${PERT_TIME}"
echo PLATE_MAP_PATH = "${PLATE_MAP_PATH}"
echo CSV_FILEPATH   = "${CSV_FILEPATH}"
echo DAVEPOOL_ID_CSV_FILEPATH_PAIRS = "${DAVEPOOL_ID_CSV_FILEPATH_PAIRS}"
echo OUTFILE = "${OUTFILE}"
# Activate conda environment
source activate merino

cd /cmap/

python merino/setup_merino.py develop
if [ -z "$DAVEPOOL_ID_CSV_FILEPATH_PAIRS" ];
then
python /cmap/merino/assemble/assemble.py -config_filepath ${CONFIG_FILEPATH} -assay_type ${ASSAY_TYPE} -pert_time ${PERT_TIME} -pmp ${PLATE_MAP_PATH} -dp_csv ${DAVEPOOL_ID_CSV_FILEPATH_PAIRS} -out ${OUTFILE}
else
python /cmap/merino/assemble/assemble.py -config_filepath ${CONFIG_FILEPATH} -assay_type ${ASSAY_TYPE} -pert_time ${PERT_TIME} -pmp ${PLATE_MAP_PATH} -csv_filepath ${CSV_FILEPATH} -out ${OUTFILE}
fi

# Deactivate conda environment
source deactivate


