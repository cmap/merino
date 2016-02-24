import logging
import setup_logger
import csv


logger = logging.getLogger(setup_logger.LOGGER_NAME)

datatype_names = ["Median", "Count"]


class DavepoolData(object):
    def __init__(self, csv_filepath=None, median_headers=None, median_data=None,
        count_headers=None, count_data=None, davepool_id=None):

        self.csv_filepath = csv_filepath
        self.median_headers = median_headers
        self.median_data = median_data
        self.count_headers = count_headers
        self.count_data = count_data
        self.davepool_id = davepool_id

    def __repr__(self):
        return " ".join(["{}:{}".format(k,v) for (k,v) in self.__dict__.items()])

def get_datatype_range(data, datatype_names):
    datatype_indexes = [i for (i,row) in enumerate(data) if len(row) > 0 and row[0] == "DataType:"]
    datatypes = [data[i][1] for i in datatype_indexes]
    datatype_indexes.append(len(data))

    r = {}
    for dn in datatype_names:
        if dn in datatypes:
            dn_index = datatypes.index(dn)
            r[dn] = (datatype_indexes[dn_index], datatype_indexes[dn_index+1])
        else:
            r[dn] = (None, None)

    return r


def read_data(csv_filepath):
    f = open(csv_filepath)
    reader = csv.reader(f)

    data = list()
    for row in reader:
        data.append(row)
    f.close()

    logger.debug("len(data):  {}".format(len(data)))

    datatype_ranges = get_datatype_range(data, datatype_names)

    for dn in datatype_names:
        if datatype_ranges[dn][0] is None:
            raise Exception("davepool_data read_data did not find expected DataType dn:  {}".format(dn))

    median_range = datatype_ranges[datatype_names[0]]
    count_range = datatype_ranges[datatype_names[1]]

    pd = DavepoolData(csv_filepath=csv_filepath)
    pd.median_headers = data[median_range[0]+1]
    pd.median_data = data[(median_range[0]+2):(median_range[1]-1)]
    logger.debug("first line of median data - pd.median_data[0]:  {}".format(pd.median_data[0]))
    logger.debug("last line of median data - pd.median_data[-1]:  {}".format(pd.median_data[-1]))

    pd.count_headers = data[count_range[0]+1]
    pd.count_data = data[(count_range[0]+2):(count_range[1]-1)]
    logger.debug("first line of count data - pd.count_data[0]:  {}".format(pd.count_data[0]))
    logger.debug("last line of count data - pd.count_data[-1]:  {}".format(pd.count_data[-1]))

    return pd
