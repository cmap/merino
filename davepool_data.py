import logging
import setup_logger
import csv


logger = logging.getLogger(setup_logger.LOGGER_NAME)

datatype_names = ["Median", "Count"]
date_header = "Date"
datatype_header = "DataType:"


class DavepoolData(object):
    '''
    contains the data associated with a davepool - davepool id, csv file information,
    relevant data read from the csv
    '''
    def __init__(self, csv_filepath=None, csv_datetime=None, median_headers=None, median_data=None,
        count_headers=None, count_data=None, davepool_id=None):

        self.csv_filepath = csv_filepath
        self.csv_datetime = csv_datetime
        self.median_headers = median_headers
        self.median_data = median_data
        self.count_headers = count_headers
        self.count_data = count_data
        self.davepool_id = davepool_id


    def __repr__(self):
        return " ".join(["{}:{}".format(k,v) for (k,v) in self.__dict__.items()])

def get_datatype_range(data, datatype_names):
    datatype_indexes = [i for (i,row) in enumerate(data) if len(row) > 0 and row[0] == datatype_header]
    datatypes = [data[i][1] for i in datatype_indexes]
    datatype_indexes.append(len(data)+1)

    r = {}
    for dn in datatype_names:
        if dn in datatypes:
            dn_index = datatypes.index(dn)
            r[dn] = (datatype_indexes[dn_index], datatype_indexes[dn_index+1])
        else:
            r[dn] = (None, None)

    return r


def get_datetime_from_header_rows(header_rows, csv_filepath):
    datetime_row = [x for x in header_rows if len(x) > 0 and x[0] == date_header]
    assert len(datetime_row) == 1, "expected to find only one datetime_row, found datetime_row:  {}  csv_filepath:  {}".format(
        datetime_row, csv_filepath)

    datetime_row = datetime_row[0]

    r = " ".join(datetime_row[1:]) if len(datetime_row) > 1 else None

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

    datetime = get_datetime_from_header_rows(data[0:datatype_ranges[datatype_names[0]][0]], csv_filepath)

    for dn in datatype_names:
        if datatype_ranges[dn][0] is None:
            raise Exception("davepool_data read_data did not find expected DataType dn:  {}".format(dn))

    median_range = datatype_ranges[datatype_names[0]]
    count_range = datatype_ranges[datatype_names[1]]

    pd = DavepoolData(csv_filepath=csv_filepath, csv_datetime=datetime)
    pd.median_headers = data[median_range[0]+1]
    logger.info("len(pd.median_headers):  {}".format(len(pd.median_headers)))
    
    pd.median_data = data[(median_range[0]+2):(median_range[1]-1)]
    logger.info("len(pd.median_data):  {}".format(len(pd.median_data)))
    logger.debug("first line of median data - pd.median_data[0]:  {}".format(pd.median_data[0]))
    logger.debug("last line of median data - pd.median_data[-1]:  {}".format(pd.median_data[-1]))

    pd.count_headers = data[count_range[0]+1]
    logger.info("len(pd.count_headers):  {}".format(len(pd.count_headers)))

    pd.count_data = data[(count_range[0]+2):(count_range[1]-1)]
    logger.info("len(pd.count_data):  {}".format(len(pd.count_data)))
    logger.debug("first line of count data - pd.count_data[0]:  {}".format(pd.count_data[0]))
    logger.debug("last line of count data - pd.count_data[-1]:  {}".format(pd.count_data[-1]))

    return pd
