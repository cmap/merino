import setup_logger
import logging
import ConfigParser
import prism_pipeline

logger = logging.getLogger(setup_logger.LOGGER_NAME)


_config_file_function_column_section = "prism_cell functional column headers"
_config_file_copy_column_section = "prism_cell copy column headers"


class PrismCell(object):
    def __init__(self, pool_id=None, analyte_id=None):
        self.pool_id = pool_id
        self.analyte_id = analyte_id

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


def read_prism_cell_from_file(filepath, config_filepath = prism_pipeline.default_config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    (headers, data) = _read_data(filepath)

    header_map = _generate_header_map(headers, cp)

    return _parse_data(header_map, data)


def _parse_data(header_map, data):
    r = []
    for row in data:
        pc = PrismCell()
        r.append(pc)
        for (h,i) in header_map.items():
            pc.__dict__[h] = row[i]

    return r


def _generate_header_map(headers, config):
    columns = list()
    columns.extend(config.items(_config_file_function_column_section))
    columns.extend(config.items(_config_file_copy_column_section))

    header_map = {}
    for (c_key, c_header_name) in columns:
        if c_header_name in headers:
            header_map[c_key] = headers.index(c_header_name)

    return header_map


def _read_data(tsv_filepath):
    f = open(tsv_filepath)
    raw_data = f.read().strip().split("\n")
    f.close()

    split_raw_data = [x.split("\t") for x in raw_data]

    headers = [x.lower() for x in split_raw_data[0]]
    data = split_raw_data[1:]
    return (headers, data)
