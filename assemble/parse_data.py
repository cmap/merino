import merino.setup_logger as setup_logger
import logging


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def parse_data(header_map, data, BuildClass):
    r = []
    for row in data:
        bc = BuildClass()
        r.append(bc)
        for (h,i) in header_map.items():
            if len(row) > i:
                raw_value = row[i]
                val = parse_raw_value(raw_value)

                bc.__dict__[h] = val

    return r


def parse_raw_value(raw_value):
    val = raw_value
    if val == "":
        val = None
    else:
        try:
            val = int(val)
        except ValueError:
            try:
                val = float(val)
            except ValueError:
                pass

    return val


def generate_header_map(headers, internal_header_file_header_pairs, do_keep_all):
    reverse_header_map = {}
    for (c_key, c_header_name) in internal_header_file_header_pairs:
        reverse_header_map[c_header_name] = c_key

    header_map = {}
    for (i, h) in enumerate(headers):
        if h in reverse_header_map:
            c_key = reverse_header_map[h]
            header_map[c_key] = i
        elif do_keep_all:
            header_map[h] = i

    if len(header_map) == 0:
        raise Exception("prism_metadata _generate_header_map header_map has no entries, possible mismatch between expected and actual columns (e.g. check config file for mapping")

    return header_map


def read_data(tsv_filepath):
    f = open(tsv_filepath, "U")
    raw_data = f.read().strip().split("\n")
    f.close()

    split_raw_data = [x.split("\t") for x in raw_data]

    headers = [x.lower() for x in split_raw_data[0]]
    logger.debug("headers:  {}".format(headers))

    data = split_raw_data[1:]
    return (headers, data)
