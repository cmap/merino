import logging
import setup_logger
import ConfigParser

logger = logging.getLogger(setup_logger.LOGGER_NAME)


class AssayPlate(object):
    def __init__(self, assay_plate_barcode=None, pool_id=None):
        self.assay_plate_barcode = assay_plate_barcode
        self.pool_id = pool_id
        self.compound_plate = {}

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


def read_assay_plates(tsv_filepaths, config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    has_problems = False
    r = {}

    for tsv_filepath in tsv_filepaths:
        (headers, data) = _read_data(tsv_filepath)

        header_map = _generate_header_map(headers, cp, "column headers")

        for (i, row) in enumerate(data):
            apb = row[header_map["assay_plate_barcode"]]
            pool_id = row[header_map["pool_id"]]
            compound_plate_map_name = row[header_map["compound_plate_map_name"]]
            poscon_plate_map_name = row[header_map["poscon_plate_map_name"]]
            if apb not in r:
                new_ap = AssayPlate(assay_plate_barcode=apb, pool_id=pool_id)
                new_ap.compound_plate[compound_plate_map_name] = poscon_plate_map_name
                r[apb] = new_ap

            ap = r[apb]

            if ap.pool_id != pool_id:
                logger.error("assay plate with 2 different pool_id - i:  {}  ap.assay_plate_barcode:  {}  ap.pool_id:  {}  "
                                "pool_id:  {}".format(i, ap.assay_plate_barcode, ap.pool_id, pool_id))
                has_problems = True

            if compound_plate_map_name in ap.compound_plate:
                if ap.compound_plate[compound_plate_map_name] != poscon_plate_map_name:
                    logger.error("assay plate with 2 different compound_plate_map_name for the same compound_plate_barcode - "
                                "i:  {}  ap.assay_plate_barcode:  {}  ap.compound_plate:  {}  compound_plate_map_name:  {}  "
                                    "poscon_plate_map_name:  {}".format(i, ap.assay_plate_barcode, ap.compound_plate,
                                                                      compound_plate_map_name, poscon_plate_map_name))
                    has_problems = True
            else:
                ap[compound_plate_map_name] = poscon_plate_map_name

    if has_problems:
        raise Exception("assay_plate read_assay_plates problems encountered, see error messages")

    return r.values()

def _read_data(tsv_filepath):
    f = open(tsv_filepath)
    raw_data = f.read().strip().split("\n")
    f.close()

    split_raw_data = [x.split("\t") for x in raw_data]

    headers = [x.lower() for x in split_raw_data[0]]
    logger.debug("headers:  {}".format(headers))

    data = split_raw_data[1:]
    return (headers, data)


def _generate_header_map(headers, config, config_section):
    columns = config.items(config_section)

    header_map = {}
    for (c_key, c_header_name) in columns:
        if c_header_name in headers:
            header_map[c_key] = headers.index(c_header_name)

    logger.debug("header_map:  {}".format(header_map))
    return header_map

def _parse_data(header_map, data, BuildClass):
    r = []
    for row in data:
        bc = BuildClass()
        r.append(bc)
        for (h,i) in header_map.items():
            if len(row) > i:
                val = row[i]
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

                bc.__dict__[h] = val

    return r
