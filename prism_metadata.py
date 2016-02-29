import setup_logger
import logging
import ConfigParser
import prism_pipeline
import json


logger = logging.getLogger(setup_logger.LOGGER_NAME)

_prism_cell_config_file_section = "PrismCell column headers"
_perturbagen_config_file_section = "Perturbagen column headers"


class PrismCell(object):
    def __init__(self, pool_id=None, analyte_id=None, davepool_id=None, id=None):
        self.pool_id = pool_id
        self.analyte_id = analyte_id
        self.davepool_id = davepool_id
        self.id = id

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


class Perturbagen(object):
    def __init__(self):
        self.well_id = None

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


def read_prism_cell_from_file(config_filepath = prism_pipeline.default_config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    filepath = cp.get("PrismCell database file", "prism_cell_database_filepath")
    (headers, data) = _read_data(filepath)

    header_map = _generate_header_map(headers, cp, _prism_cell_config_file_section)

    return _parse_data(header_map, data, PrismCell)


def read_perturbagen_from_file(filepath, config_filepath = prism_pipeline.default_config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    (headers, data) = _read_data(filepath)

    header_map = _generate_header_map(headers, cp, _perturbagen_config_file_section)

    perturbagens = _parse_data(header_map, data, Perturbagen)
    _build_additional_perturbagen_info(cp, perturbagens)

    return perturbagens


def _build_additional_perturbagen_info(config, perturbagens):
    pert_type_mapping = json.loads(config.get("Perturbagen values", "pert_type_mapping"))
    logger.debug("pert_type_mapping:  {}".format(pert_type_mapping))

    pert_dose_unit = config.get("Perturbagen values", "pert_dose_unit_uM")
    pert_type_vehicle = config.get("Perturbagen values", "pert_type_vehicle")

    pert_time = config.get("Perturbagen values", "default_pert_time")
    pert_time_unit = config.get("Perturbagen values", "default_pert_time_unit")

    pert_id_DMSO = config.get("Perturbagen values", "pert_id_DMSO")

    for (i,p) in enumerate(perturbagens):
        #multiply by 1000 for unit conversion from mM to uM, then apply diluation factor
        if p.compound_well_mmoles_per_liter is not None:
            try:
                p.pert_dose = 1000.0 * float(p.compound_well_mmoles_per_liter) / float(p.dilution_factor)
            except ValueError as e:
                msg = "the concentration or dilution factors should be numbers, they are not.  p:  {}".format(p)
                logger.exception(msg)
                raise ValueError("prism_metadata _build_additional_perturbagen_info " + msg)
        else:
            p.pert_dose = None

        p.pert_dose_unit = pert_dose_unit
        p.pert_idose = str(p.pert_dose) + " " + p.pert_dose_unit

        p.pert_id = p.pert_mfc_id[0:13] if p.pert_mfc_id is not None else pert_id_DMSO

        assert hasattr(p, "pert_type"), "pert_type attribute is missing from perturbagen i:  {}  p:  {}".format(i, p)
        if p.pert_type in pert_type_mapping:
            p.pert_type = pert_type_mapping[p.pert_type]
        else:
            p.pert_type = pert_type_vehicle

        p.pert_iname = p.pert_mfc_desc if "pert_mfc_desc" in p.__dict__ else None

        p.pert_time = pert_time
        p.pert_time_unit = pert_time_unit
        p.pert_itime = p.pert_time + " " + p.pert_time_unit


def _parse_data(header_map, data, BuildClass):
    r = []
    for row in data:
        bc = BuildClass()
        r.append(bc)
        for (h,i) in header_map.items():
            if len(row) > i:
                bc.__dict__[h] = row[i] if row[i] != "" else None

    return r


def _generate_header_map(headers, config, config_section):
    columns = config.items(config_section)

    header_map = {}
    for (c_key, c_header_name) in columns:
        if c_header_name in headers:
            header_map[c_key] = headers.index(c_header_name)

    logger.debug("header_map:  {}".format(header_map))
    return header_map


def _read_data(tsv_filepath):
    f = open(tsv_filepath)
    raw_data = f.read().strip().split("\n")
    f.close()

    split_raw_data = [x.split("\t") for x in raw_data]

    headers = [x.lower() for x in split_raw_data[0]]
    logger.debug("headers:  {}".format(headers))

    data = split_raw_data[1:]
    return (headers, data)
