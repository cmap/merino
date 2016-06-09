import setup_logger
import logging
import ConfigParser
import prism_pipeline
import json
import copy


logger = logging.getLogger(setup_logger.LOGGER_NAME)

_prism_cell_config_file_section = "PrismCell column headers"
_perturbagen_config_file_section = "Perturbagen column headers"
_assay_plate_config_file_section = "Assay Plate column headers"

cmap_pert_well = "pert_well"


class PrismCell(object):
    def __init__(self, pool_id=None, analyte_id=None, davepool_id=None, id=None):
        self.pool_id = pool_id
        self.analyte_id = analyte_id
        self.davepool_id = davepool_id
        self.id = id

        self.ignore = False

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


class Perturbagen(object):
    def __init__(self, well_id=None):
        self.well_id = well_id

    def __repr__(self):
        return " ".join(["{}:{}".format(str(k),str(v)) for (k,v) in self.__dict__.items()])

    def __str__(self):
        return self.__repr__()


class AssayPlate(object):
    def __init__(self, assay_plate_barcode=None, det_plate=None, pool_id=None):
        self.assay_plate_barcode = assay_plate_barcode
        self.det_plate = det_plate
        self.pool_id = pool_id

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


def read_perturbagen_from_CM_file(filepath, config_filepath = prism_pipeline.default_config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    (headers, data) = _read_data(filepath)

    header_map = _generate_header_map(headers, cp, _perturbagen_config_file_section)

    perturbagens = _parse_data(header_map, data, Perturbagen)
    _build_additional_perturbagen_info(cp, perturbagens)

    return perturbagens


def read_perturbagen_from_CMap_file(filepath):
    (headers, data) = _read_data(filepath)

    r = []

    for row in data:
        p = Perturbagen()
        r.append(p)

        for (i, h) in enumerate(headers):
            if len(row) > i:
                raw_value = row[i]
                val = _parse_raw_value(raw_value)

                p.__dict__[h] = val

        p.well_id = p.__dict__[cmap_pert_well]

    return r


def read_assay_plate_from_file(filepath, config_filepath = prism_pipeline.default_config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    (headers, data) = _read_data(filepath)

    header_map = _generate_header_map(headers, cp, _assay_plate_config_file_section)

    return _parse_data(header_map, data, AssayPlate)


def _build_additional_perturbagen_info(config, perturbagens):
    pert_type_mapping = json.loads(config.get("Perturbagen values", "pert_type_mapping"))
    logger.debug("pert_type_mapping:  {}".format(pert_type_mapping))

    pert_dose_unit = config.get("Perturbagen values", "pert_dose_unit_uM")
    pert_type_vehicle = config.get("Perturbagen values", "pert_type_vehicle")

    pert_time = config.get("Perturbagen values", "default_pert_time")
    pert_time_unit = config.get("Perturbagen values", "default_pert_time_unit")

    pert_id_DMSO = config.get("Perturbagen values", "pert_id_DMSO")

    for (i,p) in enumerate(perturbagens):
        #multiply by 1000 for unit conversion from mM to uM, then apply dilution factor
        if p.compound_well_mmoles_per_liter is not None:
            try:
                p.pert_dose = 1000.0 * float(p.compound_well_mmoles_per_liter) / float(p.dilution_factor)
                p.pert_dose_unit = pert_dose_unit
                p.pert_idose = ("%.2f" % p.pert_dose) + " " + p.pert_dose_unit
            except ValueError as e:
                msg = "the concentration or dilution factors should be numbers, they are not.  p:  {}".format(p)
                logger.exception(msg)
                raise ValueError("prism_metadata _build_additional_perturbagen_info " + msg)
        else:
            p.pert_dose = None
            p.pert_dose_unit = None
            p.pert_idose = None

        p.pert_id = p.pert_mfc_id[0:13] if p.pert_mfc_id is not None else pert_id_DMSO

        assert hasattr(p, "pert_type"), "pert_type attribute is missing from perturbagen i:  {}  p:  {}".format(i, p)
        if p.pert_type in pert_type_mapping:
            p.pert_type = pert_type_mapping[p.pert_type]
        else:
            p.pert_type = pert_type_vehicle

        p.pert_iname = p.pert_mfc_desc if "pert_mfc_desc" in p.__dict__ and p.pert_mfc_desc is not None else p.pert_id

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
                raw_value = row[i]
                val = _parse_raw_value(raw_value)

                bc.__dict__[h] = val

    return r


def _parse_raw_value(raw_value):
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

def validate_perturbagens(perturbagens):
    well_pert_map = {}
    mismatches = {}
    for p in perturbagens:

        well_id = p.well_id

        if not well_id in well_pert_map:
            well_pert_map[well_id] = p
        else:
            prev_p = well_pert_map[well_id]

            p_comp = (p.pert_id, p.pert_idose)
            prev_p_comp = (prev_p.pert_id, prev_p.pert_idose)
            if p_comp != prev_p_comp:
                logger.debug("mismatch well_id:  {}  p.assay_plate_barcode:  {}  p_comp:  {}  prev_p.assay_plate_barcode:  {}  "
                             "prev_p_comp:  {}".format(well_id, p.assay_plate_barcode, p_comp, prev_p.assay_plate_barcode,
                                                       prev_p_comp))

                if well_id not in mismatches:
                    mismatches[well_id] = []
                mismatches[well_id].append(p.assay_plate_barcode)

    if len(mismatches) > 0:
        mismatch_wells = list(mismatches.keys())
        mismatch_wells.sort()
        msg = "the perturbagens provided contain different compounds in the same wells of different assasy plates - mismatch_wells():  {}".format(mismatch_wells)
        logger.error(msg)
        raise Exception("prism_metadata validate_perturbagens " + msg)
    else:
        return well_pert_map
