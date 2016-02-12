import setup_logger
import logging
import poolset_data
import csv


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def parse_location_to_well(location):
    split = location.split(",")
    right_paren_index = split[1].index(")")
    return split[1][0:right_paren_index]


def main():
    # output_file = "requirements_artifacts/PCAL002_CS1_ASSEMBLE_MEDIAN.gct"
    # inputs = [("requirements_artifacts/PCAL002_P1_X1.csv", "PS1"),
    #     ("requirements_artifacts/PCAL002_P2_X1.csv", "PS2"),
    #     ("requirements_artifacts/PCAL002_P3_X1.csv", "PS3")]
    # output_file = "requirements_artifacts/PCAL001_CS1_ASSEMBLE_MEDIAN.gct"
    # inputs = [("requirements_artifacts/PCAL001_P1_X1.csv", "PS1"),
    #     ("requirements_artifacts/PCAL001_P2_X1.csv", "PS2"),
    #     ("requirements_artifacts/PCAL001_P3_X1.csv", "PS3")]
    # output_file = "requirements_artifacts/PCAL001_CS1_ASSEMBLE_MEDIAN_from_P7.gct"
    # inputs = [("requirements_artifacts/PCAL001_P7_X1.csv", "PS7")]
    output_file = "requirements_artifacts/PCAL003_CS1_ASSEMBLE_MEDIAN_from_P7.gct"
    inputs = [("requirements_artifacts/PCAL003_P7_X1.csv", "PS7")]

    poolset_data_objects = []
    for (csvfile, ps_id) in inputs:
        pd = poolset_data.read_data(csvfile)
        pd.poolset_id = ps_id
        poolset_data_objects.append(pd)


    f = open("requirements_artifacts/CalicoTranche1PrimaryMetaData_01212016.txt")
    cell_data = f.read().strip().split("\n")
    f.close()
    cell_data = [x.split("\t") for x in cell_data]
    logger.debug("cell_data headers - cell_data[0]:  {}".format(cell_data[0]))
    cd_pool_id_ind = cell_data[0].index("Pool_ID")
    cd_analyte_ind = cell_data[0].index("Analyte")
    cd_name_ind = cell_data[0].index("StrippedName")
    cd_minipool_ind = cell_data[0].index("MiniPool_ID")
    cd_arbitrary_id_ind = cell_data[0].index("Arbitrary_ID")
    cd_ccle_ind = cell_data[0].index("CCLE_Name")
    cd_lua_ind = cell_data[0].index("LUA")

    f = open("requirements_artifacts/pool_to_poolset_mapping.txt")
    raw_pool_to_poolset_mapping = f.read().strip().split("\n")
    f.close()
    raw_pool_to_poolset_mapping = [x.split("\t") for x in raw_pool_to_poolset_mapping]
    logger.debug("raw_pool_to_poolset_mapping[0]:  {}".format(raw_pool_to_poolset_mapping[0]))
    p2pmap = {}
    for row in raw_pool_to_poolset_mapping[1:]:
        pool_id = row[0]
        poolset_id = row[1]
        if poolset_id not in p2pmap:
            p2pmap[poolset_id] = []
        p2pmap[poolset_id].append(pool_id)
    logger.debug("p2pmap:  {}".format(p2pmap))

    prev_wells = None
    all_r = []
    for pd in poolset_data_objects:
        pools = p2pmap[pd.poolset_id]
        logger.debug("pools:  {}".format(pools))

        pools_cd = [x for x in cell_data if x[cd_pool_id_ind] in pools]
        logger.debug("pools_cd:  {}".format(pools_cd))

        analytes = [x[cd_analyte_ind] for x in pools_cd]

        md_analyte_inds = [i for (i,x) in enumerate(pd.median_headers) if x in analytes]
        logger.debug("md_analyte_inds:  {}".format(md_analyte_inds))

        md_to_cd_map = {}
        for i in md_analyte_inds:
            analyte_name = pd.median_headers[i]
            cd_ind = analytes.index(analyte_name)
            md_to_cd_map[i] = cd_ind

        logger.debug("md_to_cd_map:  {}".format(md_to_cd_map))


        r = [[x[cd_name_ind], x[cd_ccle_ind], x[cd_analyte_ind], x[cd_lua_ind],
            x[cd_minipool_ind], x[cd_pool_id_ind], pd.poolset_id, pd.csv_file]
            for x in pools_cd]

        wells = []
        for md in pd.median_data:
            w = parse_location_to_well(md[0])
            wells.append(w)

            for i in md_analyte_inds:
                datum = md[i]
                r_ind = md_to_cd_map[i]
                r[r_ind].append(datum)

        logger.debug("r:  {}".format(r))

        if prev_wells is not None:
            assert prev_wells == wells
        prev_wells = wells

        all_r.extend(r)

    header = [None for i in range(8)]
    header.extend(prev_wells)

    f = open(output_file, "w")

    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    for x in all_r:
        writer.writerow(x)
    f.close()

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    main()
