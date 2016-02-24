import logging
import setup_logger
import unittest
import davepool_data

logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestDavepoolData(unittest.TestCase):
    def test___init__(self):
        r = davepool_data.DavepoolData()
        assert hasattr(r, "median_headers")

    def test_read_data(self):
        csv_filepath = "requirements_artifacts/PCAL001_P1_X1.csv"
        r = davepool_data.read_data(csv_filepath)
        assert r is not None
        assert r.csv_filepath == csv_filepath, r.csv_filepath
        assert r.median_headers is not None
        assert r.median_data is not None
        assert r.count_headers is not None
        assert r.count_data is not None

    def test_get_datatype_range(self):
        data = [["a"], range(3), ["DataType:","my type"], range(4), range(5)]

        r = davepool_data.get_datatype_range(data, ["my type"])

        assert r is not None
        logger.debug("r:  {}".format(r))
        assert "my type" in r
        r = r["my type"]
        assert r[0] == 2, r[0]
        assert r[1] == 5, r[1]

        r = davepool_data.get_datatype_range(data, ["unfounded"])
        assert r is not None
        logger.debug("r:  {}".format(r))
        assert "unfounded" in r
        r = r["unfounded"]
        assert r[0] is None, r[0]
        assert r[1] is None, r[1]

        data.append(["DataType:", "another type"])
        data.append(range(6))
        r = davepool_data.get_datatype_range(data, ["my type"])

        assert r is not None
        logger.debug("r:  {}".format(r))
        assert "my type" in r
        r = r["my type"]
        assert r[0] == 2, r[0]
        assert r[1] == 5, r[1]

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
