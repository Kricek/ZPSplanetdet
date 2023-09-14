import numpy as np
import pandas as pd
from os.path import dirname, abspath, join


class _PlanetDataGeneral(object):
    """
    Common parts for different planet data files
    """
    def __init__(self):
        sefl._read_file()


class PlanetDataSuzuki18(_PlanetDataGeneral):
    """
    Suzuki et al. 2018 planet data
    """
    def _read_file(self):
        """Read input file"""
        package_path = dirname(dirname(abspath(__file__)))
        file_path = join(package_path, 'data', 'planets_data.csv')

        data_read = pd.read_csv(file_path, header=None)
        data_numpy = data_read.to_numpy()

        self._values_s = data_numpy[1:,2].astype(float)
        self._values_q = data_numpy[1:,1].astype(float) * 1.e-3
        self._weights = data_numpy[1:,3].astype(float)

