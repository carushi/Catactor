#!/usr/bin/env python
"""
__init__.py for Catactor
"""

__author__ = "Kawaguchi RK"
__copyright__ = "Copyright 2021"
__credits__ = ["Kawaguchi RK"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Kawaguchi RK"
__email__ = "rkawaguc@cshl.edu"
__status__ = "Development"


import pandas as pd

def read_marker_table(marker):
    if '.tsv' in marker:
        marker_table = pd.read_csv(marker, sep='\t', )
    elif '.csv' in marker:
        marker_table = pd.read_csv(marker, index_col=0)
        if marker_table.columns[0] == 'Unnamed: 0':
            marker_table = marker_table.set_index(marker_table.column[0])
    else: # including ".txt"
        marker_table = pd.read_csv(marker, header=None)
        marker_table.columns = 'marker'
    return marker_table
