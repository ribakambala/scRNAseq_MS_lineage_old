# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc

print("hello world")

sc.set_figure_params(dpi=100)  # low dpi (dots per inch) yields small inline figures
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

results_file = '../results/R6875_scRNA_v1_20181106/planaria_extended.h5ad'

