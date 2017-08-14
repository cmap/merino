import sys
import pandas as pd
import seaborn as sns
import string
from broadinstitute_cmap.io.pandasGEXpress import parse as pe
import matplotlib.pyplot as plt
import os


def count_dist(count_path, title, output):
    count_gct = pe.parse(count_path)
    well_medians = count_gct.data_df.median(axis=0)
    n,bins,patch = plt.hist(well_medians.dropna(), 50, color='r')
    plt.title(title)
    plt.savefig(output)
    plt.show()
    plt.clf()