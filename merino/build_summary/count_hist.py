import sys
import pandas as pd
import seaborn as sns
import string
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os


def count_dist(count_path, title, output):
    count_gct = pe.parse(count_path)
    well_medians = count_gct.data_df.median(axis=0)
    n,bins,patch = plt.hist(well_medians.dropna(), 50, color='orange')
    plt.title(title)
    plt.savefig(output)
    plt.clf()