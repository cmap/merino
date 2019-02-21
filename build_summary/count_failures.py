import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
import glob

def count_analysis(project_folder, dets=None, perts=None):

    if dets is None:
        dets = list(set(os.path.basename(f)[0:20] for f in glob.glob(os.path.join(project_folder, 'assemble', 'P*'))))

    if perts is None:
        perts = list(
            set(os.path.basename(f)[0:6] for f in glob.glob(os.path.join(project_folder, 'assemble', 'P*'))))


    for pert in perts:
        if not os.path.exists(os.path.join(project_folder, 'QC', 'Replicate_QC', pert)):
            os.mkdir(os.path.join(project_folder, 'QC', 'Replicate_QC', pert))
        pert_outfile = os.path.join(project_folder, 'QC', 'Replicate_QC', pert)
        count_paths = glob.glob(os.path.join(project_folder, 'assemble', pert + '*', pert + '*COUNT*.gct'))
        replicate_failures(count_paths, os.path.join(pert_outfile, pert))


    well_failures = []

    for plate in dets:
        print plate
        if not os.path.exists(os.path.join(project_folder, 'QC', plate)):
            os.mkdir(os.path.join(project_folder, 'QC', plate))
            os.mkdir(os.path.join(project_folder, 'QC', plate, 'distributions'))

        base_outfile = os.path.join(project_folder, 'QC', plate)
        count_path = glob.glob(os.path.join(project_folder, 'assemble', plate + '*', plate + '*_COUNT*.gct'))[0]
        count_distributions(count_path, os.path.join(base_outfile, 'distributions'))
        well_failures.append(count_failures(count_path))
        print 'Number of Failed Wells for {} = {}'.format(plate, count_failures(count_path))

    failure_distributions(well_failures, os.path.join(project_folder, 'QC'))

def count_failures(count_path):
    count_gct = pe.parse(count_path)
    median_df = count_gct.data_df.median()
    failures = median_df[median_df < 20].count()
    return failures

def failure_distributions(failures, outfile):

    n, bins, patches = plt.hist(failures, 50, facecolor='brown', alpha=1)

    plt.xlabel('Number of Failed Wells (Bead Count < 20) By Plate, n={}'.format(len(failures)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Failed Wells By Plate')

    plt.savefig(os.path.join(outfile, 'number_of_failures_by_plate.png'))

    plt.clf()

def replicate_failures(gct_list, outfile):

    failures_df = pd.DataFrame()
    for count_path in gct_list:
        name = os.path.basename(count_path)
        print name
        count_gct = pe.parse(count_path)
        median_df = count_gct.data_df.median()
        failures_df[name] = median_df.values

    failures_per_well = failures_df[failures_df < 20].count(axis=1)


    n, bins, patches = plt.hist(failures_per_well, bins=(0, 1, 2, 3, 4), facecolor='red', alpha=1)

    plt.xlabel('Number of Replicates Failed Per Well Wells (Bead Count < 20) By Replicate Set, n={}'.format(len(failures_per_well)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Failed Wells By Replicate Set')
    ax = plt.gca()
    ax.set_xticks([0, 1, 2, 3])
    plt.savefig(os.path.join(outfile + 'replicate_failures.png'))

    plt.clf()

def count_distributions(count_path, outfile):
    ###############################################################################################
    # Making distribution of ALL count values

    count_gct = pe.parse(count_path)

    count_df = count_gct.data_df
    count_data = count_df.unstack()
    count_data.dropna(inplace=True)

    count_bins = np.linspace(0, max(count_data), 50)
    n, bins, patches = plt.hist(count_data, count_bins, facecolor='yellow', alpha=1)

    plt.xlabel('Bead Count Data, n={}'.format(len(count_data)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Bead Count Data')

    plt.savefig(os.path.join(outfile, 'bead_count_dist.png'))

    plt.clf()

    ####################################################################################################
    # Making Dist of Median Count Values per well
    median_well_count = count_df.median()
    median_well_count.dropna(inplace=True)
    median_bins = np.linspace(0, max(median_well_count), 50)
    n, bins, patches = plt.hist(median_well_count, median_bins, facecolor='orange', alpha=1)

    plt.xlabel('Median Bead Count by Well, n={}'.format(len(median_well_count)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Median Count By Well')

    plt.savefig(os.path.join(outfile, 'median_count_dist.png'))

    plt.clf()