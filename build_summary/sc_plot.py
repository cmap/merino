import matplotlib.pyplot as plt
import seaborn as sns

def sc_plot(sig, title, outfile, cellspace=500):
    _types = ['ctl_vehicle', 'trt_cp', 'trt_poscon']
    data_df = sig[(sig.pert_type.isin(_types))]
    size_order = ['20 uM', '15 uM', '10 uM', '5 uM', '3.33 uM', '1.11 uM', '0.5556 uM', '-666']
    ax = sns.relplot(data=sig, x='cc_q75', y='ss_ltn2', size='pert_idose', hue='pert_type', col='pert_type', \
                     height=6, palette=['g', 'r', 'b'], sizes=(40, 100), col_order=_types, size_order=size_order)
    ax.set(xlabel='Replicate Correlation (CC_Q75)', ylabel='Num. Sens. Cell Lines (SS_ltn2)')
    plt.xlim(0, 1)
    plt.ylim(0, cellspace)
    plt.xlabel('cc_q75')
    plt.ylabel('signal strength (cells killed)')
    plt.title(title)
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0.85, 0.85, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(outfile)
    plt.clf()
