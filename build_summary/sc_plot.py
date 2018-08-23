import matplotlib.pyplot as plt

def sc_plot(sig, title, outfile, cellspace=500):
    pos = sig[sig['pert_type'] == 'trt_poscon']
    ctl = sig[sig['pert_type'] == 'ctl_vehicle']
    trt = sig[sig['pert_type'] == 'trt_cp']
    plt.scatter(trt['cc_q75'], trt['ss_ltn2'], label='trt_cp', color='g')
    plt.scatter(ctl['cc_q75'], ctl['ss_ltn2'], label='ctl_vehicle', color='b')
    plt.scatter(pos['cc_q75'], pos['ss_ltn2'], label='trt_poscon', color='r')
    plt.xlim(0, 1)
    plt.ylim(0, cellspace)
    plt.xlabel('cc_q75')
    plt.ylabel('signal strength (cells killed)')
    plt.title(title)
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0.85, 0.85, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(outfile)
    plt.clf()
