import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_enrichment_score(wtcs, cumsum, title, outfile):
#PLOT_ES Generate enrichment plot
#   PLOT_ES(WTKS, MAXI, RS, R)
# See COMPUTE_WTKS

    N = len(cumsum)
    plt.plot(range(0,N), cumsum, 'g', linewidth= 2)
    plt.plot([0, N], [0, 0], 'k-')
    plt.xlim(1, N)
    plt.ylim(-1, 1)
    #plot(r, 0, 'k+');
    plt.xlabel ('Rank in Ordered Dataset')
    plt.ylabel ('Running Enrichment Score')
    #plot(maxi, wtks, 'ro', 'markerfacecolor', 'r')
    plt.title(title + ' score = {}'.format(wtcs))
    plt.savefig(outfile)
    plt.clf()
