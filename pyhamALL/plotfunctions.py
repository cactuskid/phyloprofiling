
from scipy import sparse
from matplotlib import pylab
from networkx import draw_networkx
def plot_rows(row_mat , labels , name = ''):
    x,y,vals = find(result_mat)
    #%pylab inline
    pylab.rcParams['figure.figsize'] = (10, 20)
    fig = plt.figure()
    ax = fig.add_subplot(111 )
    ## TODO: add colors for different biological events
    ax.scatter(x,y, s = .4)
    ax.set_yticks(range(result_mat.shape[1]))
    ax.set_yticklabels(list(labels))
    plt.savefig(config.working_dir+name + 'rows.png')

def plot_distmat(distmat, labels , name = ''):
    #plot a generic distance matrix
    pylab.rcParams['figure.figsize'] = (10, 20)
    fig = plt.figure()
    ax1= fig.add_subplot(111 )
    ax.imshow(distmat)
    ax.set_yticks(range(result_mat.shape[1]))
    ax.set_yticklabels(list(labels))
    ax.set_xticks(range(result_mat.shape[1]))
    ax.set_xticklabels(list(labels))

    ax2 = fig.add_subplot(112 )
    ax.colorbar()
    plt.savefig(config.working_dir + name+'distmat.png')

def plot_network(nxobject, name= ''):
    #generate a semi decent plot using network labels.abs
    #should be a small plot_network
    draw_networkx(nxobject, pos = nx.spring_layout(nxobject))
    limits=plt.axis('off')
    plt.savefig(config.working_dir + name + 'network.png')
