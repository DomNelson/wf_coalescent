import sys, os
import numpy as np
import pandas as pd
import scipy.sparse
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse


class IBDtoGraph:
    def __init__(self, ibd_file=None, ibd_array=None):
        assert not (ibd_file and ibd_array)
        assert (ibd_file or ibd_array)
        
        self.pairwise_total_ibd_array = None
        self.ninds = None
        self.pos = None
        self.graph_params = None
        
        self.ibd_array = ibd_array
        if ibd_file is not None:
            self.ibd_array = np.load(ibd_file)['ibd_array']
            
        self.ibd_df = pd.DataFrame(self.ibd_array, columns=['ind1', 'ind2', 'start', 'end'])
        self.ibd_df['len'] = self.ibd_df['end'] - self.ibd_df['start']
        
        # Could check inds are all ints first
        self.inds = list(set(self.ibd_df[['ind1', 'ind2']].values.astype(int).ravel()))
        self.ninds = max(self.inds)
        
    def set_pairwise_total_IBD(self, maxrows=None, input_file=None, output_file=None):
        self.pairwise_total_ibd_array = scipy.sparse.lil_matrix((self.ninds + 1, self.ninds + 1))

        if input_file:
            assert output_file is not None
            self.pairwise_total_ibd_array = scipy.sparse.load_npz(input_file).tolil()
            return
        
        for i, row in tqdm(self.ibd_df.iterrows(), total=self.ibd_df.shape[0]):
            ind1, ind2 = sorted([int(row.ind1), int(row.ind2)])
            self.pairwise_total_ibd_array[ind1, ind2] += row.len
            if maxrows is not None and i > maxrows:
                break

        if output_file:
            scipy.sparse.save_npz(output_file, self.pairwise_total_ibd_array.tocoo())
                
    def build_graph(self, k=None, iterations=None, outfile=None):
        G = nx.Graph()
        coo = self.pairwise_total_ibd_array.tocoo()
        max_len = coo.data.max()
        for r, c, d in zip(coo.row, coo.col, coo.data):
            G.add_edge(r, c, weight=d / max_len)
            
        if outfile is not None:
            nx.nx_pydot.write_dot(G, os.path.expanduser(outfile))
            
        self.graph_params = {'k': k, 'iterations': iterations}
        pos = nx.spring_layout(G, k=k, iterations=iterations)
        self.pos = np.vstack(list(pos.values()))
    
    def plot_graph(self, outfile=None, **plot_args):
        if self.pos is None:
            raise ValueError("Must run 'build_graph' before plotting!")
        
        p_args = {'s': 1, 'alpha': 0.5}
        for k, v in plot_args.items():
            p_args[k] = v
            
        plt.scatter(self.pos[:, 0], self.pos[:, 1], **p_args)
        plt.title(str(self.graph_params['iterations']) + ' iterations')
        
        if outfile:
            plt.savefig(os.path.expanduser(outfile))
            
        try:
            plt.show()
        except:
            pass
        

def main(args):
    ibd_file = os.path.expanduser(args.ibd_file)
    plot_file = os.path.expanduser(args.plot_file)

    pairwise_total_ibd_input = None
    if args.pairwise_total_ibd_input:
        pairwise_total_ibd_input = os.path.expanduser(args.pairwise_total_ibd_input)

    pairwise_total_ibd_output = None
    if args.pairwise_total_ibd_output:
        pairwise_total_ibd_output = os.path.expanduser(args.pairwise_total_ibd_output)

    I = IBDtoGraph(ibd_file=ibd_file)
    I.set_pairwise_total_IBD(maxrows=args.max_rows, input_file=pairwise_total_ibd_input,
            output_file=pairwise_total_ibd_output)
    I.build_graph(iterations=args.iterations, k=args.point_spread)
    I.plot_graph(alpha=args.plot_alpha, outfile=plot_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ibd_file", required=True)
    parser.add_argument("-o", "--plot_file", required=True)
    parser.add_argument("-p", "--pairwise_total_ibd_input")
    parser.add_argument("-s", "--pairwise_total_ibd_output")
    parser.add_argument("-a", "--plot_alpha", type=float, default=0.5)
    parser.add_argument("-m", "--max_rows", type=int)
    parser.add_argument("-k", "--point_spread", type=float)
    parser.add_argument("-t", "--iterations", type=int, default=50)

    args = parser.parse_args()
    main(args)
