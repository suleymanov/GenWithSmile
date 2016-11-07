# coding=utf-8
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from contextlib import contextmanager

from eden.graph import Vectorizer
from sklearn import metrics
from scipy.sparse import vstack


gk_def = {
    'p': 0.99999, 'complexity': 4, 'r': None, 'd': None, 
    'min_r': 0,'min_d': 0, 'nbits': 20, 'rjob': 1
}


def is_isomorph(graph1, graph2):
    """
    :param graph1: nx.Graph
    :param graph2: nx.Graph
    return: True, if graphs are isomorphic, otherwise False
    """
    if nx.faster_could_be_isomorphic(graph1, graph2):
        node_match = iso.categorical_node_match('label', 'C')
        edge_match = iso.categorical_node_match(['weight', 'label'], [1, '-'])
        return iso.is_isomorphic(graph1, graph2, node_match=node_match, edge_match=edge_match)
    return False


def is_isomorph_gk(graph1_kernel_vect, graph2_kernel_vect, gk_params=gk_def):
    """
    :param graph1_kernel_vect: np.ndarray
    :param graph2_kernel_vect: np.ndarray
    :return: bool
    """
    k = metrics.pairwise.pairwise_kernels(graph1_kernel_vect, graph2_kernel_vect, metric='cosine')
    return len(np.where(k[:, 0] > gk_params['p'])[0]) >= 1


def vectorize_mol_graphs(mol_graphs_list, gk_params=gk_def):
    vectorizer = Vectorizer(
        complexity=gk_params['complexity'],
        r=gk_params['r'],
        d=gk_params['d'],
        min_r=gk_params['min_r'],
        min_d=gk_params['min_d'],
        nbits=gk_params['nbits']
    )
    return vectorizer.transform(mol_graphs_list)


def rdkitmol2graph(mol):
    """
    Converts molecule to graph.
    :param mol: rdkit.Chem.rdchem.Mol
    :return: nx.Graph
    """
    mol_atoms = mol.GetAtoms()
    rings = mol.GetRingInfo()
    num_rings = len(rings.BondRings())
    charges = np.array([atom.GetFormalCharge() for atom in mol_atoms] + [0] * num_rings)
    atoms = np.array([atom.GetSymbol() for atom in mol_atoms] + ['R'] * num_rings)
    bonds = np.array([[bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, int(bond.GetBondType())]
        for bond in mol.GetBonds()])
    num_vertex = np.amax(bonds[:, [0, 1]])
    gr = np.zeros((num_vertex + num_rings, num_vertex + num_rings), dtype=int)
    bonds[:, [0, 1]] = bonds[:, [0, 1]] - 1
    for i in xrange(len(bonds)):
        gr[bonds[i, 0], bonds[i, 1]] = bonds[i, 2]
        gr[bonds[i, 1], bonds[i, 0]] = bonds[i, 2]
    graph = nx.Graph()
    for i, symbol in enumerate(atoms):
        graph.add_node(i, label=symbol, entity=charges[i])
    edge_type_to_label = {1: '-', 2: '=', 3: '#', 12: '||', 13: 'RING'}
    for i in xrange(len(atoms)):
        for j in xrange(i + 1, len(atoms)):
            edge_type = gr[i, j]
            if edge_type != 0:
                label = edge_type_to_label.get(edge_type, '')
                graph.add_edge(i, j, weight=edge_type, label=label)
    return graph


def get_nonisomorphic_positions(graph, positions):
    """
    :param graph: nx.Graph
    :param positions: list of int
    :return: list of int
    """
    def _get_graph_combinations():
        graph1, graph2 = graph, graph.copy()
        for i in xrange(len(positions) - 1):
            if not is_isomorph_flag[i]:
                with modify_vertex(graph1, positions[i]):
                    for j in xrange(i + 1, len(positions)):
                        if not is_isomorph_flag[j]:
                            with modify_vertex(graph2, positions[j]):
                                yield (i, graph1), (j, graph2)

    is_isomorph_flag = np.zeros(len(positions), dtype=bool)

    for (_, graph1), (j, graph2) in _get_graph_combinations():
        is_isomorph_flag[j] = is_isomorph(graph1, graph2)
    return list(np.where(~is_isomorph_flag)[0])


@contextmanager
def modify_vertex(graph, node_index):
    """
    :param graph: nx.Graph
    :param node_index: int
    :return: nx.Graph
    """
    new_label = 'new_label'
    old_label = graph.node[node_index]['label']
    graph.node[node_index]['label'] = new_label
    yield
    graph.node[node_index]['label'] = old_label
