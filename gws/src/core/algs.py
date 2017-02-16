"""
Module with algorithmic utilities.
"""

from contextlib import contextmanager

import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso


def is_isomorph(graph1, graph2):
    """
    Check if molecular graphs are isomorphic.
    :param graph1: nx.Graph
    :param graph2: nx.Graph
    return: bool
    """
    if nx.faster_could_be_isomorphic(graph1, graph2):
        node_match = iso.categorical_node_match('label', 'C')
        edge_match = iso.categorical_edge_match(['weight', 'label'], [1, '-'])
        return iso.is_isomorphic(graph1, graph2, node_match=node_match, edge_match=edge_match)
    return False


def rdkitmol2graph(mol):
    """
    Converts rdkit molecule to graph.
    :param mol: rdkit.Chem.rdchem.Mol
    :return: nx.Graph
    """
    atoms = np.array(map(lambda x: x.GetSymbol(), mol.GetAtoms()))
    bonds = np.array(map(lambda x: [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), int(x.GetBondType())], 
        mol.GetBonds()))
    num_vertex = len(atoms)
    gr = np.zeros((num_vertex, num_vertex), dtype=int)
    for i in xrange(len(bonds)):
        gr[bonds[i, 0], bonds[i, 1]] = bonds[i, 2]
        gr[bonds[i, 1], bonds[i, 0]] = bonds[i, 2]
    graph = nx.Graph()
    for i, symbol in enumerate(atoms):
        graph.add_node(i, label=symbol)
    edge_type_to_label = {1: '-', 2: '=', 3: '#', 12: '||'}
    for i in xrange(len(atoms)):
        for j in xrange(i + 1, len(atoms)):
            edge_type = gr[i, j]
            if edge_type != 0:
                label = edge_type_to_label.get(edge_type, '')
                graph.add_edge(i, j, weight=edge_type, label=label)
    return graph


def get_unique_coords(mol_graph, num_points, positions):
    """
    Get positions which would provide unique modification result.
    (TODO: don't forget that coordinates should be provided explicitly.)
    :param mol_graph: nx.Graph
    :param num_points: int
    :param positions: list of list of int
    :return: list of list of int
    """
    assert num_points == 1 or num_points == 2, \
        'Three and more point modifications await their implementation.'
    assert all(map(lambda x: len(x) == num_points, positions))

    if mol_graph is None:
        raise ValueError('Cannot get unique indices of None.')
    non_iso_inds = _get_nonisomorphic_positions(mol_graph, positions)
    return map(lambda ind: positions[ind], non_iso_inds)


def get_unique_mols(mol_graphs):
    """
    Get indices of unique (non-isomorphic) molecular graphs from list.
    :param mol_graphs: list of nx.Graph
    :return: list of int
    """
    if len(mol_graphs) == 0:
        return []
    unique_inds = [0]
    if len(mol_graphs) == 1:
        return unique_inds
    for i, mol_graph in enumerate(mol_graphs[1:]):
        if any(is_isomorph(mol_graph, mol_graphs[ind]) for ind in unique_inds):
            continue
        unique_inds.append(i + 1)
    return unique_inds


def conn_index(mol, expon):
    """
    Return molecular graph connectivity index for provided exponent value.
    :param mol: rdkit.Chem.rdchem.Mol
    :param expon: float
    :return: float
    """
    return sum(map(
        lambda b: (b.GetBeginAtom().GetDegree() * b.GetEndAtom().GetDegree()) ** expon,
        mol.GetBonds()))


def _get_nonisomorphic_positions(graph, positions):
    def _get_graph_combinations():
    	graph1, graph2 = graph, graph.copy()
    	for i in xrange(len(positions) - 1):
    		if not is_isomorph_flag[i]:
    			with _modify_vertices(graph1, positions[i]):
    				for j in xrange(i + 1, len(positions)):
    					if not is_isomorph_flag[j]:
    						with _modify_vertices(graph2, positions[j]):
    							yield (i, graph1), (j, graph2)

    is_isomorph_flag = np.zeros(len(positions), dtype=bool)
    for (_, graph1), (j, graph2) in _get_graph_combinations():
    	is_isomorph_flag[j] = is_isomorph(graph1, graph2)
    return list(np.where(~is_isomorph_flag)[0])


@contextmanager
def _modify_vertices(graph, nodes_indices):
	new_labels = ['new_label'] * len(nodes_indices)
	old_labels = map(lambda ind: graph.node[ind]['label'], nodes_indices)
	for i, ind in enumerate(nodes_indices):
		graph.node[ind]['label'] = new_labels[i]
	yield
	for i, ind in enumerate(nodes_indices):
		graph.node[ind]['label'] = old_labels[i]
