"""
Module with algorithmic utilities.
"""

import numpy as np

import networkx as nx
import networkx.algorithms.isomorphism as iso

from contextlib import contextmanager

from sklearn import metrics
# from eden.graph import Vectorizer


gk_def = {
    'p': 0.99999, 'complexity': 4, 'r': None, 'd': None, 
    'min_r': 0,'min_d': 0, 'nbits': 20, 'rjob': 1
}


def is_isomorph(graph1, graph2):
    """
    Check if graphs are isomorphic.
    :param graph1: nx.Graph
    :param graph2: nx.Graph
    return: bool
    """
    if nx.faster_could_be_isomorphic(graph1, graph2):
        node_match = iso.categorical_node_match('label', 'C')
        edge_match = iso.categorical_edge_match(['weight', 'label'], [1, '-'])
        return iso.is_isomorphic(graph1, graph2, node_match=node_match, edge_match=edge_match)
    return False


def is_isomorph_gk(graph1_kernel_vect, graph2_kernel_vect, gk_params=gk_def):
    """
    Check if graphs are isomorphic using graph kernels.
    :param graph1_kernel_vect: np.ndarray
    :param graph2_kernel_vect: np.ndarray
    :return: bool
    """
    k = metrics.pairwise.pairwise_kernels(graph1_kernel_vect, graph2_kernel_vect, metric='cosine')
    return len(np.where(k[:, 0] > gk_params['p'])[0]) >= 1


# def vectorize_mol_graphs(mol_graphs_list, gk_params=gk_def):
#     vectorizer = Vectorizer(
#         complexity=gk_params['complexity'],
#         r=gk_params['r'],
#         d=gk_params['d'],
#         min_r=gk_params['min_r'],
#         min_d=gk_params['min_d'],
#         nbits=gk_params['nbits']
#     )
#     return vectorizer.transform(mol_graphs_list)


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
    bonds = np.array([[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(bond.GetBondType())]
        for bond in mol.GetBonds()])
    num_vertex = mol.GetNumAtoms()
    gr = np.zeros((num_vertex + num_rings, num_vertex + num_rings), dtype=int)
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


def get_unique_coords(molecule, num_points=1):
	"""
	Get positions which would provide unique modification result.
	:param molecule: molecule.Molecule
	:param num_points: int
	:return: list of list of int
	"""
	mol_rdkit = molecule.mol_rdkit
	mol_graph = molecule.mol_graph
	if mol_graph is None or mol_rdkit is None:
		raise ValueError('Cannot get unique indices of None.')
	if mol_rdkit.GetNumAtoms() == 1:
		if num_points == 1:
			return [[0]]
	if num_points == 1:
		positions = map(lambda pos: [pos], list(xrange(len(mol_rdkit.GetAtoms()))))
	elif num_points == 2:
		# positions = map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), mol_rdkit.GetBonds())
		positions = map(
            lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            filter(lambda b: b.IsInRing(), mol_rdkit.GetBonds()))
	else:
		raise ValueError('Three and more point modifications await their implementation.')
	non_iso_inds = get_nonisomorphic_positions(mol_graph, positions)
	return map(lambda ind: positions[ind], non_iso_inds)


def get_unique_mols(mol_list, use_gk=False):
	"""
	Get all unique (non-isomorphic) molecules from list.
	:param mol_list: list of molecule.Molecule
	:param use_gk: bool
	:return: list of molecule.Molecule
	"""
	if not mol_list:
		return []
	unique_inds = [0]
	mol_graphs = map(lambda x: x.mol_graph, mol_list)
	if use_gk:
		pass
		# vectorized_graphs = vectorize_mol_graphs(mol_graphs)
		# for i, vectorized_mol in enumerate(vectorized_graphs[1:]):
		# 	if any(is_isomorph_gk(vectorized_mol, vectorized_graphs[ind]) for ind in unique_inds):
		# 		continue
		# 	unique_inds.append(i + 1)
	else:
		for i, mol_graph in enumerate(mol_graphs[1:]):
			if any(is_isomorph(mol_graph, mol_graphs[ind]) for ind in unique_inds):
				continue
			unique_inds.append(i + 1)
	return map(lambda ind: mol_list[ind], unique_inds)


def get_nonisomorphic_positions(graph, positions):
    """
    :param graph: nx.Graph
    :param positions: list of list of int
    :return: list of int
    """
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
