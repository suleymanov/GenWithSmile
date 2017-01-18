"""
Module with algorithmic utilities.
"""

import numpy as np

import networkx as nx
import networkx.algorithms.isomorphism as iso

from contextlib import contextmanager


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


def get_unique_coords(molecule, num_points=1, positions_to_check=None):
    """
    Get positions which would provide unique modification result.
    :param molecule: modifier.Molecule
    :param num_point: int
    :param use_positions: list of int
    :return: list of list of int
    """
    mol_rdkit = molecule.rdkit
    mol_graph = molecule.graph
    if mol_graph is None or mol_rdkit is None:
        raise ValueError('Cannot get unique indices of None.')
    if mol_rdkit.GetNumAtoms() == 1:
        if num_points == 1:
            return [[0]]
    if num_points == 1:
        positions = map(lambda pos: [pos], list(xrange(mol_rdkit.GetNumAtoms())))
    elif num_points == 2:
        positions = map(
            lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            filter(lambda b: b.IsInRing(), mol_rdkit.GetBonds()))
    else:
        raise ValueError('Three and more point modifications await their implementation.')
    input_positions = (
        filter(lambda coord: all(x in positions_to_check for x in coord), positions)
        if positions_to_check else positions)
    non_iso_inds = get_nonisomorphic_positions(mol_graph, input_positions)
    return map(lambda ind: input_positions[ind], non_iso_inds)


def get_unique_mols(mol_list):
    """
    Get all unique (non-isomorphic) molecules from list.
    :param mol_list: list of modifier.Molecule
    :return: list of modifier.Molecule
    """
    if not mol_list:
        return []
    unique_inds = [0]
    mol_graphs = map(lambda x: x.mol.graph, mol_list)
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
